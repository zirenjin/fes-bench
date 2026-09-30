"""Train and evaluate T3 representation regressors on frozen temp_extrap data.

The E1 configuration is read through ``training_config_audit.csv``; only the
frame lists are replaced by the frozen split's train rows.  This script never
creates a split.  It materializes exact train/test DeepMD systems, asserts
that their ``(system, phase, T_index)`` keys are disjoint, trains one model per
system/basis/seed, calibrates one constant per system on train rows, and writes
raw provenance plus train-vs-test metrics.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import os
import shutil
import subprocess
from collections import defaultdict
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np

from fes_bench.eval.run import _root


SYSTEM_PHASES = {
    "sio2": ("quartz_beta", "cristobalite_beta", "tridymite_p63mmc", "tridymite_c2221"),
    "hf": ("hcp", "bcc"),
}
SOURCE_PHASES = {
    "sio2": {name: name for name in SYSTEM_PHASES["sio2"]},
    "hf": {"hcp": "Hf_hcp", "bcc": "Hf_bcc"},
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def parse_config(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def parse_rows(path: Path) -> dict[str, list[dict[str, Any]]]:
    split = json.loads(path.read_text(encoding="utf-8"))
    return {key: list(split[key]) for key in ("train", "test")}


def audit_row(root: Path, basis: str, seed: int) -> dict[str, str]:
    path = root / "result/experiments/crossing_reevaluation/training_config_audit.csv"
    with path.open(encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    matches = [r for r in rows if r["basis"] == basis and int(r["seed"]) == seed and r["split"] == "temp_extrap"]
    if len(matches) != 1:
        raise ValueError(f"expected one E1 audit row for {basis}/{seed}, got {len(matches)}")
    return matches[0]


def resolve_template(root: Path, cfg: dict[str, Any], audit: dict[str, str], seed: int) -> Path:
    candidates = [root / audit["config_path"]]
    fallback = str(cfg["runtime_e1_config_fallback"]).replace("{seed}", str(seed))
    candidates.append(Path(fallback))
    for candidate in candidates:
        if candidate.exists():
            if sha256(candidate) != audit["config_sha256"]:
                raise ValueError(f"E1 config hash mismatch: {candidate}")
            return candidate
    raise FileNotFoundError(f"E1 config unavailable; tried {candidates}")


def phase_source(dataset_root: Path, system: str, phase: str) -> Path:
    # The canonical DeepMD roots contain phase directories directly (the
    # system is represented by the selected root, not an extra path segment).
    return dataset_root / phase


def load_phase(source: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    sets = sorted(source.glob("set.*"))
    if not sets:
        raise FileNotFoundError(f"no DeepMD sets under {source}")
    arrays: dict[str, list[np.ndarray]] = {key: [] for key in ("coord", "box", "fparam", "free_energy")}
    for set_path in sets:
        for key in arrays:
            arrays[key].append(np.load(set_path / f"{key}.npy"))
    coord = np.concatenate(arrays["coord"], axis=0)
    box = np.concatenate(arrays["box"], axis=0)
    fparam = np.concatenate(arrays["fparam"], axis=0)
    energy = np.concatenate(arrays["free_energy"], axis=0)
    temperatures = np.asarray(fparam[:, 0], dtype=float)
    return coord, box, fparam, energy, temperatures


def materialize_phase(source: Path, target: Path, indexes: list[int]) -> None:
    coord, box, fparam, energy, _ = load_phase(source)
    target.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source / "type.raw", target / "type.raw")
    set_path = target / "set.000"
    set_path.mkdir(exist_ok=True)
    take = np.asarray(indexes, dtype=int)
    np.save(set_path / "coord.npy", coord[take])
    np.save(set_path / "box.npy", box[take])
    np.save(set_path / "fparam.npy", fparam[take])
    np.save(set_path / "free_energy.npy", energy[take])


def index_map(rows: list[dict[str, Any]], system: str) -> dict[str, list[int]]:
    out: dict[str, list[int]] = defaultdict(list)
    for row in rows:
        if row["system"] == system:
            out[str(row["phase"])].append(int(row["T_index"]))
    return dict(out)


def materialize_dataset(root: Path, dataset_root: Path, out: Path, system: str, rows: list[dict[str, Any]]) -> dict[str, list[int]]:
    selected = index_map(rows, system)
    # Frozen splits may intentionally omit a phase (e.g. C2221 was removed
    # from the active SiO2 benchmark).  Never materialize an empty DeepMD
    # system: the phase list is the split's phase list, not a hard-coded
    # inventory.
    for phase in sorted(selected):
        source_phase = SOURCE_PHASES[system][phase]
        materialize_phase(phase_source(dataset_root, system, source_phase), out / phase, selected.get(phase, []))
    return selected


def configure(template: dict[str, Any], train_dirs: list[Path], *, seed: int, out: Path) -> dict[str, Any]:
    config = copy.deepcopy(template)
    out = out.resolve()
    train_dirs = [p.resolve() for p in train_dirs]
    config["training"]["seed"] = seed
    config["model"]["fitting_net"]["seed"] = template["model"]["fitting_net"].get("seed", 17)
    config["training"]["training_data"] = {"systems": [str(p) for p in train_dirs], "pair_systems": [str(p) for p in train_dirs], "batch_size": template["training"]["training_data"].get("batch_size", 2)}
    config["training"]["validation_data"] = {"systems": [str(p) for p in train_dirs], "pair_systems": [str(p) for p in train_dirs], "batch_size": template["training"].get("validation_data", {}).get("batch_size", 2), "numb_btch": len(train_dirs)}
    config["training"]["stat_file"] = str(out / "fes_stats.hdf5")
    config["training"]["disp_file"] = str(out / "train.out")
    config["training"]["save_ckpt"] = str(out / "model.ckpt")
    return config


def run_train(dp: str, repo: Path, config_path: Path, out: Path, init_model: Path) -> Path:
    env = os.environ.copy()
    env["PYTHONPATH"] = str(repo / "src/lib")
    env["TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD"] = "1"
    log = out / "console.log"
    with log.open("w", encoding="utf-8") as handle:
        # Fine-tuning applies the selected pretrained branch (rather than
        # treating the multi-head checkpoint as a single-task init model).
        # This preserves the Domains_Alloy representation and reinitializes
        # the FES fitting head described by the E1 config.
        command = [dp, "--pt", "train", "config.json", "--finetune", str(init_model), "--model-branch", "Domains_Alloy"]
        result = subprocess.run(command, cwd=out, env=env, stdout=handle, stderr=subprocess.STDOUT, check=False)
    (out / "exit_code.txt").write_text(f"{result.returncode}\n", encoding="utf-8")
    if result.returncode:
        raise RuntimeError(f"training failed ({result.returncode}); see {log}")
    checkpoints = sorted(
        out.glob("model.ckpt-*.pt"),
        key=lambda path: int(path.stem.rsplit("-", 1)[-1]),
    )
    if not checkpoints:
        raise FileNotFoundError(f"no checkpoint produced in {out}")
    return checkpoints[-1]


def run_test(dp: str, checkpoint: Path, dataset: Path, out_prefix: Path) -> np.ndarray:
    checkpoint = checkpoint.resolve()
    dataset = dataset.resolve()
    out_prefix = out_prefix.resolve()
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    command = [dp, "--pt", "test", "-m", str(checkpoint), "-s", str(dataset), "-d", str(out_prefix), "-n", "1000000"]
    env = os.environ.copy(); env["TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD"] = "1"
    result = subprocess.run(command, capture_output=True, text=True, check=False, env=env)
    (out_prefix.parent / f"{out_prefix.name}.log").write_text(result.stdout + "\n" + result.stderr, encoding="utf-8")
    if result.returncode:
        raise RuntimeError(f"evaluation failed ({result.returncode}) for {dataset}")
    candidates = [out_prefix.with_suffix(".e_peratom.out"), Path(str(out_prefix) + ".e_peratom.out")]
    detail = next((p for p in candidates if p.exists()), None)
    if detail is not None:
        return np.loadtxt(detail, comments="#", ndmin=2)[:, :2]
    # FES heads are exposed by DeepMD's ``property`` evaluator.  PT writes
    # one ``.property.out.<frame>`` file per frame; property values are total
    # cell values, so convert them to the per-atom convention used by the
    # reference tables.
    property_files = sorted(
        out_prefix.parent.glob(out_prefix.name + ".property.out.*"),
        key=lambda path: int(path.name.rsplit(".", 1)[-1]),
    )
    if not property_files:
        raise FileNotFoundError(f"missing property/e_peratom detail output for {out_prefix}")
    rows = [np.loadtxt(path, comments="#", ndmin=2)[:, :2] for path in property_files]
    natoms = sum(1 for line in (dataset / "type.raw").read_text(encoding="utf-8").splitlines() if line.strip())
    return np.concatenate(rows, axis=0) / float(natoms)


def pair_metrics(ref: np.ndarray, pred: np.ndarray, temperatures: np.ndarray) -> dict[str, Any]:
    roots_ref = _root(ref, temperatures)
    roots_pred = _root(pred, temperatures)
    degenerate = bool(np.all(np.abs(pred) <= 1e-12))
    return {
        "delta_G_MAE_eV_per_atom": float(np.mean(np.abs(pred - ref))),
        "zero_floor_delta_G_MAE_eV_per_atom": float(np.mean(np.abs(ref))),
        "delta_G_RMSE_eV_per_atom": float(np.sqrt(np.mean((pred - ref) ** 2))),
        "sign_accuracy": float(np.mean(np.sign(pred) == np.sign(ref))),
        "reference_Tc_K": roots_ref,
        "predicted_Tc_K": ["n/a:degenerate_prediction"] if degenerate else roots_pred,
        "Tc_error_K": ["n/a:degenerate_prediction"] * len(roots_ref) if degenerate else [float(a - b) for a, b in zip(roots_pred, roots_ref)],
        "false_crossings": "n/a:degenerate_prediction" if degenerate else max(0, len(roots_pred) - len(roots_ref)),
        "missed_crossings": "n/a:degenerate_prediction" if degenerate else max(0, len(roots_ref) - len(roots_pred)),
        "n_evaluation_points": int(len(ref)),
        "status": "degenerate_prediction" if degenerate else "ok",
    }


def reference_by_phase(root: Path, system: str, phase: str) -> dict[float, float]:
    processed = {"sio2": "sio2", "hf": "hf"}[system]
    canonical = {"hf": {"hcp": "hcp", "bcc": "bcc"}}.get(system, {}).get(phase, phase)
    path = root / "data/processed" / processed / canonical / "reference_G.csv"
    with path.open(encoding="utf-8", newline="") as handle:
        return {float(row["T_K"]): float(row["G_eV_per_atom"]) for row in csv.DictReader(handle)}


def evaluate_region(root: Path, phases: dict[str, dict[str, dict[str, np.ndarray]]]) -> dict[str, Any]:
    systems = sorted(phases)
    pairs: dict[str, Any] = {}
    ranking_values: list[bool] = []
    ranking_by_system: dict[str, list[bool]] = defaultdict(list)
    scalar_errors: list[float] = []
    for system in systems:
        phase_names = list(phases[system])
        temps_by_phase = {p: phases[system][p]["T"] for p in phase_names}
        refs = {p: reference_by_phase(root, system, p) for p in phase_names}
        for left, right in combinations(phase_names, 2):
            left_data, right_data = phases[system][left], phases[system][right]
            common = sorted(set(left_data["T"]) & set(right_data["T"]))
            li = [int(np.where(left_data["T"] == t)[0][0]) for t in common]
            ri = [int(np.where(right_data["T"] == t)[0][0]) for t in common]
            ref = np.asarray([refs[left][float(t)] - refs[right][float(t)] for t in common], dtype=float)
            pred = left_data["P"][li] - right_data["P"][ri]
            pair = pair_metrics(ref, pred, np.asarray(common, dtype=float))
            pairs[f"{system}:{left}_minus_{right}"] = pair
        for p in phase_names:
            scalar_errors.extend(np.abs(phases[system][p]["P"] - phases[system][p]["Y"]).tolist())
        common_all = sorted(set.intersection(*(set(v) for v in temps_by_phase.values())))
        for t in common_all:
            values = [phases[system][p]["P"][np.where(phases[system][p]["T"] == t)[0][0]] for p in phase_names]
            correct = all(np.sign(values[i] - values[j]) == np.sign(refs[phase_names[i]][float(t)] - refs[phase_names[j]][float(t)]) for i, j in combinations(range(len(phase_names)), 2))
            ranking_values.append(correct); ranking_by_system[system].append(correct)
    return {"G_MAE_eV_per_atom": float(np.mean(scalar_errors)) if scalar_errors else None, "n_test_frames": int(sum(len(v["T"]) for ps in phases.values() for v in ps.values())), "pairs": pairs, "ranking_accuracy": float(np.mean(ranking_values)) if ranking_values else None, "ranking_accuracy_by_system": {k: float(np.mean(v)) for k, v in ranking_by_system.items() if v}}


def aggregate_regions(regions: list[dict[str, Any]]) -> dict[str, Any]:
    total = sum(int(item.get("n_test_frames", 0)) for item in regions)
    g = sum(float(item["G_MAE_eV_per_atom"]) * int(item["n_test_frames"]) for item in regions if item.get("G_MAE_eV_per_atom") is not None) / total
    pairs = {key: value for item in regions for key, value in item.get("pairs", {}).items()}
    ranking = sum(float(item["ranking_accuracy"]) * int(item["n_test_frames"]) for item in regions if item.get("ranking_accuracy") is not None) / total
    return {"G_MAE_eV_per_atom": g, "n_test_frames": total, "pairs": pairs, "ranking_accuracy": ranking}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", type=Path, default=Path("."))
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--dataset-root", type=Path, required=True)
    ap.add_argument("--init-model", type=Path, required=True)
    ap.add_argument("--dp", default="dp")
    ap.add_argument("--output-root", type=Path, default=Path("result/experiments/t3_temp_extrap"))
    ap.add_argument("--seeds", nargs="+", type=int)
    ap.add_argument("--systems", nargs="+", choices=sorted(SYSTEM_PHASES))
    args = ap.parse_args(); root = args.repo_root.resolve(); config_path = args.config if args.config.is_absolute() else root / args.config; cfg = parse_config(config_path); split_path = root / cfg["split"]; split_rows = parse_rows(split_path)
    seeds = args.seeds or [int(seed) for seed in cfg["seeds"]]
    systems = args.systems or list(cfg["systems"])
    train_keys = {(r["system"], r["phase"], int(r["T_index"])) for r in split_rows["train"]}; test_keys = {(r["system"], r["phase"], int(r["T_index"])) for r in split_rows["test"]}
    overlap = train_keys & test_keys
    if overlap:
        raise AssertionError(f"frozen train/test overlap: {sorted(overlap)[:5]}")
    split_hash = sha256(split_path); cfg_hash = sha256(config_path); args.output_root.mkdir(parents=True, exist_ok=True)
    all_system_metrics: dict[str, dict[str, Any]] = {}
    for seed in seeds:
        audit = audit_row(root, cfg["basis"], int(seed)); template_path = resolve_template(root, cfg, audit, int(seed)); template = parse_config(template_path)
        for system in systems:
            run = args.output_root / "runs" / cfg["basis"] / system / f"seed_{seed}"; run.mkdir(parents=True, exist_ok=True)
            train_data = run / "data" / "train"; test_data = run / "data" / "test"
            train_selected = materialize_dataset(root, args.dataset_root, train_data, system, split_rows["train"])
            test_selected = materialize_dataset(root, args.dataset_root, test_data, system, split_rows["test"])
            for phase, indexes in train_selected.items():
                if not set(indexes) <= {k[2] for k in train_keys if k[0] == system and k[1] == phase}:
                    raise AssertionError(f"train materialization mismatch for {system}/{phase}")
            active_phases = sorted(train_selected)
            model_cfg = configure(template, [train_data / p for p in active_phases], seed=int(seed), out=run)
            model_cfg["model"]["type_map"] = ["Si", "O"] if system == "sio2" else ["Hf"]
            (run / "config.json").write_text(json.dumps(model_cfg, indent=2) + "\n", encoding="utf-8")
            provenance = {"basis": cfg["basis"], "seed": seed, "system": system, "split": "temp_extrap", "split_sha256": split_hash, "predictor_config": str(config_path.relative_to(root)), "predictor_config_sha256": cfg_hash, "e1_config_path": audit["config_path"], "e1_config_sha256": audit["config_sha256"], "e1_template_runtime_path": str(template_path), "train_frame_count": sum(len(v) for v in train_selected.values()), "test_frame_count": sum(len(v) for v in test_selected.values()), "train_test_overlap": len(overlap), "checkpoint_path": None, "checkpoint_sha256": None, "timestamp_utc": datetime.now(timezone.utc).isoformat()}
            (run / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
            checkpoint = run_train(args.dp, root, run / "config.json", run, args.init_model)
            checkpoint = checkpoint.resolve()
            provenance["checkpoint_path"] = str(checkpoint)
            provenance["checkpoint_sha256"] = sha256(checkpoint)
            (run / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
            pred_by_region: dict[str, dict[str, dict[str, np.ndarray]]] = {"train": {system: {}}, "test": {system: {}}}
            for region, data_dir, selected in (("train", train_data, train_selected), ("test", test_data, test_selected)):
                for phase in active_phases:
                    detail = run_test(args.dp, checkpoint, data_dir / phase, run / f"{region}_{phase}")
                    temperatures = np.load(data_dir / phase / "set.000" / "fparam.npy")[:, 0]
                    pred_by_region[region][system][phase] = {"T": temperatures, "Y": detail[:, 0], "P": detail[:, 1]}
            # Fit one scalar calibration per system from train rows only and
            # apply that fixed offset to both regions.  Pair differences are
            # invariant to this offset, while scalar G diagnostics are not.
            train_values = [entry["Y"] - entry["P"] for entry in pred_by_region["train"][system].values()]
            calibration = float(np.concatenate(train_values).mean()) if train_values else 0.0
            for region in ("train", "test"):
                for entry in pred_by_region[region][system].values():
                    entry["P"] = entry["P"] + calibration
            all_system_metrics.setdefault(str(seed), {})[system] = {"run": str(run), "checkpoint": str(checkpoint), "pred": pred_by_region, "train_selected": train_selected, "test_selected": test_selected}
    # Persist the compact handoff and normalized raw runs.  The aggregation
    # below is purely post-processing and never touches checkpoints.
    all_seed_rows: list[dict[str, Any]] = []
    for seed, systems_data in all_system_metrics.items():
        train_eval = evaluate_region(root, {system: data["pred"]["train"][system] for system, data in systems_data.items()})
        test_eval = evaluate_region(root, {system: data["pred"]["test"][system] for system, data in systems_data.items()})
        seed_row = {"basis": cfg["basis"], "seed": int(seed), "train": train_eval, "test": test_eval}
        all_seed_rows.append(seed_row)
        raw_dir = args.output_root / "raw_runs" / f"repr_regression_{cfg['basis']}" / f"seed_{seed}"
        raw_dir.mkdir(parents=True, exist_ok=True)
        try:
            commit = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
        except (OSError, subprocess.CalledProcessError):
            commit = "unknown"
        provenance = {"git_commit": commit, "split": "temp_extrap", "split_sha256": split_hash, "predictor_config": str(config_path.relative_to(root)), "predictor_config_sha256": cfg_hash, "checkpoint_path": [data["checkpoint"] for data in systems_data.values()], "checkpoint_sha256": [sha256(Path(data["checkpoint"])) for data in systems_data.values()], "evaluator_version": "train_t3_temp_extrap.py", "timestamp_utc": datetime.now(timezone.utc).isoformat(), "calibration": "one constant per system fitted on train rows"}
        payload = {"provenance": provenance, "metrics": {"method": f"repr_regression_{cfg['basis']}", "folds": {"all": test_eval}, "train_region": train_eval}}
        (raw_dir / "metrics.json").write_text(json.dumps(payload, indent=2, default=lambda x: x.tolist() if isinstance(x, np.ndarray) else x) + "\n", encoding="utf-8")
    output_rows = []
    for row in all_seed_rows:
        for region_name in ("train", "test"):
            region = row[region_name]
            numeric_pairs = [p for p in region["pairs"].values() if isinstance(p.get("delta_G_MAE_eV_per_atom"), (int, float))]
            pair_mae = float(np.mean([p["delta_G_MAE_eV_per_atom"] for p in numeric_pairs])) if numeric_pairs else None
            zero_mae = float(np.mean([p["zero_floor_delta_G_MAE_eV_per_atom"] for p in numeric_pairs])) if numeric_pairs else None
            sign = float(np.mean([p["sign_accuracy"] for p in numeric_pairs])) if numeric_pairs else None
            skill = None if pair_mae is None or not zero_mae else 1.0 - pair_mae / zero_mae
            output_rows.append({"basis": cfg["basis"], "seed": row["seed"], "region": region_name, "G_MAE_eV_per_atom": region["G_MAE_eV_per_atom"], "delta_G_MAE_eV_per_atom": pair_mae, "skill_score": skill, "sign_accuracy": sign, "ranking_accuracy": region.get("ranking_accuracy")})
    with (args.output_root / "id_vs_ood.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["basis", "seed", "region", "G_MAE_eV_per_atom", "delta_G_MAE_eV_per_atom", "skill_score", "sign_accuracy", "ranking_accuracy"]); writer.writeheader(); writer.writerows(output_rows)
    (args.output_root / "training_index.json").write_text(json.dumps({"config": str(args.config), "split_sha256": split_hash, "runs": all_system_metrics}, indent=2, default=lambda x: x.tolist() if isinstance(x, np.ndarray) else x) + "\n", encoding="utf-8")
    print(args.output_root / "training_index.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
