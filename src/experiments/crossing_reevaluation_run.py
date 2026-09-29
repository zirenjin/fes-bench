"""Produce E1 checkpoint re-evaluation metrics with the unified evaluator — plan experiment E1.

The script is intentionally an adapter, not a trainer.  It evaluates frozen
TorchScript heads on their original static structures, writes temperature-keyed
curves, and then invokes ``fes_bench.eval.run`` with its ``precomputed:``
predictor mode.  Thus the benchmark evaluator, rather than this adapter,
owns root finding, slope conversion, and false/missed-crossing accounting.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


CANONICAL_PHASES = ("quartz_beta", "cristobalite_beta", "tridymite_p63mmc")
LEGACY_PHASES = CANONICAL_PHASES + ("tridymite_c2221",)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def repository_commit(repo: Path) -> str:
    """Return the source revision, including for a detached staging copy."""

    if supplied := os.environ.get("FES_BENCH_GIT_COMMIT"):
        return supplied
    try:
        return subprocess.check_output(
            ["git", "-C", str(repo), "rev-parse", "HEAD"], text=True, stderr=subprocess.DEVNULL
        ).strip()
    except subprocess.CalledProcessError:
        return "unknown"


def external_path(path: Path, root: Path) -> str:
    """Retain a stable external placeholder without publishing host paths."""

    return "external/checkpoints/" + path.relative_to(root).as_posix()


def read_csv(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return (
        np.array([float(row["T_K"]) for row in rows]),
        np.array([float(row["G_eV_per_atom"]) for row in rows]),
    )


def raw_phase(legacy_root: Path, phase: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    phase_root = legacy_root / phase
    atype = np.loadtxt(phase_root / "type.raw", dtype=np.int64)
    parts = []
    for directory in sorted(phase_root.glob("set.*")):
        parts.append(
            (
                np.load(directory / "fparam.npy"),
                np.load(directory / "coord.npy"),
                np.load(directory / "box.npy"),
                np.load(directory / "free_energy.npy"),
            )
        )
    fparam, coord, box, free_energy = (np.concatenate([part[i] for part in parts]) for i in range(4))
    order = np.argsort(fparam[:, 0])
    return fparam[order], coord[order], box[order], free_energy[order], atype, len(atype)


def reference_rows(data_root: Path, phase: str) -> tuple[np.ndarray, np.ndarray]:
    return read_csv(data_root / "sio2" / phase / "reference_G.csv")


def evaluate_model(model_path: Path, legacy_root: Path, device: str, batch_size: int) -> dict[str, dict[str, list[float]]]:
    import torch

    torch.serialization.add_safe_globals([slice])
    from deepmd.infer import DeepPot

    wrapper = DeepPot(str(model_path))
    model = wrapper.deep_eval.dp.model["Default"].eval()
    output: dict[str, dict[str, list[float]]] = {}
    for phase in CANONICAL_PHASES:
        fparam, coord, box, _labels, atype, natoms = raw_phase(legacy_root, phase)
        predictions = []
        for start in range(0, len(fparam), batch_size):
            stop = min(start + batch_size, len(fparam))
            coord_t = torch.tensor(coord[start:stop].reshape(stop - start, natoms, 3), dtype=torch.float64, device=device)
            box_t = torch.tensor(box[start:stop].reshape(stop - start, 3, 3), dtype=torch.float64, device=device)
            atype_t = torch.tensor(np.tile(atype, (stop - start, 1)), dtype=torch.long, device=device)
            fparam_t = torch.tensor(fparam[start:stop], dtype=torch.float64, device=device)
            with torch.no_grad():
                values = model(coord_t, atype_t, box_t, fparam_t)["free_energy"]
            predictions.extend((values.detach().cpu().numpy().reshape(-1) / natoms).tolist())
        output[f"sio2:{phase}"] = {"T_K": fparam[:, 0].astype(float).tolist(), "G_eV_per_atom": predictions}
    return output


def checkpoint_type_map(model_path: Path) -> list[str]:
    """Read the element mapping embedded in the frozen TorchScript model."""

    import torch

    torch.serialization.add_safe_globals([slice])
    from deepmd.infer import DeepPot

    return list(DeepPot(str(model_path)).get_type_map())


def rows_for_temperatures(data_root: Path, temperatures: set[float]) -> list[dict[str, object]]:
    rows = []
    for phase in CANONICAL_PHASES:
        source_t, _ = reference_rows(data_root, phase)
        for index, temperature in enumerate(source_t):
            if float(temperature) in temperatures:
                rows.append({"system": "sio2", "phase": phase, "T_index": index})
    return rows


def calibration(predictions: dict[str, dict[str, list[float]]], data_root: Path, temperatures: set[float]) -> float:
    residuals: list[float] = []
    for phase in CANONICAL_PHASES:
        reference_t, reference_g = reference_rows(data_root, phase)
        predicted = predictions[f"sio2:{phase}"]
        by_t = dict(zip(predicted["T_K"], predicted["G_eV_per_atom"]))
        for temperature, value in zip(reference_t, reference_g):
            if float(temperature) in temperatures:
                residuals.append(float(value) - float(by_t[float(temperature)]))
    if not residuals:
        raise ValueError("calibration window has no shared frames")
    return float(np.mean(residuals))


def shifted(predictions: dict[str, dict[str, list[float]]], constant: float) -> dict[str, dict[str, list[float]]]:
    return {
        key: {"T_K": curve["T_K"], "G_eV_per_atom": [float(value) + constant for value in curve["G_eV_per_atom"]]}
        for key, curve in predictions.items()
    }


def invoke_eval(repo: Path, curves: Path, split: Path, data_root: Path, output_root: Path, seed: int) -> dict[str, object]:
    spec = f"precomputed:{curves}"
    command = [
        sys.executable,
        "-m",
        "fes_bench.eval.run",
        "--predictor",
        spec,
        "--split",
        str(split),
        "--data-root",
        str(data_root),
        "--seeds",
        str(seed),
        "--output-root",
        str(output_root),
    ]
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(repo / "src" / "lib") + os.pathsep + environment.get("PYTHONPATH", "")
    subprocess.run(command, cwd=repo, env=environment, check=True, capture_output=True, text=True)
    slug = spec.replace(":", "_").replace("/", "_")
    metrics_path = output_root / slug / split.stem / "metrics.json"
    metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
    metrics["predictor"] = "precomputed:external/curves/" + curves.name
    return metrics


def sanitize_predictor_paths(payload: dict[str, Any]) -> None:
    """Ensure a promoted E1 result never records a host-specific curve path."""

    for record in payload.get("checkpoint_records", []):
        for window in record.get("metrics", {}).values():
            for metrics in window.values():
                predictor = metrics.get("predictor")
                if isinstance(predictor, str) and predictor.startswith("precomputed:"):
                    metrics["predictor"] = "precomputed:external/curves/" + Path(predictor.removeprefix("precomputed:")).name


def frozen_temperatures(data_root: Path, split_path: Path) -> tuple[set[float], set[float]]:
    """Read calibration regions from, but never modify, the frozen split."""

    split = json.loads(split_path.read_text(encoding="utf-8"))
    train, test = split.get("train"), split.get("test")
    if not isinstance(train, list) or not isinstance(test, list):
        raise ValueError(f"{split_path} must contain train and test rows")
    regions: dict[str, dict[str, set[int]]] = {phase: {"train": set(), "test": set()} for phase in CANONICAL_PHASES}
    for label, rows in (("train", train), ("test", test)):
        for row in rows:
            if not isinstance(row, dict) or row.get("system") != "sio2" or row.get("phase") not in regions:
                continue
            regions[str(row["phase"])][label].add(int(row["T_index"]))
    reference_t, _ = reference_rows(data_root, CANONICAL_PHASES[0])
    expected = {phase: regions[phase] for phase in CANONICAL_PHASES}
    if any(not item["train"] or not item["test"] for item in expected.values()):
        raise ValueError("frozen split lacks train or test frames for an E1 canonical phase")
    if len({tuple(sorted(item["train"])) for item in expected.values()}) != 1 or len({tuple(sorted(item["test"])) for item in expected.values()}) != 1:
        raise ValueError("frozen split does not give the same temperature indices to all E1 phases")
    return (
        {float(reference_t[index]) for index in expected[CANONICAL_PHASES[0]]["train"]},
        {float(reference_t[index]) for index in expected[CANONICAL_PHASES[0]]["test"]},
    )


def comparable_record(record: dict[str, Any]) -> dict[str, Any]:
    """Remove host-specific fields before checking numerical E1 reproducibility."""

    result = dict(record)
    result.pop("checkpoint", None)
    result.pop("config_training_systems", None)
    result.pop("provenance", None)
    normalized_metrics: dict[str, Any] = {}
    for label, window in result.get("metrics", {}).items():
        if not isinstance(window, dict):
            continue
        normalized_windows: dict[str, Any] = {}
        for scope, metrics in window.items():
            if isinstance(metrics, dict):
                normalized = dict(metrics)
                normalized.pop("predictor", None)
                normalized_windows[scope] = normalized
        normalized_metrics[label] = normalized_windows
    result["metrics"] = normalized_metrics
    return result


def scalar_fields(value: Any, path: str = "") -> list[tuple[str, Any]]:
    """Flatten every comparable scalar, preserving list cardinality as a field."""

    if isinstance(value, dict):
        return [item for key in sorted(value) for item in scalar_fields(value[key], f"{path}/{key}")]
    if isinstance(value, list):
        rows = [(f"{path}/length", len(value))]
        for index, item in enumerate(value):
            rows.extend(scalar_fields(item, f"{path}/{index}"))
        return rows
    return [(path, value)]


def tolerance_for(field: str, tolerance_config: dict[str, Any]) -> tuple[str, float]:
    """Classify a scalar by the fixed, resolution-derived E1 tolerance policy."""

    rules = tolerance_config["absolute_tolerances"]
    components = field.split("/")
    if any(part in {"false_crossings", "missed_crossings", "length"} for part in components):
        return "discrete", float(rules["discrete"])
    if any(
        "sign_accuracy" in part or part.endswith("coverage") or "fraction" in part
        or part == "ratio" or part.endswith("_ratio")
        for part in components
    ):
        return "ratio", float(rules["ratio"])
    if "Tc" in field:
        return "temperature_K", float(rules["temperature_K"])
    if "eV_per_atom" in field:
        return "energy_eV_per_atom", float(rules["energy_eV_per_atom"])
    return "exact_other", float(rules["exact_other"])


def csv_value(value: Any) -> str:
    return json.dumps(value, sort_keys=True) if isinstance(value, (dict, list)) else str(value)


def compare_records(
    new_payload: dict[str, Any], reference_path: Path, tolerance_path: Path, output_root: Path
) -> dict[str, Any]:
    """Write a fixed-tolerance, field-level audit against the historical E1 JSON."""

    reference = json.loads(reference_path.read_text(encoding="utf-8"))
    tolerance_config = json.loads(tolerance_path.read_text(encoding="utf-8"))
    expected = {(item["basis"], int(item["seed"])): item for item in reference["checkpoint_records"]}
    observed = {(item["basis"], int(item["seed"])): item for item in new_payload["checkpoint_records"]}
    rows: list[dict[str, Any]] = []
    for key in sorted(expected.keys() | observed.keys()):
        historical = comparable_record(expected[key]) if key in expected else None
        rerun = comparable_record(observed[key]) if key in observed else None
        old = dict(scalar_fields(historical)) if historical is not None else {}
        new = dict(scalar_fields(rerun)) if rerun is not None else {}
        for field in sorted(old.keys() | new.keys()):
            category, tolerance = tolerance_for(field, tolerance_config)
            old_value, new_value = old.get(field), new.get(field)
            numeric = isinstance(old_value, (int, float)) and not isinstance(old_value, bool) and isinstance(new_value, (int, float)) and not isinstance(new_value, bool)
            difference = abs(float(old_value) - float(new_value)) if numeric else None
            passed = difference <= tolerance if numeric else old_value == new_value
            rows.append({
                "basis": key[0],
                "seed": key[1],
                "field": field.lstrip("/"),
                "category": category,
                "historical_value": csv_value(old_value),
                "rerun_value": csv_value(new_value),
                "absolute_difference": "" if difference is None else difference,
                "absolute_tolerance": tolerance,
                "status": "pass" if passed else "fail",
            })
    check_path = output_root / "reproduction_check.csv"
    with check_path.open("w", newline="", encoding="utf-8") as handle:
        fields = ["basis", "seed", "field", "category", "historical_value", "rerun_value", "absolute_difference", "absolute_tolerance", "status"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    audit = {
        "comparison": "checkpoint_record_field_tolerance",
        "reference": "legacy_support/phase5/sio2_reanalysis.json",
        "historical_json_sha256": sha256(reference_path),
        "candidate_json_sha256_before_annotation": hashlib.sha256(json.dumps(new_payload, indent=2).encode("utf-8")).hexdigest(),
        "tolerance_config": "configs/reproduction_tolerance.yaml",
        "tolerance_config_sha256": sha256(tolerance_path),
        "expected_records": len(expected),
        "observed_records": len(observed),
        "checked_fields": len(rows),
        "passed_fields": sum(row["status"] == "pass" for row in rows),
        "failed_fields": sum(row["status"] == "fail" for row in rows),
        "status": "pass" if rows and all(row["status"] == "pass" for row in rows) else "changed",
        "ignored_host_specific_fields": ["checkpoint", "config_training_systems", "provenance", "metrics.*.*.predictor"],
        "check_csv": "reproduction_check.csv",
    }
    (output_root / "reproduction_compare.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    return audit


def write_canonical_meta(candidate_path: Path, reference_path: Path, tolerance_path: Path, output_root: Path, audit: dict[str, Any]) -> None:
    """Write non-self-referential provenance for the promoted canonical E1 JSON."""

    meta = {
        "canonical_json": "result/experiments/crossing_reevaluation/sio2_reanalysis.json",
        "canonical_json_sha256": sha256(candidate_path),
        "historical_json": "result/experiments/legacy_support/phase5/sio2_reanalysis.json",
        "historical_json_sha256": sha256(reference_path),
        "tolerance_config": "configs/reproduction_tolerance.yaml",
        "tolerance_config_sha256": sha256(tolerance_path),
        "reproduction_check": "result/experiments/crossing_reevaluation/reproduction_check.csv",
        "reproduction_check_sha256": sha256(output_root / "reproduction_check.csv"),
        "reproduction_status": audit["status"],
        "checked_fields": audit["checked_fields"],
        "failed_fields": audit["failed_fields"],
    }
    (output_root / "sio2_reanalysis.meta.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint-root", type=Path)
    parser.add_argument("--legacy-data-root", type=Path)
    parser.add_argument("--data-root", type=Path)
    parser.add_argument("--split", type=Path, help="frozen temp_extrap split; read only")
    parser.add_argument("--predictor-config", type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--compare-to", type=Path, help="historical sio2_reanalysis.json used for equality audit")
    parser.add_argument("--tolerance-config", type=Path, help="fixed E1 reproduction tolerance configuration")
    parser.add_argument("--candidate-json", type=Path, help="check an already-produced E1 JSON without evaluating checkpoints")
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--batch-size", type=int, default=16)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[2]
    args.output_root = args.output_root.resolve()
    args.output_root.mkdir(parents=True, exist_ok=True)
    if args.candidate_json is not None:
        if args.compare_to is None or args.tolerance_config is None:
            parser.error("--candidate-json requires --compare-to and --tolerance-config")
        payload = json.loads(args.candidate_json.resolve().read_text(encoding="utf-8"))
        sanitize_predictor_paths(payload)
        audit = compare_records(payload, args.compare_to.resolve(), args.tolerance_config.resolve(), args.output_root)
        payload.setdefault("provenance", {})["reproduction_check"] = audit
        args.candidate_json.resolve().write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        write_canonical_meta(args.candidate_json.resolve(), args.compare_to.resolve(), args.tolerance_config.resolve(), args.output_root, audit)
        if audit["status"] != "pass":
            raise RuntimeError(f"E1 re-evaluation exceeded fixed tolerances: {audit}")
        print(args.candidate_json.resolve())
        return 0
    missing = [name for name in ("checkpoint_root", "legacy_data_root", "data_root", "split", "predictor_config") if getattr(args, name) is None]
    if missing:
        parser.error("checkpoint evaluation requires " + ", ".join("--" + item.replace("_", "-") for item in missing))
    args.checkpoint_root = args.checkpoint_root.resolve()
    args.legacy_data_root = args.legacy_data_root.resolve()
    args.data_root = args.data_root.resolve()
    args.split = args.split.resolve()
    args.predictor_config = args.predictor_config.resolve()
    predictor_config = json.loads(args.predictor_config.read_text(encoding="utf-8"))
    if tuple(predictor_config.get("phases", [])) != CANONICAL_PHASES or tuple(predictor_config.get("bases", [])) != ("polynomial", "tlog_polynomial"):
        raise ValueError("predictor config does not declare the canonical E1 phases and bases")
    if predictor_config.get("frozen_split") != "data/processed/splits/temp_extrap.json":
        raise ValueError("E1 predictor config must name the frozen temp_extrap split")

    protocol = json.loads((args.legacy_data_root / "static_protocol.json").read_text(encoding="utf-8"))
    legacy_input_type_map = (args.legacy_data_root / "type_map.raw").read_text(encoding="utf-8").split()
    first_checkpoint = next((args.checkpoint_root / "polynomial").glob("seed*/model.pth"))
    type_map = checkpoint_type_map(first_checkpoint)
    if type_map != legacy_input_type_map:
        raise ValueError(f"checkpoint type_map {type_map!r} disagrees with legacy inputs {legacy_input_type_map!r}")
    map_rows = []
    for index, legacy_phase in enumerate(LEGACY_PHASES):
        fparam, _coord, _box, labels, _atype, natoms = raw_phase(args.legacy_data_root, legacy_phase)
        item: dict[str, object] = {
            "checkpoint_phase_index": index,
            "legacy_phase": legacy_phase,
            "canonical_phase": legacy_phase if legacy_phase in CANONICAL_PHASES else None,
            "n_atoms": natoms,
            "n_frames": len(fparam),
            "temperature_range_K": [float(fparam[0, 0]), float(fparam[-1, 0])],
            "type_map": type_map,
            "legacy_protocol_n_atoms": protocol["phases"][legacy_phase]["n_atoms"],
        }
        if legacy_phase in CANONICAL_PHASES:
            canonical_t, canonical_g = reference_rows(args.data_root, legacy_phase)
            raw_per_atom = labels.reshape(-1) / natoms
            raw_by_t = dict(zip(fparam[:, 0].astype(float), raw_per_atom))
            canonical_by_t = dict(zip(canonical_t.astype(float), canonical_g))
            common = sorted(set(raw_by_t) & set(canonical_by_t))
            item["reference_G_validation"] = {
                "n_common_temperatures": len(common),
                "max_abs_difference_eV_per_atom": float(max(abs(raw_by_t[t] - canonical_by_t[t]) for t in common)),
                "raw_reference_range_eV_per_atom": [float(min(raw_by_t.values())), float(max(raw_by_t.values()))],
                "canonical_reference_range_eV_per_atom": [float(min(canonical_by_t.values())), float(max(canonical_by_t.values()))],
            }
        map_rows.append(item)
    phase_map = {
        "scope": "historic four-phase checkpoints evaluated only on the named canonical three-phase subset",
        "checkpoint_type_map": type_map,
        "legacy_input_type_map": legacy_input_type_map,
        "checkpoint_phase_order_source": "training_data.systems order in adjacent config.json; model type_map contains species, not phase IDs",
        "entries": map_rows,
        "excluded_from_reporting": "tridymite_c2221",
    }
    (args.output_root / "sio2_checkpoint_phase_map.json").write_text(json.dumps(phase_map, indent=2) + "\n", encoding="utf-8")

    canonical_temperatures, _ = reference_rows(args.data_root, CANONICAL_PHASES[0])
    for phase in CANONICAL_PHASES[1:]:
        phase_temperatures, _ = reference_rows(args.data_root, phase)
        if not np.array_equal(phase_temperatures, canonical_temperatures):
            raise ValueError("canonical three-phase tables do not share a temperature grid")
    train_temperatures, held_temperatures = frozen_temperatures(args.data_root, args.split)
    all_temperatures = train_temperatures | held_temperatures
    if all_temperatures != {float(value) for value in canonical_temperatures}:
        raise ValueError("frozen split does not cover the full shared canonical temperature grid")
    full_split = {"name": "sio2_three_phase_full_range", "train": [], "test": rows_for_temperatures(args.data_root, all_temperatures)}
    subset_split = {
        "name": "sio2_temp_extrap_subset",
        "train": rows_for_temperatures(args.data_root, train_temperatures),
        "test": rows_for_temperatures(args.data_root, held_temperatures),
    }
    full_split_path = args.output_root / "sio2_three_phase_full_range.json"
    subset_split_path = args.output_root / "sio2_temp_extrap_subset.json"
    full_split_path.write_text(json.dumps(full_split, indent=2) + "\n", encoding="utf-8")
    subset_split_path.write_text(json.dumps(subset_split, indent=2) + "\n", encoding="utf-8")

    records: list[dict[str, object]] = []
    for basis in ("polynomial", "tlog_polynomial"):
        for seed_dir in sorted((args.checkpoint_root / basis).glob("seed*")):
            seed = int(seed_dir.name.removeprefix("seed"))
            config = json.loads((seed_dir / "config.json").read_text(encoding="utf-8"))
            raw = evaluate_model(seed_dir / "model.pth", args.legacy_data_root, args.device, args.batch_size)
            formal_c = calibration(raw, args.data_root, train_temperatures)
            held_c = calibration(raw, args.data_root, held_temperatures)
            variants = {"formal_train_window": formal_c, "heldout_window_audit": held_c}
            result: dict[str, object] = {
                "basis": basis,
                "seed": seed,
                "checkpoint": external_path(seed_dir / "model.pth", args.checkpoint_root),
                "config_training_systems": ["external/training_data/" + Path(item).name for item in config["training"]["training_data"]["systems"]],
                "model_internal_c_system_field": config["model"]["fitting_net"].get("c_system"),
                "calibration_constants_eV_per_atom": variants,
                "metrics": {},
                "provenance": {
                    "git_commit": repository_commit(repo),
                    "split": "data/processed/splits/temp_extrap.json",
                    "split_sha256": sha256(args.split),
                    "predictor_config": str(args.predictor_config.relative_to(repo)),
                    "predictor_config_sha256": sha256(args.predictor_config),
                    "checkpoint_path": external_path(seed_dir / "model.pth", args.checkpoint_root),
                    "checkpoint_sha256": sha256(seed_dir / "model.pth"),
                    "evaluator_version": "fes_bench.eval.run:" + sha256(repo / "src/lib/fes_bench/eval/run.py"),
                    "timestamp_utc": datetime.now(timezone.utc).isoformat(),
                },
            }
            for label, constant in variants.items():
                curve_path = args.output_root / "curves" / f"{basis}_seed{seed}_{label}.json"
                curve_path.parent.mkdir(parents=True, exist_ok=True)
                curve_path.write_text(json.dumps({"predictions": shifted(raw, constant)}, indent=2) + "\n", encoding="utf-8")
                result["metrics"][label] = {
                    "full_range": invoke_eval(repo, curve_path, full_split_path, args.data_root, args.output_root / "eval", seed),
                    "temp_extrap": invoke_eval(repo, curve_path, subset_split_path, args.data_root, args.output_root / "eval", seed),
                }
            records.append(result)
            print(f"completed {basis} seed {seed}", flush=True)
    payload = {
        "scope": "four-phase historic SiO2 FES checkpoints; reported metrics restricted to canonical quartz/cristobalite/p63mmc subset",
        "calibration_windows_K": {
            "formal_train": [min(train_temperatures), max(train_temperatures)],
            "heldout_audit": [min(held_temperatures), max(held_temperatures)],
            "source": "frozen data/processed/splits/temp_extrap.json",
        },
        "provenance": {
            "git_commit": repository_commit(repo),
            "split": "data/processed/splits/temp_extrap.json",
            "split_sha256": sha256(args.split),
            "predictor_config": str(args.predictor_config.relative_to(repo)),
            "predictor_config_sha256": sha256(args.predictor_config),
            "checkpoint_path": "external/checkpoints/continuous_sio2/<basis>/seed<seed>/model.pth",
            "checkpoint_sha256": "per checkpoint_record provenance",
            "evaluator_version": "fes_bench.eval.run:" + sha256(repo / "src/lib/fes_bench/eval/run.py"),
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        },
        "checkpoint_records": records,
    }
    output = args.output_root / "sio2_reanalysis.json"
    if args.compare_to is not None:
        if args.tolerance_config is None:
            parser.error("--compare-to requires --tolerance-config")
        payload["provenance"]["reproduction_check"] = compare_records(
            payload, args.compare_to.resolve(), args.tolerance_config.resolve(), args.output_root
        )
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    if args.compare_to is not None:
        write_canonical_meta(output, args.compare_to.resolve(), args.tolerance_config.resolve(), args.output_root, payload["provenance"]["reproduction_check"])
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
