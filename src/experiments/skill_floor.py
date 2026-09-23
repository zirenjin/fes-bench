"""Audit skill relative to each explicitly named floor — plan experiment E4."""

from __future__ import annotations

import argparse
import csv
import json
from itertools import combinations
from pathlib import Path

from _common import write_audit


def folds(body: dict) -> dict[str, dict]:
    return body.get("folds", {"all": body})


def e1_constant_floor(root: Path) -> tuple[dict[str, float], set[str], list[Path]]:
    """Derive pairwise training-mean floors from the configured reference grid."""
    sample = next(root.glob("result/experiments/crossing_reevaluation/raw_runs/*/seed_*/metrics.json"))
    sample_metrics = json.loads(sample.read_text(encoding="utf-8"))["metrics"]
    pair_names = list(sample_metrics["metrics"]["formal_train_window"]["full_range"]["pairs"])
    system_id = pair_names[0].split(":", 1)[0]
    config_path = next(
        path for path in sorted((root / "configs/systems").glob("*.yaml"))
        if json.loads(path.read_text(encoding="utf-8")).get("system") == system_id
    )
    config = json.loads(config_path.read_text(encoding="utf-8"))
    system_json_path = root / config["system_json"]
    system_meta = json.loads(system_json_path.read_text(encoding="utf-8"))
    phases = list(system_meta["phases"])
    manifest = json.loads((root / "configs/raw_runs.json").read_text(encoding="utf-8"))
    e1_spec = next(item for item in manifest["runs"] if item.get("split") == "e1_full_grid")
    floor_split_name = e1_spec["floor_split"]
    split_path = next(
        path for path in sorted((root / "configs/splits").glob("*.yaml"))
        if json.loads(path.read_text(encoding="utf-8")).get("split", json.loads(path.read_text(encoding="utf-8")).get("name")) == floor_split_name
    )
    split_config = json.loads(split_path.read_text(encoding="utf-8"))
    split = json.loads((root / split_config["frozen_json"]).read_text(encoding="utf-8"))
    train = {}
    for item in split["train"]:
        train.setdefault((item["system"], item["phase"]), set()).add(int(item["T_index"]))
    values = {}
    inputs = [config_path, system_json_path, split_path, root / split_config["frozen_json"], root / "configs/raw_runs.json"]
    for phase in phases:
        path = root / "data/processed" / system_id / phase / "reference_G.csv"
        inputs.append(path)
        with path.open(encoding="utf-8", newline="") as handle:
            values[phase] = [float(row["G_eV_per_atom"]) for row in csv.DictReader(handle)]
    floors = {}
    for left, right in combinations(phases, 2):
        pair = f"{system_id}:{left}_minus_{right}"
        delta = [a - b for a, b in zip(values[left], values[right])]
        indices = [i for i, _ in enumerate(delta) if i in train.get((system_id, left), set()) and i in train.get((system_id, right), set())]
        constant = sum(delta[i] for i in indices) / len(indices)
        floor = sum(abs(item - constant) for item in delta) / len(delta)
        floors[pair] = floor
        floors[f"{system_id}:{right}_minus_{left}"] = floor
    return floors, set(floors), inputs


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--repo-root", type=Path, default=Path(".")); args = parser.parse_args(); root = args.repo_root.resolve()
    rows, inputs = [], []
    split_dirs = sorted(path for path in (root / "result/experiments/external_baselines").glob("raw_*") if path.is_dir())
    for split_dir in split_dirs:
        split = split_dir.name.removeprefix("raw_")
        constants = split_dir / "constant_delta_g/seed_none/metrics.json"; global_floor = split_dir / "global_mean_delta_g/seed_none/metrics.json"
        if not constants.exists() or not global_floor.exists(): continue
        inputs.extend([constants, global_floor]); constant_folds = folds(json.loads(constants.read_text(encoding="utf-8"))["metrics"]); global_folds = folds(json.loads(global_floor.read_text(encoding="utf-8"))["metrics"])
        for path in sorted(split_dir.glob("*/seed_none/metrics.json")):
            predictor = path.parts[-3]
            if predictor in {"constant_delta_g", "global_mean_delta_g"}: continue
            inputs.append(path); model_folds = folds(json.loads(path.read_text(encoding="utf-8"))["metrics"]); denominator_name, floor_folds = ("constant_delta_g", constant_folds) if predictor in {"zero", "reference", "reference_noise_0.005", "constant_offset_0.005", "root_stability_iid_0.005"} else ("global_mean_delta_g", global_folds)
            for fold, model_fold in sorted(model_folds.items()):
                floor_fold = floor_folds.get(fold, {})
                model_pairs = model_fold.get("pairs", {}); floor_pairs = floor_fold.get("pairs", {}); matched = [(float(item["delta_G_MAE_eV_per_atom"]), float(floor_pairs[pair]["delta_G_MAE_eV_per_atom"])) for pair, item in model_pairs.items() if item.get("delta_G_MAE_eV_per_atom") is not None and floor_pairs.get(pair, {}).get("delta_G_MAE_eV_per_atom") is not None]
                model_value = sum(item[0] for item in matched) / len(matched) if matched else None; floor_value = sum(item[1] for item in matched) / len(matched) if matched else None
                skill = None if model_value is None or not floor_value else 1.0 - model_value / floor_value
                rows.append({"split": split, "fold": fold, "predictor": predictor, "floor_predictor": denominator_name, "aggregation": "matched_pairs", "model_MAE_eV_per_atom": model_value, "floor_MAE_eV_per_atom": floor_value, "skill": skill, "skill_status": "NEGATIVE" if skill is not None and skill < 0 else "non-negative" if skill is not None else "unavailable"})

    # E1 has ten frozen checkpoints. Recompute the requested pair-level
    # constant floor from the configured reference grid, including both
    # all-pair and no-reference-crossing-pair aggregates.
    floors, _, floor_inputs = e1_constant_floor(root)
    inputs.extend(floor_inputs)
    for path in sorted((root / "result/experiments/crossing_reevaluation/raw_runs").glob("*/seed_*/metrics.json")):
        inputs.append(path)
        record = json.loads(path.read_text(encoding="utf-8"))["metrics"]
        basis = record["basis"]; seed = record["seed"]
        pairs = record["metrics"]["formal_train_window"]["full_range"]["pairs"]
        matched = [(pair, float(item["delta_G_MAE_eV_per_atom"]), floors[pair]) for pair, item in pairs.items() if pair in floors and item.get("delta_G_MAE_eV_per_atom") is not None]
        no_reference = [(pair, model, floor) for pair, model, floor in matched if pairs[pair].get("reference_Tc_K")]
        for aggregation, selected in (("include_all_pairs", matched), ("exclude_no_reference_crossing", no_reference)):
            model_value = sum(item[1] for item in selected) / len(selected) if selected else None
            floor_value = sum(item[2] for item in selected) / len(selected) if selected else None
            skill = None if model_value is None or not floor_value else 1.0 - model_value / floor_value
            rows.append({"split": "e1_full_grid", "fold": str(seed), "predictor": basis, "floor_predictor": "constant_delta_g", "aggregation": aggregation, "model_MAE_eV_per_atom": model_value, "floor_MAE_eV_per_atom": floor_value, "skill": skill, "skill_status": "NEGATIVE" if skill is not None and skill < 0 else "non-negative" if skill is not None else "unavailable"})
    write_audit(root, "skill_floor", ["split", "fold", "predictor", "floor_predictor", "aggregation", "model_MAE_eV_per_atom", "floor_MAE_eV_per_atom", "skill", "skill_status"], rows, inputs, [])
    return 0


if __name__ == "__main__": raise SystemExit(main())
