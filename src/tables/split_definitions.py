"""Frozen split definitions and frame counts — plan Table 3."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from common import csv_write, inventory, meta_write


def _crossing_membership(root: Path, split: str, fold: str, fold_data: dict) -> str:
    """Return a reproducible, human-readable status for every reference crossing."""
    labels = {}
    allowed_systems = {str(item["system"]) for subset in ("train", "test") for item in fold_data.get(subset, [])}
    for subset in ("train", "test"):
        for item in fold_data.get(subset, []):
            labels[(item["system"], item["phase"], int(item["T_index"]))] = subset
    statuses = []
    for crossing_path in sorted((root / "data/processed").glob("*/reference_crossings.json")):
        crossing_data = json.loads(crossing_path.read_text(encoding="utf-8"))
        system = crossing_data.get("system", crossing_path.parent.name)
        if system not in allowed_systems:
            continue
        grid = crossing_data.get("grid", {})
        t_min = float(grid.get("T_min_K", 0.0)); step = float(grid.get("step_K", 1.0))
        for pair, pair_data in sorted(crossing_data.get("pairs", {}).items()):
            for crossing in pair_data.get("crossings", []):
                temperature = float(crossing["T_K"])
                index = int(round((temperature - t_min) / step))
                phases = (pair_data.get("left"), pair_data.get("right"))
                observed = [labels.get((system, phase, index)) for phase in phases if phase]
                observed = [value for value in observed if value is not None]
                if "test" in observed:
                    status = "test"
                elif observed and all(value == "train" for value in observed):
                    status = "train"
                else:
                    # A phase/system held out by construction has no train row;
                    # its crossing is therefore in the test portion.
                    status = "test" if not observed else "boundary"
                statuses.append(f"{system}:{pair}@{temperature:g}K={status}")
    if split == "temp_extrap":
        unexpected = [item for item in statuses if not item.endswith("=test")]
        if unexpected:
            raise SystemExit("temp_extrap crossing_in_test violation: " + ", ".join(unexpected))
    return ";".join(statuses)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--output-root", type=Path, default=Path("result/tables"))
    args = parser.parse_args()
    root = args.repo_root.resolve(); output = args.output_root if args.output_root.is_absolute() else root / args.output_root
    split_configs = {}
    for config_path in sorted((root / "configs/splits").glob("*.yaml")):
        import json
        config = json.loads(config_path.read_text(encoding="utf-8"))
        split_configs[config.get("split", config.get("name", config_path.stem))] = config
    rows = []
    inputs = [root / "result/experiments/data_prep/raw_inventory/inventory/seed_none/metrics.json"]
    for split, data in sorted(inventory(root)["splits"].items()):
        frozen = root / "data/processed/splits" / f"{split}.json"
        if frozen.exists():
            split_json = json.loads(frozen.read_text(encoding="utf-8")); inputs.append(frozen)
        else:
            split_json = {}
        for crossing_path in sorted((root / "data/processed").glob("*/reference_crossings.json")):
            inputs.append(crossing_path)
        for fold, counts in sorted(data.get("folds", {}).items()):
            held_out_unit = split_configs.get(split, {}).get("held_out_unit", "temperature")
            fold_data = split_json.get("folds", {}).get(fold, {}) if "folds" in split_json else split_json
            rows.append({"split": split, "fold": fold, "held_out_unit": held_out_unit, "train_frames": counts.get("train", 0), "test_frames": counts.get("test", 0), "crossing_in_test": _crossing_membership(root, split, fold, fold_data), "sha256": data.get("sha256", ""), "commit": data.get("git_commit", "")})
        if not data.get("folds"):
            held_out_unit = split_configs.get(split, {}).get("held_out_unit", "temperature")
            rows.append({"split": split, "fold": "—", "held_out_unit": held_out_unit, "train_frames": data.get("train", 0), "test_frames": data.get("test", 0), "crossing_in_test": _crossing_membership(root, split, "—", split_json), "sha256": data.get("sha256", ""), "commit": data.get("git_commit", "")})
    missing = []
    csv_path = output / "split_definitions.csv"; csv_write(csv_path, ["split", "fold", "held_out_unit", "train_frames", "test_frames", "crossing_in_test", "sha256", "commit"], rows)
    meta_write(root, csv_path, inputs, missing)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
