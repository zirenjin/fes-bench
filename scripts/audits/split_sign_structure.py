"""Audit sign distributions induced by frozen splits — plan experiment E6."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from common import write_audit


def table(root: Path, system: str, phase: str) -> dict[int, float]:
    with (root / "data" / system / phase / "reference_G.csv").open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return {index: float(row["G_eV_per_atom"]) for index, row in enumerate(rows)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--repo-root", type=Path, default=Path(".")); args = parser.parse_args(); root = args.repo_root.resolve()
    rows, inputs = [], []
    for split_path in sorted((root / "splits").glob("*.json")):
        split = json.loads(split_path.read_text(encoding="utf-8")); inputs.append(split_path)
        folds = split.get("folds", {"all": split})
        for fold, body in sorted(folds.items()):
            all_rows = body.get("train", []) + body.get("test", [])
            keys = sorted({(item["system"], item["phase"]) for item in all_rows})
            for system in sorted({item[0] for item in keys}):
                phases = [phase for item_system, phase in keys if item_system == system]
                curves = {phase: table(root, system, phase) for phase in phases}
                for left_index, left in enumerate(phases):
                    for right in phases[left_index + 1 :]:
                        shared = sorted(set(curves[left]) & set(curves[right]))
                        train_idx = {item["T_index"] for item in body.get("train", []) if item["system"] == system and item["phase"] in {left, right}}
                        test_idx = {item["T_index"] for item in body.get("test", []) if item["system"] == system and item["phase"] in {left, right}}
                        for subset, indices in (("train", train_idx), ("test", test_idx)):
                            values = [curves[left][index] - curves[right][index] for index in shared if index in indices]
                            if not values: continue
                            rows.append({"split": split_path.stem, "fold": fold, "system": system, "pair": f"{left}_minus_{right}", "subset": subset, "n_points": len(values), "positive_fraction": sum(value > 0 for value in values) / len(values), "negative_fraction": sum(value < 0 for value in values) / len(values), "zero_fraction": sum(value == 0 for value in values) / len(values)})
    write_audit(root, "split_sign_structure", ["split", "fold", "system", "pair", "subset", "n_points", "positive_fraction", "negative_fraction", "zero_fraction"], rows, inputs, [])
    return 0


if __name__ == "__main__": raise SystemExit(main())
