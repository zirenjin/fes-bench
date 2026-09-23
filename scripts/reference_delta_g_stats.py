"""Summarize reference ΔG ranges for every frozen split and fold."""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


def table(data_root: Path, system: str, phase: str) -> tuple[np.ndarray, np.ndarray]:
    with (data_root / system / phase / "reference_G.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return np.array([float(row["T_K"]) for row in rows]), np.array([float(row["G_eV_per_atom"]) for row in rows])


def curves(rows: list[dict[str, object]], data_root: Path) -> list[dict[str, object]]:
    indexes: dict[tuple[str, str], set[int]] = defaultdict(set)
    for row in rows:
        indexes[(str(row["system"]), str(row["phase"]))].add(int(row["T_index"]))
    result: list[dict[str, object]] = []
    for system in sorted({system for system, _ in indexes}):
        phases = sorted(phase for item_system, phase in indexes if item_system == system)
        loaded = {phase: table(data_root, system, phase) for phase in phases}
        for i, left in enumerate(phases):
            for right in phases[i + 1 :]:
                shared = sorted(indexes[(system, left)] & indexes[(system, right)])
                if not shared:
                    continue
                _, left_g = loaded[left]
                temperatures, right_g = loaded[right]
                delta = left_g[shared] - right_g[shared]
                result.append(
                    {
                        "system": system,
                        "pair": f"{left}_minus_{right}",
                        "n_frames": len(shared),
                        "temperature_range_K": [float(temperatures[shared[0]]), float(temperatures[shared[-1]])],
                        "mean_eV_per_atom": float(np.mean(delta)),
                        "std_eV_per_atom": float(np.std(delta)),
                        "min_eV_per_atom": float(np.min(delta)),
                        "max_eV_per_atom": float(np.max(delta)),
                    }
                )
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--splits-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--include-system-full-grids",
        action="store_true",
        help="add one all-grid record for systems absent from the frozen split files",
    )
    args = parser.parse_args()
    records = []
    for split_path in sorted(args.splits_root.glob("*.json")):
        split = json.loads(split_path.read_text(encoding="utf-8"))
        if "folds" in split:
            for fold_name, fold in split["folds"].items():
                for subset in ("train", "test"):
                    rows = fold.get(subset, [])
                    records.extend(
                        {"split": split_path.stem, "fold": fold_name, "subset": subset, **item}
                        for item in curves(rows, args.data_root)
                    )
                partner_rows = fold.get("train", []) + fold.get("test", [])
                records.extend(
                    {"split": split_path.stem, "fold": fold_name, "subset": "test_plus_train_partner", **item}
                    for item in curves(partner_rows, args.data_root)
                )
        else:
            for subset in ("train", "test"):
                rows = split.get(subset, [])
                records.extend(
                    {"split": split_path.stem, "fold": "—", "subset": subset, **item}
                    for item in curves(rows, args.data_root)
                )
            partner_rows = split.get("train", []) + split.get("test", [])
            records.extend(
                {"split": split_path.stem, "fold": "—", "subset": "test_plus_train_partner", **item}
                for item in curves(partner_rows, args.data_root)
            )
    if args.include_system_full_grids:
        existing_systems = {str(item["system"]) for item in records}
        for system_path in sorted(args.data_root.glob("*/system.json")):
            system = system_path.parent.name
            if system in existing_systems:
                continue
            system_meta = json.loads(system_path.read_text(encoding="utf-8"))
            rows = []
            for phase in system_meta.get("phases", []):
                with (args.data_root / system / phase / "reference_G.csv").open(newline="", encoding="utf-8") as handle:
                    n_points = sum(1 for _ in handle) - 1
                rows.extend({"system": system, "phase": phase, "T_index": index} for index in range(n_points))
            records.extend(
                {"split": "reference_full_grid", "fold": "—", "subset": "all", **item}
                for item in curves(rows, args.data_root)
            )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    payload = {"units": "eV/atom", "definition": "ΔG = G(left) − G(right); std uses population normalization (ddof=0)", "records": records}
    args.output.with_suffix(".json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    lines = ["# Reference ΔG statistics", "", payload["definition"], "", "| split | fold | subset | system | pair | n | T range (K) | std | min | max |", "|---|---|---|---|---|---:|---|---:|---:|---:|"]
    for item in records:
        lines.append(
            f"| {item['split']} | {item['fold']} | {item['subset']} | {item['system']} | {item['pair']} | {item['n_frames']} | "
            f"{item['temperature_range_K'][0]:.0f}–{item['temperature_range_K'][1]:.0f} | {item['std_eV_per_atom']:.6g} | "
            f"{item['min_eV_per_atom']:.6g} | {item['max_eV_per_atom']:.6g} |"
        )
    args.output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
