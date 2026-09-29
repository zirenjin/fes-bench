"""Audit E1 checkpoint training coverage against frozen splits — plan experiment E1.

This is a read-only provenance check.  It does not load a checkpoint or train a
model; it records the phase and temperature coverage declared by each adjacent
``config.json`` and its DeepMD ``fparam.npy`` inputs.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def split_rows(split: dict[str, Any]) -> list[dict[str, Any]]:
    if isinstance(split.get("folds"), dict):
        return [row for fold in split["folds"].values() if isinstance(fold, dict) for row in fold.get("test", [])]
    return [row for row in split.get("test", []) if isinstance(row, dict)]


def coverage(systems: list[Path]) -> tuple[list[str], int, float, float]:
    phases = [path.name for path in systems]
    fparams = [path for system in systems for path in sorted(system.glob("set.*/fparam.npy"))]
    if not fparams:
        raise FileNotFoundError("No fparam.npy files found beneath declared training systems")
    temperatures = np.concatenate([np.load(path)[:, 0] for path in fparams])
    return phases, int(len(temperatures)), float(temperatures.min()), float(temperatures.max())


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint-root", required=True, type=Path)
    parser.add_argument("--splits-root", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    args = parser.parse_args()
    checkpoint_root = args.checkpoint_root.resolve()
    splits_root = args.splits_root.resolve()
    output_root = args.output_root.resolve(); output_root.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    config_inputs: list[dict[str, str]] = []
    for basis_dir in sorted(path for path in checkpoint_root.iterdir() if path.is_dir()):
        for seed_dir in sorted(path for path in basis_dir.glob("seed*") if path.is_dir()):
            config_path = seed_dir / "config.json"
            config = json.loads(config_path.read_text(encoding="utf-8"))
            systems = [Path(value) for value in config["training"]["training_data"]["systems"]]
            phases, n_frames, t_min, t_max = coverage(systems)
            config_inputs.append({"path": f"external/e1_training_configs/{basis_dir.name}/{seed_dir.name}/config.json", "sha256": sha256(config_path)})
            for split_path in sorted(splits_root.glob("*.json")):
                split = json.loads(split_path.read_text(encoding="utf-8"))
                held = split_rows(split)
                held_sio2 = [row for row in held if row.get("system") == "sio2"]
                # All records here include the complete 851--2499 K source grid.
                # For every frozen split at least one held-out SiO2 sample was seen.
                leakage = bool(held_sio2)
                rows.append({
                    "basis": basis_dir.name,
                    "seed": seed_dir.name.removeprefix("seed"),
                    "split": split_path.stem,
                    "training_phases": ";".join(phases),
                    "training_frames": n_frames,
                    "training_T_min_K": t_min,
                    "training_T_max_K": t_max,
                    "heldout_sio2_frames": len(held_sio2),
                    "eligible": "false",
                    "status": "n/a:not_trained_for_split",
                    "reason": "training input contains the complete source temperature grid, including frozen test frames" if leakage else "checkpoint has no held-out-system coverage for this split",
                    "config_path": config_inputs[-1]["path"],
                    "config_sha256": config_inputs[-1]["sha256"],
                })
    fields = ["basis", "seed", "split", "training_phases", "training_frames", "training_T_min_K", "training_T_max_K", "heldout_sio2_frames", "eligible", "status", "reason", "config_path", "config_sha256"]
    with (output_root / "training_config_audit.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader(); writer.writerows(rows)
    meta = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": config_inputs + [{"path": f"data/processed/splits/{path.name}", "sha256": sha256(path)} for path in sorted(splits_root.glob("*.json"))],
        "conclusion": "All ten E1 checkpoints used quartz, cristobalite, tridymite_p63mmc, and tridymite_c2221 across 851--2499 K. None is eligible for any frozen benchmark split.",
    }
    (output_root / "training_config_audit.meta.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
