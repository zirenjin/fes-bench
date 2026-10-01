"""Fit reference ΔG slopes around every frozen crossing — Figure 1 prerequisite."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
from datetime import datetime, timezone

import numpy as np


WINDOW_K = 25.0


def fit_slope(temperatures: np.ndarray, delta_g: np.ndarray, crossing: float) -> dict[str, object]:
    """Fit ΔG(T) in the clipped Tc±25 K window and return auditable statistics."""
    lower = max(float(temperatures.min()), float(crossing) - WINDOW_K)
    upper = min(float(temperatures.max()), float(crossing) + WINDOW_K)
    mask = (temperatures >= lower) & (temperatures <= upper)
    x = temperatures[mask]
    y = delta_g[mask]
    x_centered = x - float(x.mean())
    design = np.column_stack((x_centered, np.ones_like(x_centered)))
    coefficient, *_ = np.linalg.lstsq(design, y, rcond=None)
    residual = y - design @ coefficient
    dof = max(1, len(x) - 2)
    residual_variance = float(residual @ residual) / dof
    denominator = float(np.sum(x_centered * x_centered))
    stderr = float(np.sqrt(residual_variance / denominator)) if denominator else float("nan")
    total = float(np.sum((y - y.mean()) ** 2))
    r_squared = 1.0 - float(residual @ residual) / total if total else 1.0
    return {
        "slope_eV_per_atom_per_K": float(coefficient[0]),
        "slope_stderr": stderr,
        "fit_window_K": [float(x.min()), float(x.max())],
        "n_points": int(len(x)),
        "r_squared": r_squared,
    }


def process_system(root: Path, system: str) -> tuple[list[dict[str, object]], list[Path]]:
    data_root = root / "data/processed" / system
    crossing_path = data_root / "reference_crossings.json"
    payload = json.loads(crossing_path.read_text(encoding="utf-8"))
    tables: dict[str, dict[float, float]] = {}
    for phase_dir in sorted(data_root.iterdir()):
        table_path = phase_dir / "reference_G.csv"
        if not table_path.exists():
            continue
        with table_path.open(encoding="utf-8", newline="") as handle:
            tables[phase_dir.name] = {float(row["T_K"]): float(row["G_eV_per_atom"]) for row in csv.DictReader(handle)}
    rows: list[dict[str, object]] = []
    for pair_name, pair in payload.get("pairs", {}).items():
        left, right = str(pair["left"]), str(pair["right"])
        temperatures = np.array(sorted(set(tables[left]) & set(tables[right])), dtype=float)
        delta_g = np.array([tables[left][temperature] - tables[right][temperature] for temperature in temperatures], dtype=float)
        for crossing_index, crossing in enumerate(pair.get("crossings", []), start=1):
            old = crossing.get("previous_slope_eV_per_atom_per_K", crossing.get("slope_eV_per_atom_per_K"))
            if old is None and crossing.get("slope_meV_per_atom_per_K") is not None:
                old = float(crossing["slope_meV_per_atom_per_K"]) * 1.0e-3
            fitted = fit_slope(temperatures, delta_g, float(crossing["T_K"]))
            crossing["previous_slope_eV_per_atom_per_K"] = old
            crossing.update(fitted)
            rows.append({"system": system, "pair": pair_name, "crossing_index": crossing_index, "previous_slope_eV_per_atom_per_K": old, **fitted})
    crossing_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return rows, [crossing_path]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    args = parser.parse_args()
    root = args.repo_root.resolve()
    rows: list[dict[str, object]] = []
    for system_dir in sorted((root / "data/processed").iterdir()):
        if (system_dir / "reference_crossings.json").exists():
            system_rows, _ = process_system(root, system_dir.name)
            rows.extend(system_rows)
    output = root / "result/experiments/reference_crossing_slopes/findings.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        fields = sorted({key for row in rows for key in row})
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)
    def sha256(path: Path) -> str:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()
    try:
        commit = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "unknown"
    low_r2 = [row for row in rows if float(row["r_squared"]) < 0.95]
    sio2_rows = [row for row in rows if row["system"] == "sio2"]
    meta = {
        "script": "src/experiments/data_prep/refit_reference_crossing_slopes.py",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": commit,
        "window_K": WINDOW_K,
        "fit": "ordinary least squares on reference Delta G versus temperature; window clipped to the reference grid",
        "inputs": [
            {"path": str((root / "data/processed" / system / "reference_crossings.json").relative_to(root)),
             "sha256": sha256(root / "data/processed" / system / "reference_crossings.json")}
            for system in sorted({str(row["system"]) for row in rows})
        ],
        "low_r_squared": low_r2,
        "sio2_old_new_comparison": [
            {"system": row["system"], "pair": row["pair"], "old": row["previous_slope_eV_per_atom_per_K"],
             "new": row["slope_eV_per_atom_per_K"],
             "relative_change": (float(row["slope_eV_per_atom_per_K"]) - float(row["previous_slope_eV_per_atom_per_K"]))
             / float(row["previous_slope_eV_per_atom_per_K"])}
            for row in sio2_rows
        ],
    }
    (output.parent / "findings.meta.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
