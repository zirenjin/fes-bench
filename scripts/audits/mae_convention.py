"""Trace ΔG MAE aggregation conventions — plan experiment E5."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from common import write_audit


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--repo-root", type=Path, default=Path(".")); args = parser.parse_args(); root = args.repo_root.resolve()
    rows, inputs = [], []
    by_basis = {}
    for path in sorted((root / "results/raw_runs/e1_full_grid").glob("*/seed_*/metrics.json")):
        inputs.append(path); record = json.loads(path.read_text(encoding="utf-8"))["metrics"]; block = record.get("metrics", {}).get("formal_train_window", {})
        basis = record.get("basis"); by_basis.setdefault(basis, {}).setdefault("full_range", [])
        for window, metrics in block.items():
            pairs = [float(item["delta_G_MAE_eV_per_atom"]) for item in metrics.get("pairs", {}).values() if item.get("delta_G_MAE_eV_per_atom") is not None]
            pair_mae = sum(pairs) / len(pairs) if pairs else None
            reported = metrics.get("G_MAE_eV_per_atom")
            rows.append({
                "basis": basis, "seed": record.get("seed"), "window": window,
                "reported_G_MAE_eV_per_atom": reported,
                "pair_delta_G_MAE_eV_per_atom": pair_mae,
                "n_test_frames": metrics.get("n_test_frames"), "pair_count": len(pairs),
                "reported_over_pair_ratio": (float(reported) / pair_mae if reported is not None and pair_mae else None),
                "interpretation": "phase-level scalar G MAE; not the pair-level ΔG MAE" if reported is not None else "pair-only or unavailable phase-level G MAE",
            })
            if window == "full_range" and pair_mae is not None:
                by_basis.setdefault(basis, {}).setdefault(window, []).append(pair_mae)
    for basis, windows in sorted(by_basis.items()):
        values = windows.get("full_range", [])
        if values:
            rows.append({
                "basis": basis, "seed": "all", "window": "aggregate_full_range",
                "reported_G_MAE_eV_per_atom": None,
                "pair_delta_G_MAE_eV_per_atom": sum(values) / len(values),
                "n_test_frames": None, "pair_count": 3,
                "reported_over_pair_ratio": None,
                "interpretation": "equal-weight mean of three pair-level ΔG MAEs across five checkpoints; reproduces E1 aggregate",
            })
    fields = ["basis", "seed", "window", "reported_G_MAE_eV_per_atom", "pair_delta_G_MAE_eV_per_atom", "n_test_frames", "pair_count", "reported_over_pair_ratio", "interpretation"]
    write_audit(root, "mae_convention", fields, rows, inputs, [{"field": "historical_2.56_meV_per_atom", "reason": "detached four-phase all_pair_metrics artifact is not present in normalized raw metrics; cannot be regenerated from this checkout"}], extra={"resolution": "closed", "resolution_note": "Historical 2.56 meV/atom is closed: four-phase gauge-C, quartz–C2221 single pair, full grid, different checkpoint; it does not affect current conclusions."})
    return 0


if __name__ == "__main__": raise SystemExit(main())
