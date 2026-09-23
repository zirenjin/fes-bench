"""Predictor comparison across the full metric set — plan Table 5.

Nine rows for actual predictors (trivial floors / external methods /
controlled ablations), one column per metric. The global-mean floor is kept
for the skill-floor audit denominator but is not a model/predictor row here.
Reads result/experiments/ only; never runs a model.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from common import aggregate, csv_write, meta_write, raw_runs


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--split", required=True)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--output-root", type=Path, default=Path("result/tables"))
    args = parser.parse_args(); root = args.repo_root.resolve(); output = args.output_root if args.output_root.is_absolute() else root / args.output_root
    rows, inputs, missing = [], [], []
    for path, body in raw_runs(root, args.split):
        inputs.append(path)
        predictor = path.parts[-3]
        if predictor == "global_mean_delta_g":
            continue
        summary = aggregate(body["metrics"])
        row = {"predictor": predictor, **summary}
        if summary["G_MAE_eV_per_atom"] is None:
            row["G_MAE_eV_per_atom"] = "n/a:pair_only_predictor"
        rows.append(row)
        if summary["G_MAE_eV_per_atom"] is None:
            missing.append({"field": f"{predictor}.G_MAE_eV_per_atom", "reason": "pair-only predictor or unavailable phase-energy metric"})
    rows.sort(key=lambda row: row["predictor"])
    fields = ["predictor", "G_MAE_eV_per_atom", "delta_G_MAE_eV_per_atom", "delta_G_RMSE_eV_per_atom", "sign_accuracy", "Tc_error_K", "false_crossings", "missed_crossings"]
    csv_path = output / f"predictor_comparison_{args.split}.csv"; csv_write(csv_path, fields, rows); meta_write(root, csv_path, inputs, missing)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
