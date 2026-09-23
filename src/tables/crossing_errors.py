"""Crossing errors by phase pair and predictor — plan Table 6."""

from __future__ import annotations

import argparse
from pathlib import Path

from common import csv_write, meta_write, pair_records, raw_runs


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--split", required=True)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--output-root", type=Path, default=Path("result/tables"))
    args = parser.parse_args(); root = args.repo_root.resolve(); output = args.output_root if args.output_root.is_absolute() else root / args.output_root
    rows, inputs, missing = [], [], []
    for path, body in raw_runs(root, args.split):
        inputs.append(path); predictor = path.parts[-3]
        for pair, record in pair_records(body["metrics"]):
            references = record.get("reference_Tc_K", [])
            errors = record.get("Tc_error_K", [])
            slopes = record.get("crossing_slope_eV_per_atom_per_K", [])
            n = max(len(references), len(errors), len(slopes), 1)
            reference_value = references if references else (
                "n/a:no_training_phase"
                if record.get("status") == "unavailable_without_training_phase"
                else "n/a:no_reference_crossing"
            )
            for index in range(n):
                rows.append({"predictor": predictor, "fold": record.get("fold", ""), "pair": pair, "crossing_index": index + 1, "reference_Tc_K": reference_value[index] if isinstance(reference_value, list) and index < len(reference_value) else (reference_value if index == 0 else ""), "Tc_error_K": errors[index] if index < len(errors) else "", "crossing_slope_eV_per_atom_per_K": slopes[index] if index < len(slopes) else "", "Tc_err_from_dG_K": (record.get("Tc_err_from_dG_K", [])[index] if index < len(record.get("Tc_err_from_dG_K", [])) else ""), "false_crossings": record.get("false_crossings", ""), "missed_crossings": record.get("missed_crossings", "")})
            if not references:
                reason = (
                    "predictor unavailable without training phase"
                    if record.get("status") == "unavailable_without_training_phase"
                    else "structural N/A: pair has no reference crossing on the frozen grid"
                )
                missing.append({"field": f"{predictor}:{pair}.reference_Tc_K", "reason": reason})
    fields = ["predictor", "fold", "pair", "crossing_index", "reference_Tc_K", "Tc_error_K", "crossing_slope_eV_per_atom_per_K", "Tc_err_from_dG_K", "false_crossings", "missed_crossings"]
    csv_path = output / f"crossing_errors_{args.split}.csv"; csv_write(csv_path, fields, rows); meta_write(root, csv_path, inputs, missing)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
