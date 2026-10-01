"""Crossing errors by phase pair and predictor — plan Table 6."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from common import csv_write, meta_write, pair_records, raw_runs


def reference_crossing(root: Path, pair: str, index: int) -> dict[str, object]:
    """Read the fitted reference-slope record for one system pair/root."""
    if ":" not in pair:
        return {}
    system, pair_name = pair.split(":", 1)
    path = root / "data/processed" / system / "reference_crossings.json"
    if not path.exists():
        return {}
    payload = json.loads(path.read_text(encoding="utf-8"))
    pair_payload = payload.get("pairs", {}).get(pair_name)
    if pair_payload is None and "_minus_" in pair_name:
        left, right = pair_name.split("_minus_", 1)
        pair_payload = payload.get("pairs", {}).get(f"{right}_minus_{left}")
    roots = (pair_payload or {}).get("crossings", [])
    return roots[index] if index < len(roots) and isinstance(roots[index], dict) else {}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--split", required=True)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--output-root", type=Path, default=Path("result/tables"))
    args = parser.parse_args(); root = args.repo_root.resolve(); output = args.output_root if args.output_root.is_absolute() else root / args.output_root
    rows, overlap_rows, inputs, missing = [], [], [], []
    reference_inputs: set[Path] = set()
    for path, body in raw_runs(root, args.split):
        inputs.append(path); predictor = path.parts[-3]
        records = [("main", pair, record) for pair, record in pair_records(body["metrics"])
                   if args.split != "system_loso" or record.get("fold", "") in {"hf", "ti", "zr"}]
        for fold_name, fold in body["metrics"].get("folds", {}).items() if isinstance(body["metrics"].get("folds"), dict) else []:
            if args.split == "system_loso" and fold_name not in {"hf", "ti", "zr"}:
                continue
            if isinstance(fold, dict) and isinstance(fold.get("overlap_T"), dict):
                for pair, record in pair_records({"folds": {fold_name: fold["overlap_T"]}}):
                    records.append(("overlap_T", pair, record))
        for subset, pair, record in records:
            references = record.get("reference_Tc_K", [])
            errors = record.get("Tc_error_K", [])
            slopes = record.get("crossing_slope_eV_per_atom_per_K", [])
            n = max(len(references), len(errors), len(slopes), 1)
            reference_value = references if references else (
                "n/a:no_training_phase"
                if record.get("status") == "unavailable_without_training_phase"
                else "n/a:no_reference_crossing"
            )
            if references:
                tc_marker = "n/a:missed_crossing" if (not errors or any(value == "n/a:missed_crossing" for value in errors)) else ""
            elif record.get("status") == "unavailable_without_training_phase":
                tc_marker = "n/a:no_training_phase"
            else:
                tc_marker = "n/a:no_reference_crossing"
            for index in range(n):
                fit = reference_crossing(root, pair, index)
                if ":" in pair:
                    reference_path = root / "data/processed" / pair.split(":", 1)[0] / "reference_crossings.json"
                    if reference_path.exists():
                        reference_inputs.add(reference_path)
                slope_value = slopes[index] if index < len(slopes) else (tc_marker or "")
                if slope_value is None or slope_value == "":
                    slope_value = tc_marker or "n/a:undefined_crossing_slope"
                dg_values = record.get("Tc_err_from_dG_K", [])
                dg_value = dg_values[index] if index < len(dg_values) else tc_marker
                if dg_value is None or dg_value == "":
                    dg_value = tc_marker or "n/a:undefined_crossing_slope"
                if fit:
                    slope_value = fit.get("slope_eV_per_atom_per_K", slope_value)
                    dg_error = record.get("delta_G_MAE_at_crossing_eV_per_atom", [])
                    if not tc_marker and isinstance(dg_error, list) and index < len(dg_error) and isinstance(dg_error[index], (int, float)):
                        dg_value = float(dg_error[index]) / float(slope_value)
                marker_value = tc_marker or "n/a:undefined_crossing_slope"
                window = fit.get("fit_window_K", [marker_value, marker_value]) if fit else [marker_value, marker_value]
                row = {"predictor": predictor, "subset": subset, "fold": record.get("fold") or "all", "pair": pair, "crossing_index": index + 1, "reference_Tc_K": reference_value[index] if isinstance(reference_value, list) and index < len(reference_value) else (reference_value if index == 0 else ""), "Tc_error_K": errors[index] if index < len(errors) else tc_marker, "crossing_slope_eV_per_atom_per_K": slope_value or marker_value, "crossing_slope_stderr_eV_per_atom_per_K": fit.get("slope_stderr", marker_value) if fit else marker_value, "crossing_fit_window_low_K": window[0], "crossing_fit_window_high_K": window[1], "crossing_slope_n_points": fit.get("n_points", marker_value) if fit else marker_value, "crossing_slope_r_squared": fit.get("r_squared", marker_value) if fit else marker_value, "Tc_err_from_dG_K": dg_value, "false_crossings": tc_marker or record.get("false_crossings", ""), "missed_crossings": tc_marker or record.get("missed_crossings", "")}
                (overlap_rows if subset == "overlap_T" else rows).append(row)
            if not references:
                reason = (
                    "predictor unavailable without training phase"
                    if record.get("status") == "unavailable_without_training_phase"
                    else "structural N/A: pair has no reference crossing on the frozen grid"
                )
                missing.append({"field": f"{predictor}:{pair}.reference_Tc_K", "reason": reason})
    inputs.extend(sorted(reference_inputs))
    fields = ["predictor", "subset", "fold", "pair", "crossing_index", "reference_Tc_K", "Tc_error_K", "crossing_slope_eV_per_atom_per_K", "crossing_slope_stderr_eV_per_atom_per_K", "crossing_fit_window_low_K", "crossing_fit_window_high_K", "crossing_slope_n_points", "crossing_slope_r_squared", "Tc_err_from_dG_K", "false_crossings", "missed_crossings"]
    csv_path = output / f"crossing_errors_{args.split}.csv"; csv_write(csv_path, fields, rows); meta_write(root, csv_path, inputs, missing)
    if overlap_rows:
        overlap_path = output / f"crossing_errors_{args.split}_overlap_T.csv"; csv_write(overlap_path, fields, overlap_rows); meta_write(root, overlap_path, inputs, missing, {"subset": "overlap_T"})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
