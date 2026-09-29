"""Re-evaluate frozen crossing records — plan experiment E1."""

from __future__ import annotations

import argparse
import csv
import json
import statistics
from pathlib import Path
from typing import Any

from _common import write_audit


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--repo-root", type=Path, default=Path(".")); args = parser.parse_args(); root = args.repo_root.resolve()
    rows, inputs = [], []
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for path in sorted((root / "result/experiments/crossing_reevaluation/raw_runs").glob("*/seed_*/metrics.json")):
        inputs.append(path); record = json.loads(path.read_text(encoding="utf-8"))["metrics"]; basis = path.parts[-3]; seed = path.parts[-2].removeprefix("seed_")
        metrics = record.get("metrics", {}).get("formal_train_window", {}).get("full_range", {})
        for pair, item in metrics.get("pairs", {}).items():
            common = {
                "basis": basis,
                "seed": seed,
                "pair": pair,
                "false_crossings": int(item.get("false_crossings", 0)),
                "missed_crossings": int(item.get("missed_crossings", 0)),
                "sign_accuracy": item.get("sign_accuracy"),
            }
            errors = item.get("Tc_error_K", [])
            references = item.get("reference_Tc_K", [])
            inferred = item.get("Tc_err_from_dG_K", [])
            if references:
                for index, reference in enumerate(references):
                    row = {**common, "crossing_index": index + 1, "reference_Tc_K": reference, "Tc_error_K": errors[index] if index < len(errors) else None, "Tc_err_from_dG_K": inferred[index] if index < len(inferred) else None}
                    rows.append(row); grouped.setdefault((basis, pair), []).append(row)
            else:
                row = {**common, "crossing_index": "", "reference_Tc_K": "", "Tc_error_K": "n/a:no_reference_crossing", "Tc_err_from_dG_K": "n/a:no_reference_crossing"}
                rows.append(row); grouped.setdefault((basis, pair), []).append(row)
    summary_rows: list[dict[str, Any]] = []
    conclusion: dict[str, str] = {}
    for (basis, pair), pair_rows in sorted(grouped.items()):
        numeric_errors = [float(row["Tc_error_K"]) for row in pair_rows if isinstance(row.get("Tc_error_K"), (float, int))]
        signs = [float(row["sign_accuracy"]) for row in pair_rows if isinstance(row.get("sign_accuracy"), (float, int))]
        summary_rows.append({
            "basis": basis,
            "pair": pair,
            "n_seeds": len({str(row["seed"]) for row in pair_rows}),
            "Tc_error_mean_K": statistics.mean(numeric_errors) if numeric_errors else "n/a:no_reference_crossing",
            "Tc_error_std_K": statistics.stdev(numeric_errors) if len(numeric_errors) > 1 else (0.0 if numeric_errors else "n/a:no_reference_crossing"),
            "false_crossings_total": sum(int(row["false_crossings"]) for row in pair_rows),
            "missed_crossings_total": sum(int(row["missed_crossings"]) for row in pair_rows),
            "sign_accuracy_mean": statistics.mean(signs) if signs else "",
            "sign_accuracy_std": statistics.stdev(signs) if len(signs) > 1 else (0.0 if signs else ""),
            "status": "reference_crossing" if numeric_errors else "n/a:no_reference_crossing",
        })
    for basis in sorted({row["basis"] for row in summary_rows}):
        crossing = [row for row in summary_rows if row["basis"] == basis and row["status"] == "reference_crossing"]
        no_reference = [row for row in summary_rows if row["basis"] == basis and row["status"] != "reference_crossing"]
        crossing_text = "; ".join(f"{row['pair']}: {float(row['Tc_error_mean_K']):.3f} ± {float(row['Tc_error_std_K']):.3f} K" for row in crossing)
        false_total = sum(int(row["false_crossings_total"]) for row in no_reference)
        conclusion[basis] = f"Two reference-crossing pairs: {crossing_text}. The no-reference pair has {false_total}/5 false crossings for this basis."
    output = root / "result/experiments/crossing_reevaluation"
    with (output / "summary.csv").open("w", newline="", encoding="utf-8") as handle:
        fields = ["basis", "pair", "n_seeds", "Tc_error_mean_K", "Tc_error_std_K", "false_crossings_total", "missed_crossings_total", "sign_accuracy_mean", "sign_accuracy_std", "status"]
        writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader(); writer.writerows(summary_rows)
    write_audit(root, "crossing_reevaluation", ["basis", "seed", "pair", "crossing_index", "reference_Tc_K", "Tc_error_K", "Tc_err_from_dG_K", "false_crossings", "missed_crossings", "sign_accuracy"], rows, inputs, [], extra={"notes": {"acceptance_wording": "The canonical E1 subset has two reference-crossing pairs, not three, because C2221 is excluded.", "no_reference_pair": "cristobalite_beta_minus_tridymite_p63mmc has no reference crossing; Tc fields use n/a:no_reference_crossing.", "acceptance_conclusion": "Across polynomial and tlog_polynomial, 10/10 checkpoint records produce one false crossing for the no-reference cristobalite–tridymite pair; acceptance criterion three is not met.", "conclusion": conclusion}, "supporting_outputs": ["result/experiments/crossing_reevaluation/summary.csv"]})
    return 0


if __name__ == "__main__": raise SystemExit(main())
