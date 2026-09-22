"""Render Table-6-style E1 and calibration-window E2 reports from raw metrics."""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


def number(value: float | None, digits: int = 3) -> str:
    if value is None or not math.isfinite(value):
        return "—"
    return f"{value:.{digits}f}"


def rel(delta: float, baseline: float) -> float | None:
    return None if baseline == 0.0 else 100.0 * delta / abs(baseline)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--e1-output", required=True, type=Path)
    parser.add_argument("--e2-output", required=True, type=Path)
    args = parser.parse_args()
    payload = json.loads(args.input.read_text(encoding="utf-8"))
    records = payload["checkpoint_records"]
    formal_pairs: dict[tuple[str, str], dict[int, dict[str, object]]] = defaultdict(dict)
    for record in records:
        basis, seed = record["basis"], int(record["seed"])
        pairs = record["metrics"]["formal_train_window"]["full_range"]["pairs"]
        for pair, metric in pairs.items():
            formal_pairs[(basis, pair)][seed] = metric

    e1 = [
        "# E1 — SiO2 Tc re-evaluation on the three-phase subset",
        "",
        "Historic four-phase checkpoints are evaluated only on quartz_beta, cristobalite_beta, and tridymite_p63mmc. C2221 is excluded. "
        "Every model was evaluated on the original 851–2499 K, 1 K grid; the predictor-level gauge uses the frozen temp-extrap train window (851–1442 K).",
        "",
        "Each seed cell is `predicted Tc (Tc error; false/missed)`. The reference slope is `|dΔG/dT|`; the slope conversion is reported as the evaluator's crossing-local `|ΔG error|/|dΔG/dT|`.",
        "",
        "## Table 6-style per-seed Tc results",
        "",
        "| predictor | pair | reference Tc (K) | |dΔG/dT| (eV/atom/K) | seed 11 | seed 23 | seed 37 | seed 51 | seed 67 | mean Tc error ± seed std (K) | ΔG→Tc (K) |",
        "|---|---|---:|---:|---|---|---|---|---|---:|---:|",
    ]
    summary_rows = []
    for (basis, pair), by_seed in sorted(formal_pairs.items()):
        first = by_seed[sorted(by_seed)[0]]
        reference = first["reference_Tc_K"]
        roots = max(len(reference), max(len(item["predicted_Tc_K"]) for item in by_seed.values()))
        for root_index in range(roots):
            cells, errors, conversions, false, missed = [], [], [], 0, 0
            for seed in (11, 23, 37, 51, 67):
                item = by_seed.get(seed)
                if item is None:
                    cells.append("missing")
                    continue
                predicted = item["predicted_Tc_K"]
                err = item["Tc_error_K"]
                if root_index < len(predicted):
                    error_value = err[root_index] if root_index < len(err) else None
                    cells.append(f"{number(predicted[root_index], 2)} ({number(error_value, 2)}; {item['false_crossings']}/{item['missed_crossings']})")
                    if root_index < len(err):
                        errors.append(float(err[root_index]))
                    converted = item["Tc_err_from_dG_K"]
                    if root_index < len(converted):
                        conversions.append(float(converted[root_index]))
                else:
                    cells.append(f"— (—; {item['false_crossings']}/{item['missed_crossings']})")
                false += int(item["false_crossings"])
                missed += int(item["missed_crossings"])
            ref = float(reference[root_index]) if root_index < len(reference) else None
            slope = first["crossing_slope_eV_per_atom_per_K"]
            slope_value = float(slope[root_index]) if root_index < len(slope) else None
            e1.append(
                "| " + " | ".join(
                    [basis, f"{pair}#{root_index + 1}", number(ref, 2), number(slope_value, 8), *cells,
                     f"{number(float(np.mean(errors)) if errors else None, 2)} ± {number(float(np.std(errors)) if errors else None, 2)}",
                     number(float(np.mean(conversions)) if conversions else None, 2)]
                ) + " |"
            )
            summary_rows.append({"basis": basis, "pair": pair, "root": root_index + 1, "reference_Tc_K": ref, "mean_Tc_error_K": float(np.mean(errors)) if errors else None, "seed_std_K": float(np.std(errors)) if errors else None, "false_crossings_total": false, "missed_crossings_total": missed})
    e1 += ["", "## Aggregate curve metrics", "", "| predictor | G MAE (eV/atom) | ΔG MAE (eV/atom), mean across pairs/seeds | sign accuracy, mean |", "|---|---:|---:|---:|"]
    for basis in sorted({str(item["basis"]) for item in records}):
        subset = [item for item in records if item["basis"] == basis]
        g_mae = [float(item["metrics"]["formal_train_window"]["full_range"]["G_MAE_eV_per_atom"]) for item in subset]
        pair_values = [float(pair["delta_G_MAE_eV_per_atom"]) for item in subset for pair in item["metrics"]["formal_train_window"]["full_range"]["pairs"].values()]
        sign = [float(pair["sign_accuracy"]) for item in subset for pair in item["metrics"]["formal_train_window"]["full_range"]["pairs"].values()]
        e1.append(f"| {basis} | {np.mean(g_mae):.6g} ± {np.std(g_mae):.2g} | {np.mean(pair_values):.6g} ± {np.std(pair_values):.2g} | {np.mean(sign):.6g} ± {np.std(sign):.2g} |")
    args.e1_output.parent.mkdir(parents=True, exist_ok=True)
    args.e1_output.write_text("\n".join(e1) + "\n", encoding="utf-8")

    e2 = [
        "# E2 — SiO2 calibration-window leakage audit",
        "",
        "Calibration is one least-squares energy-gauge constant for the SiO2 predictor: `mean(G_reference − G_checkpoint)`. No model is retrained. "
        "Formal constants use 851–1442 K (canonical temp-extrap train indexes 0–591); audit constants improperly use 1443–2499 K (held-out indexes 592–1648).",
        "",
        "The frozen historic checkpoint configurations have no `c_system` field. The two constants below are therefore post-hoc predictor-layer audit variants, not a claim about the original training protocol.",
        "",
        "| predictor/seed | c_train (eV/atom) | c_heldout (eV/atom) | held-out G MAE train-c | held-out G MAE heldout-c | ΔMAE | relative ΔMAE | ΔG MAE change | Tc change |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    e2_rows = []
    for record in sorted(records, key=lambda item: (item["basis"], int(item["seed"]))):
        formal = record["metrics"]["formal_train_window"]["temp_extrap"]
        held = record["metrics"]["heldout_window_audit"]["temp_extrap"]
        formal_mae, held_mae = float(formal["G_MAE_eV_per_atom"]), float(held["G_MAE_eV_per_atom"])
        pair_deltas, tc_deltas = [], []
        for name, pair in formal["pairs"].items():
            comparison = held["pairs"][name]
            pair_deltas.append(float(comparison["delta_G_MAE_eV_per_atom"]) - float(pair["delta_G_MAE_eV_per_atom"]))
            tc_deltas.extend([abs(float(a) - float(b)) for a, b in zip(comparison["predicted_Tc_K"], pair["predicted_Tc_K"])])
        delta = held_mae - formal_mae
        e2.append(
            f"| {record['basis']}/seed{record['seed']} | {record['calibration_constants_eV_per_atom']['formal_train_window']:.6g} | "
            f"{record['calibration_constants_eV_per_atom']['heldout_window_audit']:.6g} | {formal_mae:.6g} | {held_mae:.6g} | "
            f"{delta:.6g} | {number(rel(delta, formal_mae), 3)}% | {np.mean(pair_deltas):.3g} | {np.max(tc_deltas) if tc_deltas else 0.0:.3g} |"
        )
        e2_rows.append({"basis": record["basis"], "seed": record["seed"], "formal_G_MAE": formal_mae, "heldout_G_MAE": held_mae, "delta_G_MAE": delta, "relative_delta_percent": rel(delta, formal_mae), "mean_pair_delta_G_MAE": float(np.mean(pair_deltas)), "max_Tc_delta_K": float(np.max(tc_deltas) if tc_deltas else 0.0)})
    e2 += ["", "A common system constant cancels exactly from every ΔG curve. The reported zero ΔG/Tc changes are an expected invariance that the two evaluator runs verify, while the scalar G MAE is sensitive to leakage."]
    args.e2_output.write_text("\n".join(e2) + "\n", encoding="utf-8")
    args.e2_output.with_suffix(".json").write_text(json.dumps({"rows": e2_rows}, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
