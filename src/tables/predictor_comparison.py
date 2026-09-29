"""Predictor comparison on each frozen split's test frames — plan Table 5.

Reads normalized raw runs only. Folded splits are scored within each test fold,
then aggregated with pair test-point weights; default rows exclude pairs without
a reference crossing.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from common import csv_write, meta_write, numeric, raw_runs


PAIR_ONLY = "n/a:pair_only_predictor"
NO_TRAINING = "n/a:no_training_phase"
DEGENERATE = "n/a:degenerate_prediction"
NOT_TRAINED = "n/a:not_trained_for_split"
NO_SEED = "n/a:no_seed"


def folds(metrics: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {str(key): value for key, value in metrics.get("folds", {"all": metrics}).items() if isinstance(value, dict)}


def marker(values: list[Any]) -> str | None:
    for value in values:
        if value == DEGENERATE:
            return DEGENERATE
    for value in values:
        if value == NO_TRAINING:
            return NO_TRAINING
    return None


def weighted(values: list[tuple[float, int]]) -> float | None:
    if not values:
        return None
    total = sum(weight for _, weight in values)
    return sum(value * weight for value, weight in values) / total if total else None


def summarize(metrics: dict[str, Any], include_all_pairs: bool, system: str | None = None) -> dict[str, Any]:
    """Aggregate one predictor from its own split/fold test metrics."""

    g_values: list[tuple[float, int]] = []
    ranking_values: list[tuple[float, int]] = []
    pair_values: dict[str, list[tuple[float, int]]] = {key: [] for key in ("delta_G_MAE_eV_per_atom", "delta_G_RMSE_eV_per_atom", "sign_accuracy")}
    tc_values: list[float] = []
    slope_values: list[float] = []
    false_values: list[int] = []
    missed_values: list[int] = []
    markers: list[Any] = []
    pairs_covered: set[str] = set()
    all_pairs_seen: set[str] = set()
    for fold in folds(metrics).values():
        g = numeric(fold.get("G_MAE_eV_per_atom"))
        if g is not None:
            g_values.append((g, int(fold.get("n_test_frames", 0) or 0)))
        ranking_source = fold.get("ranking_accuracy") if system is None else fold.get("ranking_accuracy_by_system", {}).get(system)
        ranking = numeric(ranking_source)
        if ranking is not None:
            ranking_values.append((ranking, int(fold.get("n_test_frames", 0) or 0)))
        for name, pair in fold.get("pairs", {}).items():
            if not isinstance(pair, dict):
                continue
            if system is not None and not str(name).startswith(f"{system}:"):
                continue
            all_pairs_seen.add(str(name))
            if not include_all_pairs and not pair.get("reference_Tc_K", []):
                continue
            if pair.get("status") in {"not_fittable_without_pair_training_labels", "unavailable_without_training_phase"}:
                markers.append(NO_TRAINING)
                continue
            pairs_covered.add(str(name))
            weight = int(pair.get("n_evaluation_points", fold.get("n_test_frames", 0)) or 0)
            for key, values in pair_values.items():
                value = numeric(pair.get(key))
                if value is not None:
                    values.append((value, weight))
            for key, target in (("Tc_error_K", tc_values), ("Tc_err_from_dG_K", slope_values)):
                values = pair.get(key, [])
                for value in values if isinstance(values, list) else [values]:
                    if isinstance(value, str) and value.startswith("n/a:"):
                        markers.append(value)
                    elif numeric(value) is not None:
                        target.append(abs(float(value)))
            for key, target in (("false_crossings", false_values), ("missed_crossings", missed_values)):
                value = pair.get(key)
                if isinstance(value, str) and value.startswith("n/a:"):
                    markers.append(value)
                elif numeric(value) is not None:
                    target.append(int(value))
    unavailable = marker(markers)
    pair_only = str(metrics.get("method", "")) in {"zero", "constant_delta_g", "global_mean_delta_g"}
    g_mae = weighted(g_values) if system is None and g_values else (PAIR_ONLY if pair_only else "n/a:input_unavailable")
    return {
        "G_MAE_eV_per_atom": g_mae,
        "delta_G_MAE_eV_per_atom": weighted(pair_values["delta_G_MAE_eV_per_atom"]) if pair_values["delta_G_MAE_eV_per_atom"] else (unavailable or ""),
        "delta_G_RMSE_eV_per_atom": weighted(pair_values["delta_G_RMSE_eV_per_atom"]) if pair_values["delta_G_RMSE_eV_per_atom"] else (unavailable or ""),
        "sign_accuracy": weighted(pair_values["sign_accuracy"]) if pair_values["sign_accuracy"] else (unavailable or ""),
        "ranking_accuracy": weighted(ranking_values) if ranking_values else (PAIR_ONLY if g_mae == PAIR_ONLY else (unavailable or "n/a:input_unavailable")),
        "Tc_error_K": sum(tc_values) / len(tc_values) if tc_values else (unavailable or ""),
        "Tc_err_from_dG_K": sum(slope_values) / len(slope_values) if slope_values else (unavailable or ""),
        "false_crossings": sum(false_values) if false_values else (unavailable or ""),
        "missed_crossings": sum(missed_values) if missed_values else (unavailable or ""),
        "pairs_covered": len(pairs_covered),
        "all_pairs_seen": len(all_pairs_seen),
    }


def floor(rows: dict[str, dict[str, Any]]) -> tuple[str, float | None]:
    for predictor in ("constant_delta_g", "global_mean_delta_g"):
        value = numeric(rows.get(predictor, {}).get("delta_G_MAE_eV_per_atom"))
        if value is not None and value > 0:
            return predictor, value
    return "", None


def result_rows(rows: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    floor_name, floor_value = floor(rows)
    output = []
    for predictor, summary in rows.items():
        mae = numeric(summary["delta_G_MAE_eV_per_atom"])
        output.append({
            "predictor": predictor,
            **{key: summary[key] for key in ("G_MAE_eV_per_atom", "delta_G_MAE_eV_per_atom", "delta_G_RMSE_eV_per_atom", "sign_accuracy", "ranking_accuracy", "Tc_error_K", "Tc_err_from_dG_K", "false_crossings", "missed_crossings", "pairs_covered")},
            "skill_score": 1.0 - mae / floor_value if mae is not None and floor_value else "",
            "floor_predictor": floor_name if mae is not None and floor_value else "",
            "seed_mean": NO_SEED,
            "seed_std": NO_SEED,
        })
    return output


def e1_rows() -> list[dict[str, Any]]:
    return [{"predictor": name, **{field: NOT_TRAINED for field in FIELDS if field not in {"predictor", "pairs_covered"}}, "pairs_covered": 0} for name in ("repr. regression (polynomial)", "repr. regression (tlog)")]


FIELDS = ["predictor", "G_MAE_eV_per_atom", "delta_G_MAE_eV_per_atom", "delta_G_RMSE_eV_per_atom", "sign_accuracy", "ranking_accuracy", "Tc_error_K", "Tc_err_from_dG_K", "false_crossings", "missed_crossings", "pairs_covered", "skill_score", "floor_predictor", "seed_mean", "seed_std"]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--split", required=True)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--output-root", type=Path, default=Path("result/tables"))
    args = parser.parse_args(); root = args.repo_root.resolve(); output = args.output_root if args.output_root.is_absolute() else root / args.output_root
    inputs: list[Path] = []
    default: dict[str, dict[str, Any]] = {}
    all_pairs: dict[str, dict[str, Any]] = {}
    sio2: dict[str, dict[str, Any]] = {}
    fold_rows: list[dict[str, Any]] = []
    for path, body in raw_runs(root, args.split):
        inputs.append(path)
        predictor = path.parts[-3]
        metrics = body["metrics"]
        default[predictor] = summarize(metrics, include_all_pairs=False)
        all_pairs[predictor] = summarize(metrics, include_all_pairs=True)
        sio2[predictor] = summarize(metrics, include_all_pairs=False, system="sio2")
        for fold_name, fold_body in folds(metrics).items():
            fold_rows.append({"predictor": predictor, "fold": fold_name, "n_test_frames": fold_body.get("n_test_frames", ""), **summarize({"folds": {fold_name: fold_body}}, include_all_pairs=False)})
    rows = result_rows(default) + e1_rows()
    all_rows = result_rows(all_pairs) + e1_rows()
    sio2_rows = result_rows(sio2) + e1_rows()
    rows.sort(key=lambda row: row["predictor"]); all_rows.sort(key=lambda row: row["predictor"]); sio2_rows.sort(key=lambda row: row["predictor"]); fold_rows.sort(key=lambda row: (row["predictor"], row["fold"]))
    csv_path = output / f"predictor_comparison_{args.split}.csv"
    csv_write(csv_path, FIELDS, rows)
    csv_write(output / f"predictor_comparison_{args.split}_include_all_pairs.csv", FIELDS, all_rows)
    csv_write(output / f"predictor_comparison_{args.split}_sio2.csv", FIELDS, sio2_rows)
    csv_write(output / f"predictor_comparison_{args.split}_folds.csv", ["predictor", "fold", "n_test_frames", "G_MAE_eV_per_atom", "delta_G_MAE_eV_per_atom", "delta_G_RMSE_eV_per_atom", "sign_accuracy", "Tc_error_K", "Tc_err_from_dG_K", "false_crossings", "missed_crossings", "pairs_covered", "all_pairs_seen"], fold_rows)
    meta_write(root, csv_path, inputs, [], {"aggregation": {"folded_splits": "pair metrics: per-fold test region then n_evaluation_points-weighted; G MAE: n_test_frames-weighted", "default_pairs": "exclude pairs without a reference crossing", "include_all_pairs_csv": f"predictor_comparison_{args.split}_include_all_pairs.csv", "sio2_csv": f"predictor_comparison_{args.split}_sio2.csv", "fold_detail_csv": f"predictor_comparison_{args.split}_folds.csv"}, "e1": "All E1 checkpoints trained all four SiO2 phases over the full source grid and are therefore ineligible for every frozen split."})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
