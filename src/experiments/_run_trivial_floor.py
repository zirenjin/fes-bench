"""Evaluate label-free zero and train-fitted pairwise-constant ΔG floors.

Unlike phase-energy predictors, a constant ΔG is defined independently for
each unordered phase pair.  It therefore cannot be represented by a single
set of phase curves when three or more phases are present.  This script owns
the pairwise calculation explicitly and never uses a test label to fit the
constant.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np

from fes_bench.data.load import load
from fes_bench.eval.run import _root


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", default="data/processed")
    parser.add_argument("--splits", nargs="+", default=["data/processed/splits/temp_extrap.json", "data/processed/splits/phase_lopo.json", "data/processed/splits/system_loso.json"])
    parser.add_argument("--output-root", default="result/experiments/skill_floor")
    return parser


def _rows_by_phase(rows: list[dict[str, object]]) -> dict[tuple[str, str], set[int]]:
    grouped: dict[tuple[str, str], set[int]] = defaultdict(set)
    for row in rows:
        grouped[(str(row["system"]), str(row["phase"]))].add(int(row["T_index"]))
    return grouped


def _curve(data_root: Path, system: str, phase: str, indexes: set[int]) -> dict[float, float]:
    table = load(system, phase, data_root).G_table
    return {float(table[i].T_K): float(table[i].G_eV_per_atom) for i in indexes}


def _slope(reference: np.ndarray, temperatures: np.ndarray, root: float) -> float:
    nearest = int(np.argmin(np.abs(temperatures - root)))
    if nearest == 0:
        return float(abs((reference[1] - reference[0]) / (temperatures[1] - temperatures[0])))
    if nearest == len(reference) - 1:
        return float(abs((reference[-1] - reference[-2]) / (temperatures[-1] - temperatures[-2])))
    return float(abs((reference[nearest + 1] - reference[nearest - 1]) / (temperatures[nearest + 1] - temperatures[nearest - 1])))


def _pair_metrics(reference: np.ndarray, observed: np.ndarray, temperatures: np.ndarray, *, fitted: bool, constant: float | None) -> dict[str, object]:
    roots_ref = _root(reference, temperatures)
    if not fitted:
        return {
            "status": "not_fittable_without_pair_training_labels",
            "train_delta_G_mean_eV_per_atom": constant,
            "n_evaluation_points": int(len(temperatures)),
            "delta_G_MAE_eV_per_atom": "n/a:no_training_phase",
            "delta_G_RMSE_eV_per_atom": "n/a:no_training_phase",
            "sign_accuracy": "n/a:no_training_phase",
            "reference_Tc_K": roots_ref,
            "predicted_Tc_K": ["n/a:no_training_phase"],
            "Tc_error_K": ["n/a:no_training_phase"] * len(roots_ref),
            "delta_G_MAE_at_crossing_eV_per_atom": ["n/a:no_training_phase"] * len(roots_ref),
            "crossing_slope_eV_per_atom_per_K": ["n/a:no_training_phase"] * len(roots_ref),
            "Tc_err_from_dG_K": ["n/a:no_training_phase"] * len(roots_ref),
            "false_crossings": "n/a:no_training_phase",
            "missed_crossings": "n/a:no_training_phase",
        }
    roots_pred = _root(observed, temperatures)
    degenerate = bool(np.all(np.abs(observed) <= 1.0e-12))
    crossing_error = [float(abs(np.interp(root, temperatures, observed - reference))) for root in roots_ref]
    slopes = [_slope(reference, temperatures, root) for root in roots_ref]
    return {
        "status": "degenerate_prediction" if degenerate else "ok",
        "train_delta_G_mean_eV_per_atom": constant,
        "n_evaluation_points": int(len(temperatures)),
        "delta_G_MAE_eV_per_atom": float(np.mean(np.abs(observed - reference))),
        "delta_G_RMSE_eV_per_atom": float(np.sqrt(np.mean((observed - reference) ** 2))),
        "sign_accuracy": float(np.mean(np.sign(observed) == np.sign(reference))),
        "reference_Tc_K": roots_ref,
        "predicted_Tc_K": roots_pred if not degenerate else ["n/a:degenerate_prediction"],
        "Tc_error_K": [float(pred - ref) for pred, ref in zip(roots_pred, roots_ref)] if not degenerate else ["n/a:degenerate_prediction"] * len(roots_ref),
        "delta_G_MAE_at_crossing_eV_per_atom": crossing_error if not degenerate else ["n/a:degenerate_prediction"] * len(roots_ref),
        "crossing_slope_eV_per_atom_per_K": slopes if not degenerate else ["n/a:degenerate_prediction"] * len(roots_ref),
        "Tc_err_from_dG_K": [float(error / slope) if slope else math.nan for error, slope in zip(crossing_error, slopes)] if not degenerate else ["n/a:degenerate_prediction"] * len(roots_ref),
        "false_crossings": max(0, len(roots_pred) - len(roots_ref)) if not degenerate else "n/a:degenerate_prediction",
        "missed_crossings": max(0, len(roots_ref) - len(roots_pred)) if not degenerate else "n/a:degenerate_prediction",
    }


def _evaluate_fold(fold: dict[str, object], data_root: Path, method: str) -> dict[str, object]:
    train = _rows_by_phase(list(fold.get("train", [])))
    test = _rows_by_phase(list(fold["test"]))
    all_rows = {key: set(test.get(key, set())) | set(train.get(key, set())) for key in set(train) | set(test)}
    # Cache source tables once per phase.  The ranking audit visits every
    # temperature and otherwise repeatedly reparses the same CSVs.
    tables = {(system, phase): load(system, phase, data_root).G_table for system, phase in all_rows}

    def cached_curve(system: str, phase: str, indexes: set[int]) -> dict[float, float]:
        table = tables[(system, phase)]
        return {float(table[i].T_K): float(table[i].G_eV_per_atom) for i in indexes}

    pairs: dict[str, object] = {}
    ranking_values: list[bool] = []
    ranking_by_system: dict[str, list[bool]] = defaultdict(list)
    ranking_unavailable = False
    for system in sorted({key[0] for key in test}):
        phases = sorted(phase for candidate_system, phase in all_rows if candidate_system == system)
        for left, right in combinations(phases, 2):
            if (system, left) not in test and (system, right) not in test:
                continue
            # When both phases contribute test rows (temp extrapolation), score
            # only their shared test grid.  In LOPO/LOSO one partner is train
            # only, so retain the full shared grid needed for a physical pair.
            test_shared = test.get((system, left), set()) & test.get((system, right), set())
            eval_shared = test_shared or (all_rows[(system, left)] & all_rows[(system, right)])
            if not eval_shared:
                continue
            train_shared = train.get((system, left), set()) & train.get((system, right), set())
            left_eval, right_eval = cached_curve(system, left, eval_shared), cached_curve(system, right, eval_shared)
            temperatures = np.array(sorted(set(left_eval) & set(right_eval)))
            reference = np.array([left_eval[t] - right_eval[t] for t in temperatures])
            if method == "zero":
                observed, fitted, constant = np.zeros_like(reference), True, 0.0
            elif train_shared:
                left_train, right_train = cached_curve(system, left, train_shared), cached_curve(system, right, train_shared)
                shared_train_t = sorted(set(left_train) & set(right_train))
                constant = float(np.mean([left_train[t] - right_train[t] for t in shared_train_t]))
                observed, fitted = np.full_like(reference, constant), True
            else:
                # Keep the diagnostic zero curve and mark it unavailable; it
                # must not enter a fitted-constant aggregate.
                observed, fitted, constant = np.zeros_like(reference), False, None
            pairs[f"{system}:{left}_minus_{right}"] = _pair_metrics(reference, observed, temperatures, fitted=fitted, constant=constant)
        evaluation_temperatures = sorted({float(tables[(system, phase)][index].T_K) for (item_system, phase), indexes in test.items() if item_system == system for index in indexes})
        phase_curves = {phase: cached_curve(system, phase, all_rows[(system, phase)]) for phase in phases}
        for temperature in evaluation_temperatures:
            ok = True
            for left, right in combinations(phases, 2):
                if temperature not in phase_curves[left] or temperature not in phase_curves[right]:
                    ok = False; break
                reference = phase_curves[left][temperature] - phase_curves[right][temperature]
                if method == "zero":
                    observed = 0.0
                else:
                    train_shared = train.get((system, left), set()) & train.get((system, right), set())
                    if not train_shared:
                        ranking_unavailable = True
                        ok = False; break
                    left_train, right_train = cached_curve(system, left, train_shared), cached_curve(system, right, train_shared)
                    constant_value = float(np.mean([left_train[t] - right_train[t] for t in sorted(set(left_train) & set(right_train))]))
                    observed = constant_value
                if np.sign(observed) != np.sign(reference):
                    ok = False; break
            ranking_values.append(ok)
            ranking_by_system[system].append(ok)
    eligible = [v for v in pairs.values() if isinstance(v, dict) and v["status"] == "ok"]
    return {
        "n_test_frames": sum(len(indexes) for indexes in test.values()),
        "G_MAE_eV_per_atom": None,
        "G_MAE_note": "Not defined: both floors are ΔG-pair predictors, not phase-energy predictors.",
        "pairs": pairs,
        "aggregate_pair_MAE_eV_per_atom": float(np.mean([v["delta_G_MAE_eV_per_atom"] for v in eligible])) if eligible else None,
        "aggregate_pair_sign_accuracy": float(np.mean([v["sign_accuracy"] for v in eligible])) if eligible else None,
        "n_fittable_pairs": len(eligible),
        "ranking_accuracy": "n/a:no_training_phase" if ranking_unavailable else (float(np.mean(ranking_values)) if ranking_values else None),
        "ranking_accuracy_by_system": {system: ("n/a:no_training_phase" if ranking_unavailable else float(np.mean(values))) for system, values in ranking_by_system.items() if values},
    }


def _evaluate(split: dict[str, object], data_root: Path, method: str) -> dict[str, object]:
    if isinstance(split.get("folds"), dict):
        folds = {name: _evaluate_fold(fold, data_root, method) for name, fold in split["folds"].items()}
        return {"method": method, "folds": folds}
    return {"method": method, "folds": {"all": _evaluate_fold(split, data_root, method)}}


def _fmt(value: object) -> str:
    if value is None:
        return "—"
    if isinstance(value, str) and value.startswith("n/a:"):
        return value
    return f"{float(value):.6g}"


def _markdown(payload: dict[str, object], title: str) -> str:
    lines = [f"# {title}", "", "Skill scores in the tables use the zero floor uniformly: `1 − MAE / MAE_zero`, with `MAE_zero` equal to the mean reference |ΔG| on evaluated points. `constant_delta_g` is fitted separately for each ordered ΔG = G(left) − G(right) from pair *training* points and is unavailable when a fold has no pairwise training labels.", "", "Phase-level G MAE and coverage are not defined for these pair predictors.", ""]
    for split_name, result in payload.items():
        lines.extend([f"## {split_name}", ""])
        for fold_name, fold in result["folds"].items():
            lines.extend([f"### {fold_name}", "", f"Fittable pairs: {fold['n_fittable_pairs']}; aggregate ΔG MAE: {_fmt(fold['aggregate_pair_MAE_eV_per_atom'])} eV/atom; aggregate sign accuracy: {_fmt(fold['aggregate_pair_sign_accuracy'])}.", "", "| Pair | status | train mean ΔG (eV/atom) | ΔG MAE | ΔG RMSE | sign accuracy | reference / predicted Tc (K) | false / missed |", "|---|---|---:|---:|---:|---:|---|---:|"])
            for name, pair in fold["pairs"].items():
                lines.append(f"| {name} | {pair['status']} | {_fmt(pair['train_delta_G_mean_eV_per_atom'])} | {_fmt(pair['delta_G_MAE_eV_per_atom'])} | {_fmt(pair['delta_G_RMSE_eV_per_atom'])} | {_fmt(pair['sign_accuracy'])} | {', '.join(_fmt(x) for x in pair['reference_Tc_K']) or '—'} / {', '.join(_fmt(x) for x in pair['predicted_Tc_K']) or '—'} | {pair['false_crossings']} / {pair['missed_crossings']} |")
            lines.append("")
    lines.extend(["## Deviations from design", "", "No new model was trained. The design evaluator's phase-level G MAE/coverage cannot be applied to a pair-only ΔG predictor without inventing a phase-energy gauge, so they are explicitly N/A. In phase-LOPO and system-LOSO, a held phase/system has no pairwise training ΔG labels; `constant_delta_g` is consequently reported as unavailable rather than fitted with held-out reference values.", ""])
    return "\n".join(lines)


def _json_safe(value: object) -> object:
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_json_safe(item) for item in value]
    return value


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    data_root, output = Path(args.data_root).resolve(), Path(args.output_root).resolve()
    output.mkdir(parents=True, exist_ok=True)
    all_payload: dict[str, object] = {}
    for split_arg in args.splits:
        split_path = Path(split_arg)
        split = json.loads(split_path.read_text(encoding="utf-8"))
        split_name = split_path.stem
        all_payload[split_name] = {method: _evaluate(split, data_root, method) for method in ("zero", "constant_delta_g")}
    (output / "metrics.json").write_text(json.dumps(_json_safe(all_payload), indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    (output / "summary.md").write_text(_markdown({name: data["zero"] for name, data in all_payload.items()}, "Zero ΔG floor") + "\n\n" + _markdown({name: data["constant_delta_g"] for name, data in all_payload.items()}, "Training-mean constant ΔG floor") + "\n", encoding="utf-8")
    print(output / "metrics.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
