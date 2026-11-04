"""Evaluate a predictor specification against a frozen benchmark split.

Metric definitions:
- G_MAE | mean absolute error of phase G | eV/atom | ↓
- delta_G_MAE | mean absolute error of pair ΔG = G(left) − G(right) | eV/atom | ↓
- delta_G_RMSE | root mean squared error of pair ΔG | eV/atom | ↓
- sign_accuracy | fraction of evaluation points with matching ΔG sign | fraction | ↑
- balanced_sign_accuracy | mean of positive- and negative-reference sign recall (classes present on the evaluated grid) | fraction | ↑
- Tc_error | predicted crossing temperature minus reference crossing temperature (signed; smaller magnitude is better) | K | ↓ (magnitude)
- false_crossings | predicted crossings beyond the reference crossing count | count | ↓
- missed_crossings | reference crossings not predicted | count | ↓
- coverage_2sigma | fraction of reference values covered by seed mean ± 2σ | fraction | ↑
- ranking_accuracy | fraction of evaluation points at which every phase pair of a system has the correct ΔG sign, i.e. the full phase ordering is correct | fraction | ↑
- Tc_err_from_dG | ΔG error at the reference crossing divided by |dΔG/dT| there; converts an energy error into a temperature error | K | ↓
- skill_score | 1 − MAE / MAE_zero, where MAE_zero is the mean of reference |ΔG| over the evaluated points (the MAE of predicting ΔG = 0); negative means worse than predicting zero | fraction | ↑
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

from fes_bench.data.load import ReferencePoint, load


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--predictor",
        required=True,
        help=(
            "reference, constant_offset:<eV>, reference_iid_noise:<eV>, or "
            "precomputed:<curve-json>"
        ),
    )
    parser.add_argument("--split", required=True, help="frozen split JSON")
    parser.add_argument("--data-root", required=True)
    parser.add_argument("--seeds", default="11,23,37,51,67")
    parser.add_argument("--output-root", default="result/experiments/external_baselines")
    return parser


def _root(delta: np.ndarray, temperatures: np.ndarray, merge_within_K: float = 0.0) -> list[float]:
    roots: list[float] = []
    index = 0
    while index < len(delta):
        left = delta[index]
        if left == 0:
            end = index
            while end + 1 < len(delta) and delta[end + 1] == 0:
                end += 1
            roots.append(float((temperatures[index] + temperatures[end]) / 2.0))
            index = end + 1
            continue
        if index + 1 < len(delta) and left * delta[index + 1] < 0:
            right = delta[index + 1]
            fraction = abs(left) / (abs(left) + abs(right))
            roots.append(float(temperatures[index] + fraction * (temperatures[index + 1] - temperatures[index])))
        index += 1
    if merge_within_K <= 0 or len(roots) < 2:
        return roots
    clusters: list[list[float]] = [[roots[0]]]
    for value in roots[1:]:
        if value - clusters[-1][-1] <= merge_within_K:
            clusters[-1].append(value)
        else:
            clusters.append([value])
    return [float(np.mean(cluster)) for cluster in clusters]


def _balanced_sign_accuracy(observed: np.ndarray, reference: np.ndarray) -> float:
    """Return class-balanced sign accuracy without counting exact-zero labels.

    A non-crossing pair has only one sign class on its frozen grid.  In that
    structural case the mean is taken over the class that is present; this
    keeps the metric defined while making the convention explicit in the
    provenance and table documentation.
    """
    classes: list[float] = []
    positive = reference > 0.0
    negative = reference < 0.0
    if np.any(positive):
        classes.append(float(np.mean(observed[positive] > 0.0)))
    if np.any(negative):
        classes.append(float(np.mean(observed[negative] < 0.0)))
    return float(np.mean(classes)) if classes else 0.0


def _predictor(spec: str, seed: int, key: str, reference: np.ndarray) -> np.ndarray:
    if spec == "reference":
        return reference.copy()
    prefix = "reference_noise:"
    if spec.startswith(prefix):
        sigma = float(spec[len(prefix) :])
        digest = hashlib.sha256(f"{seed}:{key}".encode()).digest()
        rng = np.random.default_rng(int.from_bytes(digest[:8], "little"))
        # One Gaussian energy-gauge perturbation per phase/seed preserves a
        # physically smooth G(T) curve, so Tc displacement can be compared to
        # MAE divided by the local crossing slope.
        return reference + rng.normal(0.0, sigma)
    prefix = "constant_offset:"
    if spec.startswith(prefix):
        return reference + float(spec[len(prefix) :])
    prefix = "reference_iid_noise:"
    if spec.startswith(prefix):
        sigma = float(spec[len(prefix) :])
        digest = hashlib.sha256(f"iid:{seed}:{key}".encode()).digest()
        rng = np.random.default_rng(int.from_bytes(digest[:8], "little"))
        return reference + rng.normal(0.0, sigma, size=reference.shape)
    raise ValueError(
        "predictor must be reference, constant_offset:<eV>, "
        "reference_noise:<eV>, reference_iid_noise:<eV>, or precomputed:<curve-json>"
    )


def _output_slug(spec: str) -> str:
    if spec.startswith("reference_iid_noise:"):
        return "root_stability_iid_" + spec.split(":", 1)[1]
    return spec.replace(":", "_").replace("/", "_")


def _load_precomputed(spec: str) -> dict[str, dict[float, float]] | None:
    """Load per-phase predictor curves for checkpoint reanalysis.

    The JSON payload is deliberately minimal and model-agnostic::

        {"predictions": {"sio2:quartz_beta": {"T_K": [...], "G_eV_per_atom": [...]}}}

    This lets the unified evaluator own all split, root, slope, and crossing
    bookkeeping while an external inference adapter supplies a frozen model's
    values.  Temperatures are lookup keys rather than array positions so a
    caller cannot accidentally apply a checkpoint curve to a mismatched grid.
    """

    prefix = "precomputed:"
    if not spec.startswith(prefix):
        return None
    path = Path(spec[len(prefix) :]).expanduser()
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        raw_predictions = payload["predictions"]
    except (OSError, KeyError, TypeError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid precomputed predictor {path}: {exc}") from exc
    if not isinstance(raw_predictions, dict):
        raise ValueError("precomputed predictor requires an object named predictions")
    curves: dict[str, dict[float, float]] = {}
    for key, raw_curve in raw_predictions.items():
        if not isinstance(key, str) or not isinstance(raw_curve, dict):
            raise ValueError("precomputed predictor entries must be keyed curve objects")
        temperatures = raw_curve.get("T_K")
        values = raw_curve.get("G_eV_per_atom")
        if not isinstance(temperatures, list) or not isinstance(values, list) or len(temperatures) != len(values):
            raise ValueError(f"precomputed curve {key!r} requires equally sized T_K/G_eV_per_atom lists")
        curve: dict[float, float] = {}
        for temperature, value in zip(temperatures, values):
            try:
                t, g = float(temperature), float(value)
            except (TypeError, ValueError) as exc:
                raise ValueError(f"precomputed curve {key!r} has a non-numeric value") from exc
            if not math.isfinite(t) or not math.isfinite(g) or t in curve:
                raise ValueError(f"precomputed curve {key!r} has invalid or duplicate temperature {t!r}")
            curve[t] = g
        curves[key] = curve
    return curves


def _table_cache(data_root: Path, rows: list[dict[str, object]]) -> dict[tuple[str, str], tuple[ReferencePoint, ...]]:
    keys = {(str(row["system"]), str(row["phase"])) for row in rows}
    return {key: load(key[0], key[1], data_root).G_table for key in keys}


def _state_key(point: ReferencePoint) -> tuple[float, float]:
    """Return a stable (T, P) key for exact reference-state pairing."""
    return round(point.T_K, 8), round(point.P_GPa, 8)


def _load_precomputed_pt(spec: str) -> dict[str, dict[tuple[float, float | None], float]] | None:
    """Read prediction curves keyed by both temperature and pressure.

    Multi-pressure data require an explicit ``P_GPa`` list. Legacy
    temperature-only curves remain valid only when the input has a unique
    pressure for every temperature; that check is performed during lookup.
    """
    prefix = "precomputed:"
    if not spec.startswith(prefix):
        return None
    path = Path(spec[len(prefix):]).expanduser()
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        raw_predictions = payload["predictions"]
    except (OSError, KeyError, TypeError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid precomputed predictor {path}: {exc}") from exc
    if not isinstance(raw_predictions, dict):
        raise ValueError("precomputed predictor requires an object named predictions")
    curves: dict[str, dict[tuple[float, float | None], float]] = {}
    for key, raw_curve in raw_predictions.items():
        if not isinstance(key, str) or not isinstance(raw_curve, dict):
            raise ValueError("precomputed predictor entries must be keyed curve objects")
        temperatures = raw_curve.get("T_K")
        values = raw_curve.get("G_eV_per_atom")
        pressures = raw_curve.get("P_GPa")
        if not isinstance(temperatures, list) or not isinstance(values, list) or len(temperatures) != len(values):
            raise ValueError(f"precomputed curve {key!r} requires equally sized T_K/G_eV_per_atom lists")
        if pressures is not None and (not isinstance(pressures, list) or len(pressures) != len(temperatures)):
            raise ValueError(f"precomputed curve {key!r} P_GPa must match T_K length")
        curve: dict[tuple[float, float | None], float] = {}
        for index, (temperature, value) in enumerate(zip(temperatures, values)):
            try:
                t, g = float(temperature), float(value)
                p = float(pressures[index]) if pressures is not None else math.nan
            except (TypeError, ValueError) as exc:
                raise ValueError(f"precomputed curve {key!r} has a non-numeric value") from exc
            if not math.isfinite(t) or not math.isfinite(g) or (pressures is not None and not math.isfinite(p)):
                raise ValueError(f"precomputed curve {key!r} has a non-finite state or value")
            state = (round(t, 8), round(p, 8)) if pressures is not None else (round(t, 8), None)
            if state in curve:
                raise ValueError(f"precomputed curve {key!r} has a duplicate state {state!r}")
            curve[state] = g
        curves[key] = curve
    return curves


def _evaluate_pressure_grid(
    spec: str,
    test: list[dict[str, object]],
    train: list[dict[str, object]],
    cache: dict[tuple[str, str], tuple[ReferencePoint, ...]],
    seeds: list[int],
) -> dict[str, object]:
    """Evaluate a multi-pressure table without pairing unlike P states."""
    precomputed = _load_precomputed_pt(spec)
    pair_rows = test + train
    by_phase: dict[tuple[str, str], set[int]] = defaultdict(set)
    test_by_phase: dict[tuple[str, str], set[int]] = defaultdict(set)
    for row in pair_rows:
        by_phase[(str(row["system"]), str(row["phase"]))].add(int(row["T_index"]))
    for row in test:
        test_by_phase[(str(row["system"]), str(row["phase"]))].add(int(row["T_index"]))

    point_errors: list[float] = []
    coverage_values: list[bool] = []
    ranking_values: list[bool] = []
    pair_metrics: dict[str, dict[str, object]] = {}
    systems = sorted({system for system, _ in test_by_phase})
    for system in systems:
        phases = sorted(phase for item_system, phase in by_phase if item_system == system)
        predictions: dict[str, tuple[list[tuple[float, float]], np.ndarray, np.ndarray]] = {}
        references: dict[str, dict[tuple[float, float], float]] = {}
        for phase in phases:
            table = cache[(system, phase)]
            indexes = sorted(by_phase[(system, phase)])
            states = [_state_key(table[index]) for index in indexes]
            if len(set(states)) != len(states):
                raise ValueError(f"duplicate (T_K, P_GPa) reference states for {system}:{phase}")
            reference = np.asarray([table[index].G_eV_per_atom for index in indexes], dtype=float)
            references[phase] = dict(zip(states, reference.tolist()))
            if precomputed is not None:
                curve_key = f"{system}:{phase}"
                if curve_key not in precomputed:
                    raise ValueError(f"precomputed predictor has no curve for {curve_key!r}")
                curve = precomputed[curve_key]
                missing = []
                values = []
                for t, p in states:
                    key = (t, p)
                    if key not in curve and (t, None) in curve:
                        matches = [g for (curve_t, _), g in curve.items() if curve_t == t]
                        if len(matches) == 1:
                            values.append(matches[0])
                            continue
                        missing.append((t, p))
                        continue
                    if key not in curve:
                        missing.append((t, p))
                    else:
                        values.append(curve[key])
                if missing:
                    raise ValueError(f"precomputed curve {curve_key!r} misses or ambiguously keys states {missing[:3]!r}")
                mean = np.asarray(values, dtype=float)
                seed_predictions = np.stack([mean.copy() for _ in seeds])
            elif spec == "reference":
                mean = reference.copy()
                seed_predictions = np.stack([mean.copy() for _ in seeds])
            else:
                seed_predictions = np.stack(
                    [_predictor(spec, seed, f"{system}:{phase}", reference) for seed in seeds]
                )
                mean = seed_predictions.mean(axis=0)
            predictions[phase] = (states, mean, seed_predictions)
            test_indexes = sorted(test_by_phase.get((system, phase), set()))
            if test_indexes:
                test_states = [_state_key(table[index]) for index in test_indexes]
                prediction_map = dict(zip(states, mean.tolist()))
                test_reference = np.asarray([table[index].G_eV_per_atom for index in test_indexes], dtype=float)
                test_mean = np.asarray([prediction_map[state] for state in test_states], dtype=float)
                point_errors.extend(np.abs(test_mean - test_reference).tolist())
                std_map = dict(zip(states, seed_predictions.std(axis=0).tolist()))
                coverage_values.extend(
                    np.abs(test_mean - test_reference) <= 2.0 * np.asarray([std_map[state] for state in test_states]) + 1.0e-14
                )

        # Ranking is scored independently at each shared (T, P) state.
        evaluation_states = sorted({
            _state_key(cache[(system, phase)][index])
            for (item_system, phase), indexes in test_by_phase.items() if item_system == system
            for index in indexes
        })
        for state in evaluation_states:
            pairs = [(left, right) for i, left in enumerate(phases) for right in phases[i + 1:]
                     if state in references[left] and state in references[right]]
            if not pairs:
                continue
            correct = True
            for left, right in pairs:
                left_states, left_mean, _ = predictions[left]
                right_states, right_mean, _ = predictions[right]
                left_value = dict(zip(left_states, left_mean.tolist()))[state]
                right_value = dict(zip(right_states, right_mean.tolist()))[state]
                if np.sign(left_value - right_value) != np.sign(references[left][state] - references[right][state]):
                    correct = False
                    break
            ranking_values.append(correct)

        for left_index, left in enumerate(phases):
            for right in phases[left_index + 1:]:
                if (system, left) not in test_by_phase and (system, right) not in test_by_phase:
                    continue
                common_states = set(references[left]) & set(references[right])
                test_shared = {
                    _state_key(cache[(system, left)][index])
                    for index in test_by_phase.get((system, left), set())
                } & {
                    _state_key(cache[(system, right)][index])
                    for index in test_by_phase.get((system, right), set())
                }
                if test_shared:
                    common_states &= test_shared
                pressures = sorted({p for _, p in common_states})
                left_states, left_g, left_seed_g = predictions[left]
                right_states, right_g, right_seed_g = predictions[right]
                left_map = dict(zip(left_states, left_g.tolist()))
                right_map = dict(zip(right_states, right_g.tolist()))
                left_seed_maps = [dict(zip(left_states, row.tolist())) for row in left_seed_g]
                right_seed_maps = [dict(zip(right_states, row.tolist())) for row in right_seed_g]
                for pressure in pressures:
                    temperatures = np.asarray(sorted(t for t, p in common_states if p == pressure), dtype=float)
                    states = [(float(t), pressure) for t in temperatures]
                    if not states:
                        continue
                    observed = np.asarray([left_map[state] - right_map[state] for state in states])
                    reference = np.asarray([references[left][state] - references[right][state] for state in states])
                    roots_reference = _root(reference, temperatures)
                    merge_within = 3.0 if np.any(reference == 0) else 0.0
                    roots_predicted = _root(observed, temperatures, merge_within_K=merge_within)
                    seed_roots = [
                        _root(np.asarray([lm[state] - rm[state] for state in states]), temperatures,
                              merge_within_K=merge_within)
                        for lm, rm in zip(left_seed_maps, right_seed_maps)
                    ]
                    crossing_mae: list[float] = []
                    slopes: list[float] = []
                    for root in roots_reference:
                        crossing_mae.append(abs(float(np.interp(root, temperatures, observed - reference))))
                        nearest = int(np.argmin(np.abs(temperatures - root)))
                        if reference[nearest] == 0:
                            lower, upper = nearest, nearest
                            while lower > 0 and reference[lower] == 0:
                                lower -= 1
                            while upper < len(reference) - 1 and reference[upper] == 0:
                                upper += 1
                            slope = ((reference[upper] - reference[lower]) / (temperatures[upper] - temperatures[lower])
                                     if upper != lower else math.nan)
                        elif nearest == 0:
                            slope = (reference[1] - reference[0]) / (temperatures[1] - temperatures[0])
                        elif nearest == len(temperatures) - 1:
                            slope = (reference[-1] - reference[-2]) / (temperatures[-1] - temperatures[-2])
                        else:
                            slope = (reference[nearest + 1] - reference[nearest - 1]) / (temperatures[nearest + 1] - temperatures[nearest - 1])
                        slopes.append(float(abs(slope)))
                    degenerate = bool(np.all(np.abs(observed) <= 1.0e-12))
                    missed = (not degenerate) and len(roots_predicted) < len(roots_reference)
                    key = f"{system}:{left}_minus_{right}"
                    if len(pressures) > 1:
                        key += f"@P={pressure:g}GPa"
                    pair_metrics[key] = {
                        "pressure_GPa": pressure,
                        "pair_support": "test_and_train_partner" if train else "test_only",
                        "status": "degenerate_prediction" if degenerate else "ok",
                        "delta_G_MAE_eV_per_atom": float(np.mean(np.abs(observed - reference))),
                        "delta_G_RMSE_eV_per_atom": float(np.sqrt(np.mean((observed - reference) ** 2))),
                        "delta_G_MAE_at_crossing_eV_per_atom": crossing_mae if not degenerate else ["n/a:degenerate_prediction"] * len(roots_reference),
                        "sign_accuracy": float(np.mean(np.sign(observed) == np.sign(reference))),
                        "balanced_sign_accuracy": _balanced_sign_accuracy(observed, reference),
                        "reference_Tc_K": roots_reference,
                        "predicted_Tc_K": roots_predicted if not degenerate else ["n/a:degenerate_prediction"],
                        "predicted_Tc_by_seed_K": seed_roots if not degenerate else [["n/a:degenerate_prediction"] for _ in seeds],
                        "Tc_scatter_K": [float(np.nanstd([r[i] for r in seed_roots if len(r) > i])) if any(len(r) > i for r in seed_roots) else math.nan for i in range(len(roots_reference))] if not degenerate else ["n/a:degenerate_prediction"] * len(roots_reference),
                        "Tc_error_K": ["n/a:missed_crossing"] * len(roots_reference) if missed else ([pred - ref for pred, ref in zip(roots_predicted, roots_reference)] if not degenerate else ["n/a:degenerate_prediction"] * len(roots_reference)),
                        "crossing_slope_eV_per_atom_per_K": slopes if not degenerate else ["n/a:degenerate_prediction"] * len(roots_reference),
                        "Tc_err_from_dG_K": ["n/a:missed_crossing"] * len(roots_reference) if missed else ([error / slope if slope else math.nan for error, slope in zip(crossing_mae, slopes)] if not degenerate else ["n/a:degenerate_prediction"] * len(roots_reference)),
                        "false_crossings": max(0, len(roots_predicted) - len(roots_reference)) if not degenerate else "n/a:degenerate_prediction",
                        "missed_crossings": max(0, len(roots_reference) - len(roots_predicted)) if not degenerate else "n/a:degenerate_prediction",
                    }
    pair_values = [pair for pair in pair_metrics.values() if isinstance(pair.get("delta_G_MAE_eV_per_atom"), (int, float))]
    return {
        "predictor": spec,
        "n_test_frames": len(test),
        "seeds": seeds,
        "G_MAE_eV_per_atom": float(np.mean(point_errors)) if point_errors else math.nan,
        "ranking_accuracy": float(np.mean(ranking_values)) if ranking_values else math.nan,
        "balanced_sign_accuracy": float(np.mean([float(pair["balanced_sign_accuracy"]) for pair in pair_values])) if pair_values else math.nan,
        "coverage_2sigma": float(np.mean(coverage_values)) if coverage_values else math.nan,
        "pairs": pair_metrics,
    }


def evaluate(spec: str, split: dict[str, object], data_root: Path, seeds: list[int]) -> dict[str, object]:
    if "folds" in split:
        raw_folds = split["folds"]
        if not isinstance(raw_folds, dict):
            raise ValueError("split folds must be an object")
        folds = {
            str(name): evaluate(spec, fold, data_root, seeds)
            for name, fold in raw_folds.items()
            if isinstance(fold, dict) and isinstance(fold.get("test"), list)
        }
        total_frames = sum(int(result["n_test_frames"]) for result in folds.values())
        weighted_mae = (
            sum(float(result["G_MAE_eV_per_atom"]) * int(result["n_test_frames"]) for result in folds.values())
            / total_frames
            if total_frames
            else math.nan
        )
        weighted_coverage = (
            sum(float(result["coverage_2sigma"]) * int(result["n_test_frames"]) for result in folds.values())
            / total_frames
            if total_frames
            else math.nan
        )
        return {
            "predictor": spec,
            "n_test_frames": total_frames,
            "seeds": seeds,
            "G_MAE_eV_per_atom": weighted_mae,
            "coverage_2sigma": weighted_coverage,
            "folds": folds,
        }
    test = split.get("test")
    if not isinstance(test, list):
        raise ValueError("split must provide a test list")
    train = split.get("train", [])
    if not isinstance(train, list):
        raise ValueError("split train entries must be a list")
    # A phase-LOPO fold evaluates only its held phase.  To calculate a physical
    # ΔG and Tc it must still be paired with predictions for the known phases
    # at the same state points.  Those partner frames remain outside scalar G
    # test metrics, but are explicitly included in pairwise curve metrics.
    pair_rows = test + train
    cache = _table_cache(data_root, pair_rows)
    # The frozen v1 evaluator is temperature-only. Keep its exact behavior for
    # existing single-pressure datasets, but route multi-pressure references
    # through a state-keyed path so equal temperatures at different pressures
    # are never merged or paired.
    pressures_by_system: dict[str, set[float]] = defaultdict(set)
    for (system, _), table in cache.items():
        pressures_by_system[system].update(point.P_GPa for point in table)
    if any(len(pressures) > 1 for pressures in pressures_by_system.values()):
        return _evaluate_pressure_grid(spec, test, train, cache, seeds)
    precomputed = _load_precomputed(spec)
    by_system_phase: dict[tuple[str, str], list[int]] = defaultdict(list)
    test_by_system_phase: dict[tuple[str, str], list[int]] = defaultdict(list)
    for row in pair_rows:
        by_system_phase[(str(row["system"]), str(row["phase"]))].append(int(row["T_index"]))
    for row in test:
        test_by_system_phase[(str(row["system"]), str(row["phase"]))].append(int(row["T_index"]))
    point_errors: list[float] = []
    pair_metrics: dict[str, dict[str, object]] = {}
    coverage_values: list[bool] = []
    ranking_values: list[bool] = []
    systems = sorted({system for system, _ in test_by_system_phase})
    for system in systems:
        phases = sorted(phase for item_system, phase in by_system_phase if item_system == system)
        predictions: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
        for phase in phases:
            indexes = sorted(by_system_phase[(system, phase)])
            reference = np.array([cache[(system, phase)][index].G_eV_per_atom for index in indexes])
            temperatures = np.array([cache[(system, phase)][index].T_K for index in indexes])
            if precomputed is not None:
                curve_key = f"{system}:{phase}"
                if curve_key not in precomputed:
                    raise ValueError(f"precomputed predictor has no curve for {curve_key!r}")
                curve = precomputed[curve_key]
                missing = [float(t) for t in temperatures if float(t) not in curve]
                if missing:
                    raise ValueError(f"precomputed predictor {curve_key!r} misses temperatures {missing[:3]!r}")
                mean = np.array([curve[float(t)] for t in temperatures])
                seed_predictions = np.stack([mean.copy() for _ in seeds])
                std = np.zeros_like(mean)
            elif spec == "reference":
                seed_predictions = np.stack([reference.copy() for _ in seeds])
                mean, std = reference.copy(), np.zeros_like(reference)
            else:
                seed_predictions = np.stack(
                    [_predictor(spec, seed, f"{system}:{phase}", reference) for seed in seeds]
                )
                mean, std = seed_predictions.mean(axis=0), seed_predictions.std(axis=0)
            predictions[phase] = (temperatures, mean, seed_predictions)
            test_indexes = sorted(test_by_system_phase.get((system, phase), []))
            if test_indexes:
                test_temperatures = np.array([cache[(system, phase)][index].T_K for index in test_indexes])
                test_reference = np.array([cache[(system, phase)][index].G_eV_per_atom for index in test_indexes])
                prediction_map = dict(zip(temperatures.tolist(), mean.tolist()))
                std_map = dict(zip(temperatures.tolist(), std.tolist()))
                test_mean = np.array([prediction_map[value] for value in test_temperatures])
                test_std = np.array([std_map[value] for value in test_temperatures])
                point_errors.extend(np.abs(test_mean - test_reference).tolist())
                coverage_values.extend(np.abs(test_mean - test_reference) <= 2.0 * test_std + 1.0e-14)
        evaluation_temperatures = sorted({float(cache[(system, phase)][index].T_K) for (item_system, phase), indexes in test_by_system_phase.items() if item_system == system for index in indexes})
        for temperature in evaluation_temperatures:
            correct = True
            for left_index, left in enumerate(phases):
                for right in phases[left_index + 1 :]:
                    left_reference = next(point.G_eV_per_atom for point in cache[(system, left)] if point.T_K == temperature)
                    right_reference = next(point.G_eV_per_atom for point in cache[(system, right)] if point.T_K == temperature)
                    left_t, left_mean, _ = predictions[left]
                    right_t, right_mean, _ = predictions[right]
                    observed = float(dict(zip(left_t.tolist(), left_mean.tolist()))[temperature] - dict(zip(right_t.tolist(), right_mean.tolist()))[temperature])
                    if np.sign(observed) != np.sign(left_reference - right_reference):
                        correct = False
                        break
                if not correct:
                    break
            ranking_values.append(correct)
        for left_index, left in enumerate(phases):
            for right in phases[left_index + 1 :]:
                if (system, left) not in test_by_system_phase and (system, right) not in test_by_system_phase:
                    continue
                left_t, left_g, left_seed_g = predictions[left]
                right_t, right_g, right_seed_g = predictions[right]
                test_shared_indices = set(test_by_system_phase.get((system, left), [])) & set(test_by_system_phase.get((system, right), []))
                test_shared = {float(cache[(system, left)][index].T_K) for index in test_shared_indices}
                shared = sorted(test_shared or (set(left_t.tolist()) & set(right_t.tolist())))
                if not shared:
                    continue
                left_map = dict(zip(left_t.tolist(), left_g.tolist()))
                right_map = dict(zip(right_t.tolist(), right_g.tolist()))
                temperatures = np.array(shared)
                observed = np.array([left_map[t] - right_map[t] for t in shared])
                ref_left = {(cache[(system, left)][i].T_K): cache[(system, left)][i].G_eV_per_atom for i in by_system_phase[(system, left)]}
                ref_right = {(cache[(system, right)][i].T_K): cache[(system, right)][i].G_eV_per_atom for i in by_system_phase[(system, right)]}
                reference = np.array([ref_left[t] - ref_right[t] for t in shared])
                roots_reference = _root(reference, temperatures)
                # Some source tables are quantized and contain a flat zero
                # plateau.  A smooth gauge perturbation can split that same
                # plateau into adjacent numerical roots; merge only these
                # near-neighbours, leaving all other crossings intact.
                roots_predicted = _root(
                    observed,
                    temperatures,
                    merge_within_K=3.0 if np.any(reference == 0) else 0.0,
                )
                seed_roots = []
                for seed_index in range(len(seeds)):
                    seed_observed = np.array(
                        [
                            dict(zip(left_t.tolist(), left_seed_g[seed_index].tolist()))[t]
                            - dict(zip(right_t.tolist(), right_seed_g[seed_index].tolist()))[t]
                            for t in shared
                        ]
                    )
                    seed_roots.append(
                        _root(
                            seed_observed,
                            temperatures,
                            merge_within_K=3.0 if np.any(reference == 0) else 0.0,
                        )
                    )
                delta_mae = float(np.mean(np.abs(observed - reference)))
                slope_by_root: list[float] = []
                crossing_mae: list[float] = []
                for root in roots_reference:
                    # The design calls for the local ΔG error at each
                    # reference crossing, not a whole-curve average.  Linear
                    # interpolation keeps this well-defined on nonuniform
                    # grids and makes the reported slope conversion auditable.
                    crossing_error = float(
                        np.interp(root, temperatures, observed - reference)
                    )
                    crossing_mae.append(abs(crossing_error))
                    nearest = int(np.argmin(np.abs(temperatures - root)))
                    if reference[nearest] == 0:
                        lower, upper = nearest, nearest
                        while lower > 0 and reference[lower] == 0:
                            lower -= 1
                        while upper < len(reference) - 1 and reference[upper] == 0:
                            upper += 1
                        slope = (reference[upper] - reference[lower]) / (temperatures[upper] - temperatures[lower])
                    elif nearest == 0:
                        slope = (reference[1] - reference[0]) / (temperatures[1] - temperatures[0])
                    elif nearest == len(temperatures) - 1:
                        slope = (reference[-1] - reference[-2]) / (temperatures[-1] - temperatures[-2])
                    else:
                        slope = (reference[nearest + 1] - reference[nearest - 1]) / (
                            temperatures[nearest + 1] - temperatures[nearest - 1]
                        )
                    slope_by_root.append(float(abs(slope)))
                degenerate = bool(np.all(np.abs(observed) <= 1.0e-12))
                missed = (not degenerate) and len(roots_predicted) < len(roots_reference)
                pair_metrics[f"{system}:{left}_minus_{right}"] = {
                    "pair_support": "test_and_train_partner" if train else "test_only",
                    "status": "degenerate_prediction" if degenerate else "ok",
                    "delta_G_MAE_eV_per_atom": delta_mae,
                    "delta_G_RMSE_eV_per_atom": float(np.sqrt(np.mean((observed - reference) ** 2))),
                    "delta_G_MAE_at_crossing_eV_per_atom": crossing_mae if not degenerate else ["n/a:degenerate_prediction"] * len(roots_reference),
                    "sign_accuracy": float(np.mean(np.sign(observed) == np.sign(reference))),
                    "balanced_sign_accuracy": _balanced_sign_accuracy(observed, reference),
                    "reference_Tc_K": roots_reference,
                    "predicted_Tc_K": roots_predicted if not degenerate else ["n/a:degenerate_prediction"],
                    "predicted_Tc_by_seed_K": seed_roots if not degenerate else [["n/a:degenerate_prediction"] for _ in seeds],
                    "Tc_scatter_K": [
                        float(np.nanstd([roots[kk] for roots in seed_roots if len(roots) > kk]))
                        if any(len(roots) > kk for roots in seed_roots)
                        else math.nan
                        for kk in range(len(roots_reference))
                    ] if not degenerate else ["n/a:degenerate_prediction"] * len(roots_reference),
                    "Tc_error_K": ["n/a:missed_crossing"] * len(roots_reference) if missed else ([pred - ref for pred, ref in zip(roots_predicted, roots_reference)] if not degenerate else ["n/a:degenerate_prediction"] * len(roots_reference)),
                    "crossing_slope_eV_per_atom_per_K": slope_by_root if not degenerate else ["n/a:degenerate_prediction"] * len(roots_reference),
                    "Tc_err_from_dG_K": ["n/a:missed_crossing"] * len(roots_reference) if missed else ([
                        error / slope if slope else math.nan
                        for error, slope in zip(crossing_mae, slope_by_root)
                    ] if not degenerate else ["n/a:degenerate_prediction"] * len(roots_reference)),
                    "false_crossings": max(0, len(roots_predicted) - len(roots_reference)) if not degenerate else "n/a:degenerate_prediction",
                    "missed_crossings": max(0, len(roots_reference) - len(roots_predicted)) if not degenerate else "n/a:degenerate_prediction",
                }
    pair_values = [pair for pair in pair_metrics.values() if isinstance(pair, dict) and isinstance(pair.get("delta_G_MAE_eV_per_atom"), (int, float))]
    return {
        "predictor": spec,
        "n_test_frames": len(test),
        "seeds": seeds,
        "G_MAE_eV_per_atom": float(np.mean(point_errors)) if point_errors else math.nan,
        "ranking_accuracy": float(np.mean(ranking_values)) if ranking_values else math.nan,
        "balanced_sign_accuracy": float(np.mean([float(pair["balanced_sign_accuracy"]) for pair in pair_values])) if pair_values else math.nan,
        "coverage_2sigma": float(np.mean(coverage_values)) if coverage_values else math.nan,
        "pairs": pair_metrics,
    }


def _number(value: object) -> str:
    if isinstance(value, (float, int)):
        return f"{value:.6g}"
    return "—"


def _summary_markdown(metrics: dict[str, object], title: str) -> str:
    """Render scalar and pairwise metrics as a portable Markdown table."""

    lines = [f"# {title}", "", f"Frames: {metrics['n_test_frames']}", "", f"G MAE: {_number(metrics['G_MAE_eV_per_atom'])} eV/atom", ""]
    folds = metrics.get("folds")
    metric_views: list[tuple[str, dict[str, object]]]
    if isinstance(folds, dict):
        metric_views = [(str(name), value) for name, value in folds.items() if isinstance(value, dict)]
    else:
        metric_views = [("all test frames", metrics)]
    for fold_name, view in metric_views:
        pairs = view.get("pairs", {})
        if not isinstance(pairs, dict) or not pairs:
            lines.extend([f"## {fold_name}", "", "No comparable phase pair is defined for this fold.", ""])
            continue
        lines.extend(
            [
                f"## {fold_name}",
                "",
                "| Pair | support | ΔG MAE (eV/atom) | ΔG RMSE (eV/atom) | ΔG MAE@Tc (eV/atom) | sign accuracy | reference Tc (K) | predicted Tc (K) | Tc error (K) | Tc from ΔG (K) | Tc scatter (K) | false / missed |",
                "|---|---|---:|---:|---|---:|---|---|---|---|---|---:|",
            ]
        )
        for pair_name, pair in sorted(pairs.items()):
            if not isinstance(pair, dict):
                continue
            crossing = lambda key: ", ".join(_number(value) for value in pair.get(key, []) if isinstance(value, (float, int))) or "—"
            lines.append(
                "| "
                + " | ".join(
                    (
                        str(pair_name),
                        str(pair.get("pair_support", "test_only")),
                        _number(pair.get("delta_G_MAE_eV_per_atom")),
                        _number(pair.get("delta_G_RMSE_eV_per_atom")),
                        crossing("delta_G_MAE_at_crossing_eV_per_atom"),
                        _number(pair.get("sign_accuracy")),
                        crossing("reference_Tc_K"),
                        crossing("predicted_Tc_K"),
                        crossing("Tc_error_K"),
                        crossing("Tc_err_from_dG_K"),
                        crossing("Tc_scatter_K"),
                        f"{pair.get('false_crossings', 0)} / {pair.get('missed_crossings', 0)}",
                    )
                )
                + " |"
            )
        lines.append("")
    return "\n".join(lines)


def run(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    split_path, data_root = Path(args.split).resolve(), Path(args.data_root).resolve()
    split = json.loads(split_path.read_text(encoding="utf-8"))
    seeds = [int(item) for item in args.seeds.split(",") if item]
    metrics = evaluate(args.predictor, split, data_root, seeds)
    out = Path(args.output_root).resolve() / _output_slug(args.predictor) / split_path.stem
    out.mkdir(parents=True, exist_ok=True)
    (out / "metrics.json").write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (out / "summary.md").write_text(
        _summary_markdown(metrics, f"{args.predictor} on {split_path.stem}"), encoding="utf-8"
    )
    print(out / "metrics.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
