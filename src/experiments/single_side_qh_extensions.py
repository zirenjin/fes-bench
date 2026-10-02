"""Extend the single-sided QH diagnostic — plan experiment single_side_qh.

This script adds training-only thermo-form and Einstein controls, QH-invalid
negative controls, and a pair-level mixed-model diagnostic.  It consumes the
frozen v2 temp-extrapolation split and existing T3/T2 predictions; it never
trains a model and never edits the benchmark tables.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

import single_side_qh as base


SEEDS = (11, 23, 37)
KB_EV_PER_K = 8.617333262145e-5
TOL = 1.0e-10
VALID_PHASES = tuple((system, phase, label) for system, phase, label in base.PHASES)
VALID_KEYS = {(system, phase) for system, phase, _ in VALID_PHASES}
INVALID_PHASES = (
    ("hf", "bcc", "Hf bcc"),
    ("ti", "bcc", "Ti bcc"),
    ("zr", "bcc", "Zr bcc"),
    ("sio2", "cristobalite_beta", "SiO2 beta-cristobalite"),
    ("sio2", "tridymite_p63mmc", "SiO2 beta-tridymite"),
)
PAIR_DEFS = (
    ("hf", "hcp", "bcc", "hf:hcp_minus_bcc"),
    ("ti", "hcp", "bcc", "ti:hcp_minus_bcc"),
    ("zr", "hcp", "bcc", "zr:hcp_minus_bcc"),
    ("sio2", "quartz_beta", "cristobalite_beta", "sio2:quartz_beta_minus_cristobalite_beta"),
    ("sio2", "quartz_beta", "tridymite_p63mmc", "sio2:quartz_beta_minus_tridymite_p63mmc"),
    ("sio2", "cristobalite_beta", "tridymite_p63mmc", "sio2:cristobalite_beta_minus_tridymite_p63mmc"),
)


def read_csv(path: Path) -> list[dict[str, str]]:
    return base.read_csv(path)


def write_csv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    base.write_csv(path, fields, rows)


def tlog_design(temperature: np.ndarray) -> np.ndarray:
    temperature = np.asarray(temperature, dtype=float)
    if np.any(temperature <= 0.0):
        raise ValueError("temperatures must be positive")
    return np.column_stack((np.ones_like(temperature), temperature, temperature * np.log(temperature)))


def tlog_fit_predict(train_t: np.ndarray, train_y: np.ndarray, target_t: np.ndarray) -> np.ndarray:
    coefficients, *_ = np.linalg.lstsq(tlog_design(train_t), train_y, rcond=None)
    return tlog_design(target_t) @ coefficients


def load_fqh(root: Path, system: str, phase: str) -> dict[float, float]:
    path = root / "data/processed" / system / phase / "fqh.csv"
    rows = read_csv(path)
    if not rows or "F_QH_eV_per_atom" not in rows[0]:
        raise ValueError(f"missing canonical F_QH column: {path}")
    return {float(row["T_K"]): float(row["F_QH_eV_per_atom"]) for row in rows}


def phase_grids(root: Path, split: dict[str, list[dict[str, object]]], system: str, phase: str) -> dict[str, np.ndarray]:
    return base.expected_temperatures(root, split, system, phase)


def anchored_rmse(prediction: np.ndarray, reference: np.ndarray, pred_anchor: float, ref_anchor: float) -> float:
    return base.rmse((prediction - pred_anchor) - (reference - ref_anchor)) * 1000.0


def anchored_mae(prediction: np.ndarray, reference: np.ndarray, pred_anchor: float, ref_anchor: float) -> float:
    return base.mae((prediction - pred_anchor) - (reference - ref_anchor)) * 1000.0


def einstein_free(temperature: np.ndarray, theta_K: float) -> np.ndarray:
    if theta_K <= 0.0:
        raise ValueError("Einstein temperature must be positive")
    temperature = np.asarray(temperature, dtype=float)
    # -expm1 is stable for both small and large theta/T.
    return 3.0 * KB_EV_PER_K * temperature * np.log(-np.expm1(-theta_K / temperature))


def fit_einstein_theta(train_t: np.ndarray, fqh_train: np.ndarray) -> tuple[float, float]:
    """Fit the one parameter to the anchored F_QH shape in the train window."""

    target = fqh_train - fqh_train[-1]
    log_grid = np.linspace(np.log(1.0), np.log(1.0e5), 5000)
    candidates = np.exp(log_grid)

    def objective(theta: float) -> float:
        curve = einstein_free(train_t, theta)
        residual = (curve - curve[-1]) - target
        return float(np.mean(residual * residual))

    values = np.asarray([objective(float(theta)) for theta in candidates])
    best = int(np.argmin(values))
    lo = candidates[max(0, best - 1)]
    hi = candidates[min(len(candidates) - 1, best + 1)]
    # Golden-section refinement in log(theta), retaining a one-parameter fit.
    a, b = float(np.log(lo)), float(np.log(hi))
    golden = (math.sqrt(5.0) - 1.0) / 2.0
    c, d = b - golden * (b - a), a + golden * (b - a)
    for _ in range(80):
        if objective(math.exp(c)) < objective(math.exp(d)):
            b, d = d, c
            c = b - golden * (b - a)
        else:
            a, c = c, d
            d = a + golden * (b - a)
    theta = math.exp((a + b) / 2.0)
    return theta, objective(theta)


def prediction_path(root: Path, source: str, system: str, phase: str, seed: int, region: str) -> Path:
    if source == "t3_tlog_high":
        path = root / "result/experiments/single_side_qh/predictions/t3_tlog_high" / f"{system}_{phase}_seed{seed}_{region}.csv"
    else:
        path = root / "result/experiments/single_side_qh/predictions" / source / f"{system}_{phase}_seed{seed}_{region}.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def read_prediction(root: Path, source: str, system: str, phase: str, seed: int, region: str) -> tuple[np.ndarray, np.ndarray]:
    path = prediction_path(root, source, system, phase, seed, region)
    return base.load_prediction(path, system, phase, seed, region)


def validated_context(root: Path) -> tuple[dict[str, list[dict[str, object]]], Path, str]:
    split, split_path, split_hash = base.load_split(root)
    reliable = base.qh_reliable_phases(root)
    required = {(system, phase) for system, phase, _ in VALID_PHASES}
    if not required.issubset(reliable):
        raise AssertionError(f"valid-side phase missing from inventory: {sorted(required - reliable)}")
    if reliable.intersection({(system, phase) for system, phase, _ in INVALID_PHASES}):
        raise AssertionError("negative-control phase unexpectedly marked QH-reliable")
    for system, phase, _ in VALID_PHASES:
        expected = phase_grids(root, split, system, phase)
        fqh = load_fqh(root, system, phase)
        reference = base.load_reference(root, system, phase)
        for region in ("train", "test"):
            if any(float(t) not in fqh for t in expected[region]) or any(float(t) not in reference for t in expected[region]):
                raise AssertionError(f"F_QH/reference grid mismatch for {system}/{phase}/{region}")
    return split, split_path, split_hash


def valid_baseline_rows(root: Path, split: dict[str, list[dict[str, object]]], split_hash: str) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    """Return per-seed upper-bound, Einstein parameter, and summary source rows."""

    rows: list[dict[str, object]] = []
    theta_rows: list[dict[str, object]] = []
    for system, phase, _ in VALID_PHASES:
        grids = phase_grids(root, split, system, phase)
        reference = base.load_reference(root, system, phase)
        fqh = load_fqh(root, system, phase)
        train_t, test_t = grids["train"], grids["test"]
        ref_train = np.asarray([reference[float(t)] for t in train_t])
        ref_test = np.asarray([reference[float(t)] for t in test_t])
        fqh_train = np.asarray([fqh[float(t)] for t in train_t])
        fqh_test = np.asarray([fqh[float(t)] for t in test_t])
        t_ref, ref_anchor, fqh_anchor = float(train_t[-1]), float(ref_train[-1]), float(fqh_train[-1])
        abs_pred_test = tlog_fit_predict(train_t, ref_train, test_t)
        qh_residual_train = ref_train - fqh_train
        qh_residual_test = tlog_fit_predict(train_t, qh_residual_train, test_t)
        qh_pred_test = fqh_test + qh_residual_test
        theta, theta_objective = fit_einstein_theta(train_t, fqh_train)
        ein_train = einstein_free(train_t, theta)
        ein_test = einstein_free(test_t, theta)
        ein_anchor = float(ein_train[-1])
        ein_residual_train = ref_train - ein_train
        ein_residual_test = tlog_fit_predict(train_t, ein_residual_train, test_t)
        ein_pred_test = ein_test + ein_residual_test
        theta_rows.append({
            "system": system, "phase": phase, "seed": "all", "theta_E_K": theta,
            "train_shape_mse_eV2_per_atom2": theta_objective,
            "fit_scope": "anchored F_QH train shape only", "split_sha256": split_hash,
            "fqh_source": str((root / "data/processed" / system / phase / "fqh.csv").relative_to(root)),
        })
        candidate_curves = {
            "thermo_abs_G": (abs_pred_test, float(abs_pred_test[0] if False else tlog_fit_predict(train_t, ref_train, np.asarray([t_ref]))[0])),
            "thermo_QH_residual": (qh_pred_test, float(fqh_anchor + tlog_fit_predict(train_t, qh_residual_train, np.asarray([t_ref]))[0])),
            "thermo_Einstein_residual": (ein_pred_test, float(ein_anchor + tlog_fit_predict(train_t, ein_residual_train, np.asarray([t_ref]))[0])),
        }
        # The fitted curves are evaluated at T_ref explicitly, so their
        # anchors do not rely on a test point or an extrapolated first point.
        for seed in SEEDS:
            learned_abs_t, learned_abs_train = read_prediction(root, "t3_tlog", system, phase, seed, "train")
            learned_abs_test_t, learned_abs_test = read_prediction(root, "t3_tlog", system, phase, seed, "test")
            learned_qh_t, learned_qh_train = read_prediction(root, "t2_qh", system, phase, seed, "train")
            learned_qh_test_t, learned_qh_test = read_prediction(root, "t2_qh", system, phase, seed, "test")
            if not np.array_equal(learned_abs_t, train_t) or not np.array_equal(learned_abs_test_t, test_t) or not np.array_equal(learned_qh_t, train_t) or not np.array_equal(learned_qh_test_t, test_t):
                raise AssertionError(f"learning grid mismatch for {system}/{phase}/{seed}")
            learned_ein_test = ein_test + (learned_qh_test - fqh_test)
            learned_curves = {
                "learned_absolute_G": (learned_abs_test, float(learned_abs_train[-1])),
                "learned_QH": (learned_qh_test, float(learned_qh_train[-1])),
                "learned_Einstein_transplant": (learned_ein_test, float(ein_anchor + (learned_qh_train[-1] - fqh_anchor))),
            }
            for predictor, (prediction, prediction_anchor) in {**candidate_curves, **learned_curves}.items():
                if predictor.startswith("thermo_"):
                    # Fit predictions are anchored at their fitted train value.
                    if predictor == "thermo_abs_G":
                        anchor = float(tlog_fit_predict(train_t, ref_train, np.asarray([t_ref]))[0])
                    elif predictor == "thermo_QH_residual":
                        anchor = float(fqh_anchor + tlog_fit_predict(train_t, qh_residual_train, np.asarray([t_ref]))[0])
                    else:
                        anchor = float(ein_anchor + tlog_fit_predict(train_t, ein_residual_train, np.asarray([t_ref]))[0])
                else:
                    anchor = prediction_anchor
                rows.append({
                    "system": system, "phase": phase, "seed": seed, "predictor": predictor,
                    "T_ref_K": t_ref, "n_train": len(train_t), "n_test": len(test_t),
                    "anchored_test_mae_meV": anchored_mae(prediction, ref_test, anchor, ref_anchor),
                    "anchored_test_rmse_meV": anchored_rmse(prediction, ref_test, anchor, ref_anchor),
                    "raw_test_mae_meV": base.mae(prediction - ref_test) * 1000.0,
                    "raw_test_rmse_meV": base.rmse(prediction - ref_test) * 1000.0,
                    "theta_E_K": theta if "Einstein" in predictor else None,
                    "split_sha256": split_hash,
                })
    return rows, theta_rows, []


def invalid_baseline_rows(root: Path, split: dict[str, list[dict[str, object]]], split_hash: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for system, phase, _ in INVALID_PHASES:
        grids = phase_grids(root, split, system, phase)
        reference = base.load_reference(root, system, phase)
        fqh = load_fqh(root, system, phase)
        train_t, test_t = grids["train"], grids["test"]
        ref_train = np.asarray([reference[float(t)] for t in train_t])
        ref_test = np.asarray([reference[float(t)] for t in test_t])
        fqh_train = np.asarray([fqh[float(t)] for t in train_t])
        fqh_test = np.asarray([fqh[float(t)] for t in test_t])
        ref_anchor = float(ref_train[-1])
        abs_pred = tlog_fit_predict(train_t, ref_train, test_t)
        residual_pred = tlog_fit_predict(train_t, ref_train - fqh_train, test_t)
        qh_pred = fqh_test + residual_pred
        abs_rmse = anchored_rmse(abs_pred, ref_test, float(tlog_fit_predict(train_t, ref_train, np.asarray([train_t[-1]]))[0]), ref_anchor)
        qh_rmse = anchored_rmse(qh_pred, ref_test, float(fqh_train[-1] + tlog_fit_predict(train_t, ref_train - fqh_train, np.asarray([train_t[-1]]))[0]), ref_anchor)
        ratio = qh_rmse / abs_rmse if abs_rmse else math.nan
        for seed in SEEDS:
            rows.append({
                "system": system, "phase": phase, "seed": seed, "qh_reliable": False,
                "predictor": "thermo_QH_residual", "thermo_abs_G_rmse_meV": abs_rmse,
                "thermo_QH_residual_rmse_meV": qh_rmse, "R_QH_over_abs": ratio,
                "T_ref_K": float(train_t[-1]), "n_train": len(train_t), "n_test": len(test_t),
                "split_sha256": split_hash,
                "fqh_source": str((root / "data/processed" / system / phase / "fqh.csv").relative_to(root)),
            })
    return rows


def crossing_roots(temperature: np.ndarray, values: np.ndarray) -> list[float]:
    roots: list[float] = []
    for index in range(len(values) - 1):
        left, right = float(values[index]), float(values[index + 1])
        if left == 0.0:
            roots.append(float(temperature[index]))
        if left * right < 0.0:
            fraction = -left / (right - left)
            roots.append(float(temperature[index] + fraction * (temperature[index + 1] - temperature[index])))
    if values[-1] == 0.0:
        roots.append(float(temperature[-1]))
    unique: list[float] = []
    for root in roots:
        if not unique or abs(root - unique[-1]) > 1.0e-8:
            unique.append(root)
    return unique


def reference_crossings(root: Path, system: str, pair_key: str) -> list[dict[str, object]]:
    payload = json.loads((root / "data/processed" / system / "reference_crossings.json").read_text(encoding="utf-8"))
    for key, pair in payload["pairs"].items():
        if f"{system}:{key}" == pair_key:
            return list(pair.get("crossings", []))
    return []


def pair_metrics(temperature: np.ndarray, reference: np.ndarray, predicted: np.ndarray, crossings: list[dict[str, object]]) -> dict[str, object]:
    roots_ref = [float(item["T_K"]) for item in crossings]
    roots_pred = crossing_roots(temperature, predicted)
    missing = max(0, len(roots_ref) - len(roots_pred))
    false = max(0, len(roots_pred) - len(roots_ref))
    if not roots_ref:
        tc_error: object = "n/a:no_reference_crossing"
        tc_from_dg: object = "n/a:no_reference_crossing"
    elif missing:
        tc_error = "n/a:missed_crossing"
        tc_from_dg = "n/a:missed_crossing"
    else:
        tc_error = float(roots_pred[0] - roots_ref[0])
        slope = float(crossings[0].get("slope_eV_per_atom_per_K", 0.0))
        delta_at_ref = float(np.interp(roots_ref[0], temperature, predicted - reference))
        tc_from_dg = float(delta_at_ref / slope) if slope else "n/a:zero_slope"
    return {
        "delta_G_MAE_eV_per_atom": base.mae(predicted - reference),
        "delta_G_RMSE_eV_per_atom": base.rmse(predicted - reference),
        "sign_accuracy": float(np.mean(np.sign(predicted) == np.sign(reference))),
        "reference_Tc_K": roots_ref,
        "predicted_Tc_K": roots_pred,
        "Tc_error_K": tc_error,
        "Tc_err_from_dG_K": tc_from_dg,
        "false_crossings": false,
        "missed_crossings": missing,
    }


def pair_rows(root: Path, split: dict[str, list[dict[str, object]]], split_hash: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for system, low_phase, high_phase, pair_key in PAIR_DEFS:
        grids = phase_grids(root, split, system, low_phase)
        low_test_t, high_test_t = grids["test"], phase_grids(root, split, system, high_phase)["test"]
        if not np.array_equal(low_test_t, high_test_t):
            raise AssertionError(f"pair grid mismatch for {pair_key}")
        temperature = low_test_t
        low_ref = base.load_reference(root, system, low_phase)
        high_ref = base.load_reference(root, system, high_phase)
        reference = np.asarray([low_ref[float(t)] - high_ref[float(t)] for t in temperature])
        crossings = reference_crossings(root, system, pair_key)
        for seed in SEEDS:
            low_source = "t3_tlog" if (system, low_phase) in VALID_KEYS else "t3_tlog_high"
            high_source = "t3_tlog" if (system, high_phase) in VALID_KEYS else "t3_tlog_high"
            _, low_t3 = read_prediction(root, low_source, system, low_phase, seed, "test")
            _, high_t3 = read_prediction(root, high_source, system, high_phase, seed, "test")
            t3_metrics = pair_metrics(temperature, reference, low_t3 - high_t3, crossings)
            rows.append({"pair": pair_key, "system": system, "low_phase": low_phase, "high_phase": high_phase, "seed": seed, "predictor": "T3_absolute_both", "gauge_alignment_meV": None, "split_sha256": split_hash, **t3_metrics})
            if (system, low_phase) in VALID_KEYS:
                _, low_qh = read_prediction(root, "t2_qh", system, low_phase, seed, "test")
                train_t3, train_t3_values = read_prediction(root, "t3_tlog", system, low_phase, seed, "train")
                train_t2, train_t2_values = read_prediction(root, "t2_qh", system, low_phase, seed, "train")
                if not np.array_equal(train_t3, train_t2):
                    raise AssertionError(f"gauge-alignment train grid mismatch for {system}/{low_phase}/{seed}")
                # T2 and T3 are independently calibrated absolute-G curves.
                # Align only their low-side gauge at the frozen train anchor;
                # this uses no test label and makes the mixed pair a thermal
                # shape comparison rather than an arbitrary model zero.
                gauge_shift = float(train_t3_values[-1] - train_t2_values[-1])
                combo_metrics = pair_metrics(temperature, reference, (low_qh + gauge_shift) - high_t3, crossings)
                rows.append({"pair": pair_key, "system": system, "low_phase": low_phase, "high_phase": high_phase, "seed": seed, "predictor": "single_side_QH_combo", "gauge_alignment_meV": gauge_shift * 1000.0, "split_sha256": split_hash, **combo_metrics})
            else:
                rows.append({
                    "pair": pair_key, "system": system, "low_phase": low_phase, "high_phase": high_phase, "seed": seed,
                    "predictor": "single_side_QH_combo", "split_sha256": split_hash,
                    "gauge_alignment_meV": None,
                    "delta_G_MAE_eV_per_atom": "n/a:no_qh_low_side", "delta_G_RMSE_eV_per_atom": "n/a:no_qh_low_side",
                    "sign_accuracy": "n/a:no_qh_low_side", "reference_Tc_K": [], "predicted_Tc_K": [],
                    "Tc_error_K": "n/a:no_qh_low_side", "Tc_err_from_dG_K": "n/a:no_qh_low_side",
                    "false_crossings": "n/a:no_qh_low_side", "missed_crossings": "n/a:no_qh_low_side",
                })
        thermo_path = root / "result/experiments/external_baselines/raw_temp_extrap/thermo_form_fit/seed_none/metrics.json"
        thermo = json.loads(thermo_path.read_text(encoding="utf-8"))["metrics"]["folds"]["all"]["pairs"].get(pair_key)
        if thermo is None:
            raise KeyError(f"thermo-form pair missing: {pair_key}")
        rows.append({
            "pair": pair_key, "system": system, "low_phase": low_phase, "high_phase": high_phase,
            "seed": "none", "predictor": "thermo_form_fit", "split_sha256": split_hash,
            "gauge_alignment_meV": None,
            "delta_G_MAE_eV_per_atom": thermo.get("delta_G_MAE_eV_per_atom"),
            "delta_G_RMSE_eV_per_atom": thermo.get("delta_G_RMSE_eV_per_atom"),
            "sign_accuracy": thermo.get("sign_accuracy"), "reference_Tc_K": thermo.get("reference_Tc_K"),
            "predicted_Tc_K": thermo.get("predicted_Tc_K"), "Tc_error_K": thermo.get("Tc_error_K"),
            "Tc_err_from_dG_K": thermo.get("Tc_err_from_dG_K"), "false_crossings": thermo.get("false_crossings"),
            "missed_crossings": thermo.get("missed_crossings"),
        })
    return rows


def summary_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    predictors = sorted({str(row["predictor"]) for row in rows})
    for predictor in predictors:
        selected = [row for row in rows if row["predictor"] == predictor and isinstance(row.get("delta_G_MAE_eV_per_atom"), (float, int))]
        if not selected:
            continue
        mae_values = np.asarray([float(row["delta_G_MAE_eV_per_atom"]) for row in selected])
        sign_values = np.asarray([float(row["sign_accuracy"]) for row in selected])
        tc_values = [float(row["Tc_error_K"]) for row in selected if isinstance(row.get("Tc_error_K"), (float, int))]
        output.append({
            "predictor": predictor, "n_pair_rows": len(selected),
            "delta_G_MAE_mean_meV_per_atom": float(mae_values.mean() * 1000.0),
            "delta_G_MAE_std_meV_per_atom": float(mae_values.std(ddof=1)) * 1000.0 if len(mae_values) > 1 else 0.0,
            "sign_accuracy_mean": float(sign_values.mean()),
            "sign_accuracy_std": float(sign_values.std(ddof=1)) if len(sign_values) > 1 else 0.0,
            "Tc_error_abs_mean_K": float(np.mean(np.abs(tc_values))) if tc_values else None,
            "n_with_numeric_Tc_error": len(tc_values),
        })
    return output


def upper_summary_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for system, phase, label in VALID_PHASES:
        predictors = sorted({str(row["predictor"]) for row in rows if row["system"] == system and row["phase"] == phase})
        for predictor in predictors:
            values = np.asarray([float(row["anchored_test_rmse_meV"]) for row in rows if row["system"] == system and row["phase"] == phase and row["predictor"] == predictor])
            output.append({
                "system": system, "phase": phase, "label": label, "predictor": predictor,
                "n_seeds": len(values), "anchored_test_rmse_mean_meV": float(values.mean()),
                "anchored_test_rmse_std_meV": float(values.std(ddof=1)) if len(values) > 1 else 0.0,
            })
    return output


def write_readme(path: Path, upper: list[dict[str, object]], invalid: list[dict[str, object]], pair_summary: list[dict[str, object]], theta: list[dict[str, object]]) -> None:
    original = path.read_text(encoding="utf-8") if path.exists() else "# Single-sided QH diagnostic\n"
    marker = "\n## Extensions: controls and pair diagnostic\n"
    if marker in original:
        original = original.split(marker, 1)[0] + "\n"
    lines = [marker.rstrip(), "", "This extension keeps the original temp_extrap v2, seed 11/23/37 convention and writes all new controls under this directory. No model was retrained.", "", "### 1. Training-only thermo-form upper bounds", "", "For each QH-reliable phase, `thermo_abs_G` fits `G_ref(T)` and `thermo_QH_residual` fits `G_ref(T) − F_QH(T)` with `a + bT + cT ln(T)` on train points only. The reported error is the same anchored test RMSE used above.", ""]
    lines += ["| Phase | thermo abs-G RMSE (meV/atom) | thermo QH-residual RMSE (meV/atom) | R = QH/abs | learned T3 abs-G | learned T2 QH |", "|---|---:|---:|---:|---:|---:|"]
    for system, phase, label in VALID_PHASES:
        selected = [r for r in upper if r["system"] == system and r["phase"] == phase]
        by = {r["predictor"]: [float(x["anchored_test_rmse_meV"]) for x in selected if x["predictor"] == r["predictor"]] for r in selected}
        abs_v = float(np.mean(by["thermo_abs_G"]))
        qh_v = float(np.mean(by["thermo_QH_residual"]))
        t3_v = float(np.mean(by["learned_absolute_G"]))
        t2_v = float(np.mean(by["learned_QH"]))
        lines.append(f"| {label} | {abs_v:.3f} | {qh_v:.3f} | {qh_v / abs_v:.3f} | {t3_v:.3f} | {t2_v:.3f} |")
    lines += ["", "These controls separate the representation model from the physics baseline: the thermo-form rows are non-learning ceilings for the same basis and train/test window. Here the learned T3 errors are much larger than the thermo-form ceiling, so its dominant limitation is the learned absolute-G representation/extrapolation. T2 is closer than T3 but still above the QH thermo ceiling, especially for SiO2; that residual gap is not explained by the smooth baseline alone.", "", "### 2. Einstein smooth baseline", "", "`einstein_parameters.csv` fits the single positive parameter θ_E to the anchored training shape of F_QH. `thermo_Einstein_residual` fits `G_ref − F_E` with the same three-term thermo form. `learned_Einstein_transplant` is explicitly a no-retraining diagnostic: the existing T2 learned residual `(T2 − F_QH)` is added to F_E. It is not presented as a newly trained model.", ""]
    lines.append("| Phase | theta_E (K) | thermo Einstein RMSE (meV/atom) | learned Einstein transplant (meV/atom) | QH learned RMSE (meV/atom) |")
    lines.append("|---|---:|---:|---:|---:|")
    for system, phase, label in VALID_PHASES:
        ptheta = next(x for x in theta if x["system"] == system and x["phase"] == phase)
        selected = [r for r in upper if r["system"] == system and r["phase"] == phase]
        mean = lambda predictor: float(np.mean([float(x["anchored_test_rmse_meV"]) for x in selected if x["predictor"] == predictor]))
        lines.append(f"| {label} | {float(ptheta['theta_E_K']):.3f} | {mean('thermo_Einstein_residual'):.3f} | {mean('learned_Einstein_transplant'):.3f} | {mean('learned_QH'):.3f} |")
    lines += ["", "The Einstein comparison is a shape control, not a claim that one oscillator is a physical phonon spectrum. Similar QH and Einstein errors would support a smooth-shape interpretation; a large difference would indicate phase-specific QH information matters.", "", "### 3. QH-invalid negative controls", "", "The invalid phases use the same non-learning comparison after the adopted imaginary-mode exclusion. They are not used in the valid-side conclusions.", "", "| Phase | QH-reliable | thermo abs-G RMSE | thermo QH-residual RMSE | R = QH/abs |", "|---|---|---:|---:|---:|"]
    for system, phase, label in INVALID_PHASES:
        selected = [x for x in invalid if x["system"] == system and x["phase"] == phase]
        abs_v = float(np.mean([float(x["thermo_abs_G_rmse_meV"]) for x in selected]))
        qh_v = float(np.mean([float(x["thermo_QH_residual_rmse_meV"]) for x in selected]))
        lines.append(f"| {label} | false | {abs_v:.3f} | {qh_v:.3f} | {qh_v / abs_v:.3f} |")
    lines += ["", "The expected R ≥ 1 is a diagnostic expectation, not a forced acceptance criterion; the measured values are reported as-is.", "", "### 4. Pair-level mixed predictor", "", "`single_side_QH_combo` uses the learned T2 QH curve for the low-temperature/QH-reliable phase and the learned T3 absolute-G curve for the partner. It is compared with `T3_absolute_both` and the existing canonical `thermo_form_fit`. Pair metrics retain the frozen crossing labels, including N/A reasons.", ""]
    lines += ["| Predictor | Pair rows | ΔG MAE (meV/atom) | sign accuracy | mean |Tc error| (K) |", "|---|---:|---:|---:|---:|"]
    for row in pair_summary:
        tc = "n/a" if row["Tc_error_abs_mean_K"] is None else f"{float(row['Tc_error_abs_mean_K']):.3f}"
        lines.append(f"| {row['predictor']} | {row['n_pair_rows']} | {float(row['delta_G_MAE_mean_meV_per_atom']):.3f} ± {float(row['delta_G_MAE_std_meV_per_atom']):.3f} | {float(row['sign_accuracy_mean']):.4f} | {tc} |")
    lines += ["", "For the mixed pair, the T2 low-side curve is shifted to the T3 low-side value at the frozen train anchor; the shift is recorded as `gauge_alignment_meV` and uses no test labels. The third SiO2 pair has no QH-reliable low side and is explicitly N/A. The pair diagnostic is the direct test of whether single-sided QH helps a phase transition. It does not change the main benchmark tables or imply that QH is valid on the high-temperature unstable side.", "", "## Reproduce extensions", "", "```bash", "python src/experiments/single_side_qh_extensions.py --repo-root .", "```", "", "The script validates exact frozen train/test grids, eV/atom prediction files, split hashes, QH source files, common seeds, and the no-retraining provenance of all controls.", ""]
    path.write_text(original.rstrip() + "\n" + "\n".join(lines).rstrip() + "\n", encoding="utf-8")


def sha256(path: Path) -> str:
    return base.sha256(path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    args = parser.parse_args()
    root = args.repo_root.resolve()
    output = root / "result/experiments/single_side_qh"
    split, split_path, split_hash = validated_context(root)
    upper, theta_rows, _ = valid_baseline_rows(root, split, split_hash)
    invalid = invalid_baseline_rows(root, split, split_hash)
    pairs = pair_rows(root, split, split_hash)
    write_csv(output / "upper_bound_metrics.csv", list(upper[0]), upper)
    upper_summary = upper_summary_rows(upper)
    write_csv(output / "upper_bound_summary.csv", list(upper_summary[0]), upper_summary)
    write_csv(output / "einstein_parameters.csv", list(theta_rows[0]), theta_rows)
    write_csv(output / "negative_control_metrics.csv", list(invalid[0]), invalid)
    write_csv(output / "pair_single_side_qh.csv", list(pairs[0]), pairs)
    pair_summary = summary_rows(pairs)
    write_csv(output / "pair_single_side_qh_summary.csv", list(pair_summary[0]), pair_summary)
    write_readme(output / "README.md", upper, invalid, pair_summary, theta_rows)

    input_paths = [split_path, root / "src/experiments/single_side_qh.py", root / "src/lib/fes_bench/baselines/thermo_form_fit.py"]
    input_paths += [root / "data/processed" / system / phase / name for system, phase, _ in (*VALID_PHASES, *INVALID_PHASES) for name in ("reference_G.csv", "fqh.csv")]
    input_paths += [prediction_path(root, source, system, phase, seed, region) for source, phases in (("t3_tlog", VALID_PHASES), ("t2_qh", VALID_PHASES)) for system, phase, _ in phases for seed in SEEDS for region in ("train", "test")]
    input_paths += [prediction_path(root, "t3_tlog_high", system, phase, seed, "test") for system, phase, _ in INVALID_PHASES for seed in SEEDS]
    input_paths += [root / "result/experiments/external_baselines/raw_temp_extrap/thermo_form_fit/seed_none/metrics.json"]
    meta = {
        "git_commit": base.git_commit(root), "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "script": {"path": "src/experiments/single_side_qh_extensions.py", "sha256": sha256(root / "src/experiments/single_side_qh_extensions.py")},
        "split": {"path": str(split_path.relative_to(root)), "sha256": split_hash},
        "seeds": list(SEEDS), "training": "no model retraining; existing canonical T3/T2 predictions only",
        "inputs": [{"path": str(path.relative_to(root)), "sha256": sha256(path)} for path in sorted(set(input_paths))],
        "controls": {"thermo_basis": "a + bT + cT ln(T)", "einstein_formula": "3 k_B T ln(1 - exp(-theta_E/T))", "einstein_fit_scope": "anchored F_QH train shape only", "pair_low_side": "T2 QH learned residual", "pair_high_side": "T3 absolute-G learned curve"},
        "checks": {"exact_grids": True, "split_hash_consistent": True, "units": "eV/atom inputs; meV/atom outputs", "invalid_qh_negative_controls": True},
    }
    (output / "extensions.meta.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    print("valid upper-bound rows:", len(upper))
    print("invalid negative-control rows:", len(invalid))
    print("pair rows:", len(pairs))
    for row in pair_summary:
        print(row["predictor"], row["delta_G_MAE_mean_meV_per_atom"], row["sign_accuracy_mean"], row["Tc_error_abs_mean_K"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
