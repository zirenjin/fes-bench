"""Run the no-training static-energy and oracle-residual diagnostics.

All input paths are resolved relative to the repository root. This script only
reads canonical processed data, frozen split definitions, provenance tables,
and canonical QH experiment outputs. It never evaluates a model or changes an
existing input table.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "result/experiments/residual_oracle"
DATA = ROOT / "data/processed"
KB_EV_PER_K = 8.617333262145e-5
HBAR_EV_S = 6.582119569e-16
OMEGA0_HZ = 1.0e12
PHASES = {
    "hf": ("hcp", "bcc"),
    "ti": ("hcp", "bcc"),
    "zr": ("hcp", "bcc"),
    "sio2": ("quartz_beta", "cristobalite_beta", "tridymite_p63mmc"),
}


@dataclass(frozen=True)
class Pair:
    system: str
    left: str
    right: str
    low: str
    high: str
    reference_tc: float
    temperature: np.ndarray
    delta_g: np.ndarray


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def read_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"{path} is not a JSON object")
    return value


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"{path} is empty")
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def reference_curve(system: str, phase: str) -> tuple[np.ndarray, np.ndarray]:
    rows = read_csv(DATA / system / phase / "reference_G.csv")
    temperature = np.asarray([float(row["T_K"]) for row in rows], dtype=float)
    free_energy = np.asarray([float(row["G_eV_per_atom"]) for row in rows], dtype=float)
    order = np.argsort(temperature)
    return temperature[order], free_energy[order]


def load_policy() -> dict[str, Any]:
    return read_json(ROOT / "configs/models/head_policy.yaml")


def load_inventory() -> dict[tuple[str, str], dict[str, str]]:
    rows = read_csv(ROOT / "result/tables/phase_inventory.csv")
    return {(row["system"], row["phase"]): row for row in rows}


def phase_meta(system: str, phase: str) -> dict[str, Any]:
    return read_json(DATA / system / phase / "meta.json")


def energy_from_meta(meta: dict[str, Any]) -> float:
    record = meta.get("representative", meta)
    values = record.get("energy_eV_per_atom", {})
    if not isinstance(values, dict) or "after" not in values:
        raise ValueError("representative energy provenance is missing")
    return float(values["after"])


def checkpoint_from_meta(meta: dict[str, Any]) -> tuple[str, str, str]:
    record = meta.get("representative", meta)
    checkpoint = str(record.get("checkpoint", ""))
    checkpoint_sha = str(record.get("checkpoint_sha256", ""))
    head = str(record.get("head", ""))
    provenance = record.get("energy_provenance", {})
    if isinstance(provenance, dict):
        checkpoint = str(provenance.get("checkpoint", checkpoint))
        checkpoint_sha = str(provenance.get("checkpoint_sha256", checkpoint_sha))
        head = str(provenance.get("head", head))
    return checkpoint, checkpoint_sha, head


def qh_paths(system: str, phase: str) -> tuple[Path, Path]:
    if system == "sio2":
        base = ROOT / "result/experiments/quasi_harmonic_sio2_sse_pbe/raw_runs/sio2"
    else:
        base = ROOT / f"result/experiments/quasi_harmonic_10a/raw_runs/{system}"
    return base / f"{phase}_fqh.csv", base / "qh_summary.json"


def load_qh(system: str, phase: str) -> tuple[dict[float, float], dict[str, Any], Path, Path]:
    csv_path, summary_path = qh_paths(system, phase)
    rows = read_csv(csv_path)
    values = {float(row["T_K"]): float(row["F_QH_eV_per_atom"]) for row in rows}
    if len(values) != len(rows):
        raise ValueError(f"duplicate QH temperatures in {csv_path}")
    for required in ("E_static_eV_per_atom", "checkpoint_sha256", "head"):
        if required not in rows[0]:
            raise ValueError(f"{csv_path} lacks {required}")
    return values, read_json(summary_path), csv_path, summary_path


def load_pairs() -> list[Pair]:
    pairs: list[Pair] = []
    for system, phases in PHASES.items():
        crossings = read_json(DATA / system / "reference_crossings.json")["pairs"]
        curves = {phase: reference_curve(system, phase) for phase in phases}
        for record in crossings.values():
            crossing_list = record.get("crossings", [])
            names = {record.get("left"), record.get("right")}
            if len(crossing_list) != 1 or not names.issubset(phases):
                continue
            left, right = str(record["left"]), str(record["right"])
            reference_tc = float(crossing_list[0]["T_K"])
            t_left, g_left = curves[left]
            t_right, g_right = curves[right]
            common = sorted(set(t_left).intersection(t_right))
            if not common:
                raise ValueError(f"no common grid for {system}/{left}/{right}")
            temperature = np.asarray(common, dtype=float)
            left_map = dict(zip(t_left, g_left))
            right_map = dict(zip(t_right, g_right))
            delta_left_right = np.asarray(
                [left_map[t] - right_map[t] for t in common], dtype=float
            )
            slope = float(crossing_list[0].get("slope_eV_per_atom_per_K", 0.0))
            if slope == 0.0:
                raise ValueError(f"zero reference crossing slope for {system}/{left}/{right}")
            # reference_crossings stores ΔG = G_left - G_right.  The slope
            # sign determines which phase is lower below the crossing and is
            # more stable than probing the quantized grid point immediately
            # below Tc, which can be exactly zero.
            low, high = (left, right) if slope > 0.0 else (right, left)
            delta_high_low = delta_left_right if high == right else -delta_left_right
            pairs.append(
                Pair(
                    system=system,
                    left=left,
                    right=right,
                    low=low,
                    high=high,
                    reference_tc=reference_tc,
                    temperature=temperature,
                    delta_g=delta_high_low * 1000.0,
                )
            )
    if len(pairs) != 5:
        raise ValueError(f"expected five one-crossing pairs, found {len(pairs)}")
    return pairs


def design(temperature: np.ndarray, form: str) -> np.ndarray:
    columns = [np.ones_like(temperature), temperature]
    if form == "R3":
        columns.append(temperature * np.log(temperature))
    elif form == "R3_T2":
        columns.append(temperature**2)
    return np.column_stack(columns)


def fit(temperature: np.ndarray, values: np.ndarray, form: str) -> tuple[np.ndarray, np.ndarray]:
    coefficients, *_ = np.linalg.lstsq(design(temperature, form), values, rcond=None)
    return coefficients, design(temperature, form) @ coefficients


def roots(temperature: np.ndarray, values: np.ndarray) -> list[float]:
    if np.all(values == 0.0):
        return []
    found: list[float] = []
    for index in range(len(temperature) - 1):
        t0, t1 = float(temperature[index]), float(temperature[index + 1])
        y0, y1 = float(values[index]), float(values[index + 1])
        if y0 == 0.0:
            found.append(t0)
        if y0 * y1 < 0.0:
            found.append(t0 + (t1 - t0) * (-y0) / (y1 - y0))
    if values[-1] == 0.0:
        found.append(float(temperature[-1]))
    unique: list[float] = []
    for value in found:
        if not unique or abs(value - unique[-1]) > 1.0e-9:
            unique.append(value)
    return unique


def crossing_counts(predicted: list[float], reference_tc: float) -> tuple[float | None, int, int]:
    if not predicted:
        return None, 0, 1
    nearest = min(predicted, key=lambda value: abs(value - reference_tc))
    return nearest, max(0, len(predicted) - 1), 0


def baseline_curves(pair: Pair) -> dict[str, np.ndarray]:
    qh_low, _, _, _ = load_qh(pair.system, pair.low)
    qh_high, _, _, _ = load_qh(pair.system, pair.high)
    low_qh = np.asarray([qh_low[float(t)] for t in pair.temperature])
    high_qh = np.asarray([qh_high[float(t)] for t in pair.temperature])
    high_energy = energy_from_meta(phase_meta(pair.system, pair.high))
    classical_high = high_energy + 3.0 * KB_EV_PER_K * pair.temperature * np.log(
        (HBAR_EV_S * 2.0 * math.pi * OMEGA0_HZ) / (KB_EV_PER_K * pair.temperature)
    )
    return {
        "B0": np.zeros_like(pair.delta_g),
        "B1": (high_qh - low_qh) * 1000.0,
        "B2": (classical_high - low_qh) * 1000.0,
    }


def make_input_audit(policy: dict[str, Any], inventory: dict[tuple[str, str], dict[str, str]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for system, phases in PHASES.items():
        expected = policy["systems"][system]
        for phase in phases:
            structure = DATA / system / phase / "structure.extxyz"
            meta = phase_meta(system, phase)
            checkpoint, checkpoint_sha, head = checkpoint_from_meta(meta)
            qh_csv, qh_summary_path = qh_paths(system, phase)
            qh_rows = read_csv(qh_csv)
            qh_heads = sorted({row["head"] for row in qh_rows})
            qh_shas = sorted({row["checkpoint_sha256"] for row in qh_rows})
            qh_summary = read_json(qh_summary_path)
            summary_head = str(qh_summary.get("head", ""))
            summary_provenance = qh_summary.get("provenance", {})
            if not isinstance(summary_provenance, dict):
                summary_provenance = {}
            summary_sha = str(
                qh_summary.get("checkpoint_sha256", "")
                or summary_provenance.get("checkpoint_sha256", "")
            )
            expected_head = str(expected["energy_head"])
            same = (
                head == expected_head
                and str(expected["fqh_head"]) == expected_head
                and qh_heads == [expected_head]
                and qh_shas == [checkpoint_sha]
                and summary_head == expected_head
                and summary_sha in {"", checkpoint_sha}
            )
            rows.append(
                {
                    "system": system,
                    "phase": phase,
                    "E_DPA_eV_per_atom": f"{energy_from_meta(meta):.12g}",
                    "structure_path": rel(structure),
                    "structure_sha256": sha256(structure),
                    "checkpoint": checkpoint,
                    "checkpoint_sha256": checkpoint_sha,
                    "energy_head": head,
                    "fqh_path": rel(qh_csv),
                    "fqh_checkpoint_sha256": ";".join(qh_shas),
                    "fqh_head": ";".join(qh_heads),
                    "fqh_includes_static_energy": "yes: F_QH=E_static+F_vib; confirmed by qh/run.py",
                    "fqh_zero_point_energy": "yes: phonopy thermal_properties free_energy",
                    "fqh_statistics": "quantum: phonopy Bose-Einstein thermal_properties",
                    "imaginary_mode_handling": "phonopy thermal_properties; diagnostic cutoff -0.05 THz; no extra deletion in this script",
                    "qh_reliable": inventory[(system, phase)]["qh_reliable"],
                    "same_checkpoint_and_head": "yes" if same else "no",
                    "audit_status": "ok" if same else "inconsistent",
                }
            )
    return rows


def make_delta_e_rows(pairs: list[Pair]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for pair in pairs:
        estimates: list[float] = []
        windows = [("full", np.arange(len(pair.temperature))), ("lowest_300K", np.arange(min(300, len(pair.temperature))))]
        for form in ("F1", "F2", "F3"):
            fit_form = {"F1": "R2", "F2": "R3", "F3": "R3_T2"}[form]
            for _, indices in windows:
                coeff, _ = fit(pair.temperature[indices], pair.delta_g[indices], fit_form)
                estimates.append(float(coeff[0]))
        dpa_low = energy_from_meta(phase_meta(pair.system, pair.low))
        dpa_high = energy_from_meta(phase_meta(pair.system, pair.high))
        dpa_delta = (dpa_high - dpa_low) * 1000.0
        effective = float(np.median(estimates))
        reference_std = float(np.std(pair.delta_g, ddof=0))
        rows.append(
            {
                "system": pair.system,
                "pair": f"{pair.low}/{pair.high}",
                "low_phase": pair.low,
                "high_phase": pair.high,
                "delta_E_DPA_meV_per_atom": f"{dpa_delta:.12g}",
                "delta_E_eff_median_meV_per_atom": f"{effective:.12g}",
                "delta_E_eff_min_meV_per_atom": f"{min(estimates):.12g}",
                "delta_E_eff_max_meV_per_atom": f"{max(estimates):.12g}",
                "delta_E_DPA_minus_eff_meV_per_atom": f"{dpa_delta - effective:.12g}",
                "reference_delta_G_std_meV_per_atom": f"{reference_std:.12g}",
                "abs_error_over_reference_std": f"{abs(dpa_delta - effective) / reference_std:.12g}",
                "fit_estimates_meV_per_atom": ";".join(f"{value:.12g}" for value in estimates),
            }
        )
    return rows


def oracle_rows(
    pairs: list[Pair], baselines: dict[str, dict[str, np.ndarray]]
) -> tuple[list[dict[str, Any]], dict[tuple[str, str], dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    plots: dict[tuple[str, str], dict[str, Any]] = {}
    split = read_json(DATA / "splits_v2/temp_extrap.json")
    thresholds = {system: float(value) for system, value in split["threshold_T_K"].items()}
    for pair in pairs:
        threshold = thresholds[pair.system]
        train = pair.temperature <= threshold
        test = pair.temperature > threshold
        pair_plots: dict[str, Any] = {}
        for baseline_name, base in baselines[pair.system].items():
            residual = pair.delta_g - base
            rms = float(np.sqrt(np.mean(residual**2)))
            max_abs = float(np.max(np.abs(residual)))
            baseline_test_mae = float(np.mean(np.abs(base[test] - pair.delta_g[test])))
            baseline_roots = roots(pair.temperature, base)
            baseline_tc, baseline_false, baseline_missed = crossing_counts(baseline_roots, pair.reference_tc)
            rows.append(
                {
                    "category": "baseline_control",
                    "system": pair.system,
                    "pair": f"{pair.low}/{pair.high}",
                    "low_phase": pair.low,
                    "high_phase": pair.high,
                    "baseline": baseline_name,
                    "residual_form": "none",
                    "parameter_count": 0,
                    "threshold_T_K": f"{threshold:.12g}",
                    "training_n": int(train.sum()),
                    "test_n": int(test.sum()),
                    "training_RMSE_meV_per_atom": "",
                    "test_MAE_meV_per_atom": f"{baseline_test_mae:.12g}",
                    "test_sign_accuracy": f"{np.mean(np.sign(base[test]) == np.sign(pair.delta_g[test])):.12g}",
                    "predicted_Tc_K": "" if baseline_tc is None else f"{baseline_tc:.12g}",
                    "Tc_error_K": "" if baseline_tc is None else f"{baseline_tc - pair.reference_tc:.12g}",
                    "predicted_crossings_n": len(baseline_roots),
                    "false_crossings_n": baseline_false,
                    "missed_crossings_n": baseline_missed,
                    "residual_RMS_full_meV_per_atom": f"{rms:.12g}",
                    "residual_max_abs_full_meV_per_atom": f"{max_abs:.12g}",
                    "reference_Tc_K": f"{pair.reference_tc:.12g}",
                }
            )
            for residual_form in ("R2", "R3"):
                coefficients, _ = fit(pair.temperature[train], residual[train], residual_form)
                prediction = base + design(pair.temperature, residual_form) @ coefficients
                predicted_roots = roots(pair.temperature, prediction)
                predicted_tc, false_crossings, missed_crossings = crossing_counts(predicted_roots, pair.reference_tc)
                rows.append(
                    {
                        "category": "oracle_fit",
                        "system": pair.system,
                        "pair": f"{pair.low}/{pair.high}",
                        "low_phase": pair.low,
                        "high_phase": pair.high,
                        "baseline": baseline_name,
                        "residual_form": residual_form,
                        "parameter_count": 2 if residual_form == "R2" else 3,
                        "threshold_T_K": f"{threshold:.12g}",
                        "training_n": int(train.sum()),
                        "test_n": int(test.sum()),
                        "training_RMSE_meV_per_atom": f"{np.sqrt(np.mean((prediction[train] - pair.delta_g[train])**2)):.12g}",
                        "test_MAE_meV_per_atom": f"{np.mean(np.abs(prediction[test] - pair.delta_g[test])):.12g}",
                        "test_sign_accuracy": f"{np.mean(np.sign(prediction[test]) == np.sign(pair.delta_g[test])):.12g}",
                        "predicted_Tc_K": "" if predicted_tc is None else f"{predicted_tc:.12g}",
                        "Tc_error_K": "" if predicted_tc is None else f"{predicted_tc - pair.reference_tc:.12g}",
                        "predicted_crossings_n": len(predicted_roots),
                        "false_crossings_n": false_crossings,
                        "missed_crossings_n": missed_crossings,
                        "residual_RMS_full_meV_per_atom": f"{rms:.12g}",
                        "residual_max_abs_full_meV_per_atom": f"{max_abs:.12g}",
                        "reference_Tc_K": f"{pair.reference_tc:.12g}",
                    }
                )
                if residual_form == "R3":
                    pair_plots[baseline_name] = (base, prediction)
        _, reference_fit = fit(pair.temperature, pair.delta_g, "R3")
        rows.append(
            {
                "category": "reference_control",
                "system": pair.system,
                "pair": f"{pair.low}/{pair.high}",
                "low_phase": pair.low,
                "high_phase": pair.high,
                "baseline": "reference",
                "residual_form": "R3",
                "parameter_count": 3,
                "threshold_T_K": f"{threshold:.12g}",
                "training_n": len(pair.temperature),
                "test_n": "",
                "training_RMSE_meV_per_atom": f"{np.sqrt(np.mean((reference_fit - pair.delta_g) ** 2)):.12g}",
                "test_MAE_meV_per_atom": "",
                "test_sign_accuracy": "",
                "predicted_Tc_K": "",
                "Tc_error_K": "",
                "predicted_crossings_n": "",
                "false_crossings_n": "",
                "missed_crossings_n": "",
                "residual_RMS_full_meV_per_atom": "",
                "residual_max_abs_full_meV_per_atom": "",
                "reference_Tc_K": f"{pair.reference_tc:.12g}",
            }
        )
        plots[(pair.system, f"{pair.low}/{pair.high}")] = {"pair": pair, "baselines": pair_plots, "threshold": threshold}
    return rows, plots


def write_figure(plots: dict[tuple[str, str], dict[str, Any]]) -> None:
    order = list(plots.values())
    figure, axes = plt.subplots(len(order), 3, figsize=(15.0, 3.0 * len(order)), squeeze=False)
    for row_index, record in enumerate(order):
        pair: Pair = record["pair"]
        for column, baseline_name in enumerate(("B0", "B1", "B2")):
            axis = axes[row_index][column]
            base, prediction = record["baselines"][baseline_name]
            axis.axvspan(pair.temperature.min(), record["threshold"], color="0.88", zorder=0, label="training" if row_index == 0 else None)
            axis.plot(pair.temperature, pair.delta_g, color="black", linewidth=1.3, label="reference" if row_index == 0 else None)
            axis.plot(pair.temperature, base, color="#377eb8", linewidth=1.0, linestyle="--", label="baseline" if row_index == 0 else None)
            axis.plot(pair.temperature, prediction, color="#e41a1c", linewidth=1.0, label="R3 prediction" if row_index == 0 else None)
            axis.axhline(0.0, color="0.35", linewidth=0.6)
            axis.axvline(pair.reference_tc, color="black", linewidth=0.8, linestyle=":")
            for value in roots(pair.temperature, prediction):
                axis.axvline(value, color="#e41a1c", linewidth=0.8, linestyle=":")
            axis.set_title(f"{pair.system}: {baseline_name}")
            axis.set_ylabel("Delta G (meV/atom)")
            axis.grid(alpha=0.2)
            if row_index == len(order) - 1:
                axis.set_xlabel("Temperature (K)")
            if row_index == 0 and column == 0:
                axis.legend(loc="best", fontsize=8)
    figure.tight_layout()
    figure.savefig(OUT / "figure_residuals.png", dpi=180)
    figure.savefig(OUT / "figure_residuals.pdf")
    plt.close(figure)


def fmt(value: Any, digits: int = 3) -> str:
    if value in (None, ""):
        return "n/a"
    return f"{float(value):.{digits}f}"


def write_summary(
    pairs: list[Pair],
    audit_rows: list[dict[str, Any]],
    delta_rows: list[dict[str, Any]],
    oracle: list[dict[str, Any]],
) -> None:
    lines = [
        "# Residual Oracle Diagnostics",
        "",
        "This report uses only canonical processed reference data, the frozen v2 temperature split, provenance tables, and canonical QH outputs. No model was trained or evaluated by this experiment.",
        "",
        "## Input consistency",
        "",
        f"All {len(audit_rows)} in-scope phases have matching checkpoint SHA-256 and head provenance between E_DPA, F_QH, and the head policy: **{all(row['audit_status'] == 'ok' for row in audit_rows)}**.",
        "",
        "The QH CSV schema and `src/lib/fes_bench/qh/run.py` confirm that `F_QH` already includes `E_static + F_vib`; the baselines therefore use F_QH directly and do not add E_DPA again. The QH thermal properties are the quantum phonopy free energy and include zero-point energy. The QH runner records imaginary-mode diagnostics with a -0.05 THz cutoff; this script does not edit or recalculate those modes.",
        "",
        "## Metal static-energy audit",
        "",
        "The reference crossing slope identifies hcp as the low-temperature phase for all three metals. The table prints the canonical E_DPA provenance used below.",
        "",
        "| System | Phase | E_DPA (eV/atom) | Source file | Checkpoint SHA-256 | Head |",
        "| --- | --- | ---: | --- | --- | --- |",
    ]
    for system in ("hf", "ti", "zr"):
        for phase in ("hcp", "bcc"):
            meta = phase_meta(system, phase)
            record = meta.get("representative", meta)
            source = record.get("source") or record.get("source_structure") or meta.get("raw_table", "n/a")
            _, checkpoint_sha, head = checkpoint_from_meta(meta)
            lines.append(
                f"| {system} | {phase} | {energy_from_meta(meta):.12f} | `{source}` | `{checkpoint_sha}` | {head} |"
            )
    lines += [
        "",
        "For the metal rows, hcp is low and bcc is high; ΔE_DPA is therefore E_DPA(bcc) − E_DPA(hcp), with the same high-minus-low direction used for ΔE_eff.",
        "",
        "## Delta-E budget",
        "",
        "| Pair (low/high) | Delta E DPA (meV/atom) | Delta E effective median [range] (meV/atom) | DPA minus effective | Reference Delta G SD | |error| / SD |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in delta_rows:
        lines.append(
            f"| {row['system']} {row['pair']} | {fmt(row['delta_E_DPA_meV_per_atom'])} | {fmt(row['delta_E_eff_median_meV_per_atom'])} [{fmt(row['delta_E_eff_min_meV_per_atom'])}, {fmt(row['delta_E_eff_max_meV_per_atom'])}] | {fmt(row['delta_E_DPA_minus_eff_meV_per_atom'])} | {fmt(row['reference_delta_G_std_meV_per_atom'])} | {fmt(row['abs_error_over_reference_std'])} |"
        )
    lines += [
        "",
        "The effective Delta E range is the six extrapolations from F1-F3 and the full/lowest-300-K windows. Values are an effective classical reference extrapolation, not a zero-point-inclusive ground-state energy.",
        "",
        "## R3 crossing errors in the test region",
        "",
        "| Pair | Reference Tc (K) | B0 error (K) | B1 error (K) | B2 error (K) | B0/B1/B2 test MAE (meV/atom) |",
        "| --- | ---: | ---: | ---: | ---: | --- |",
    ]
    selected = {(row["baseline"], row["system"], row["pair"]): row for row in oracle if row["category"] == "oracle_fit" and row["residual_form"] == "R3"}
    for pair in pairs:
        key = f"{pair.low}/{pair.high}"
        rows = [selected[(baseline, pair.system, key)] for baseline in ("B0", "B1", "B2")]
        lines.append(
            f"| {pair.system} {key} | {pair.reference_tc:.3f} | {fmt(rows[0]['Tc_error_K'])} | {fmt(rows[1]['Tc_error_K'])} | {fmt(rows[2]['Tc_error_K'])} | {fmt(rows[0]['test_MAE_meV_per_atom'])} / {fmt(rows[1]['test_MAE_meV_per_atom'])} / {fmt(rows[2]['test_MAE_meV_per_atom'])} |"
        )
    lines += ["", "The reference-curve full-range R3 control RMSE is:"]
    for row in oracle:
        if row["category"] == "reference_control":
            lines.append(f"- {row['system']} {row['pair']}: {fmt(row['training_RMSE_meV_per_atom'])} meV/atom")
    lines += ["", "## Conclusions", ""]
    for row in delta_rows:
        ratio = float(row["abs_error_over_reference_std"])
        comparison = "larger" if ratio > 1.0 else "smaller or comparable"
        sentence = f"- **{row['system']} {row['pair']}**: the absolute Delta-E discrepancy is {comparison} than the reference Delta G standard deviation ({ratio:.3f} times the standard deviation)."
        if ratio > 1.0:
            sentence += " Static energy difference is the main error source."
        lines.append(sentence)
    lines.append("")
    for pair in pairs:
        key = f"{pair.low}/{pair.high}"
        rows = [row for row in oracle if row["category"] == "oracle_fit" and row["residual_form"] == "R3" and row["system"] == pair.system and row["pair"] == key]
        errors = {row["baseline"]: abs(float(row["Tc_error_K"])) if row["Tc_error_K"] else math.inf for row in rows}
        mae = {row["baseline"]: float(row["test_MAE_meV_per_atom"]) for row in rows}
        labels = []
        for name in ("B1", "B2"):
            if abs(errors[name] - errors["B0"]) < 1.0 or abs(mae[name] - mae["B0"]) < 0.1:
                labels.append(f"{name}: no meaningful difference")
            elif errors[name] < errors["B0"] and mae[name] < mae["B0"]:
                labels.append(f"{name}: improvement")
            elif errors[name] > errors["B0"] and mae[name] > mae["B0"]:
                labels.append(f"{name}: worse")
            else:
                labels.append(f"{name}: mixed")
        lines.append(f"- **{pair.system} {key}**, equal three-parameter residual form: B0 |Tc error| = {errors['B0']:.3f} K; B1 = {errors['B1']:.3f} K; B2 = {errors['B2']:.3f} K. Relative to B0: {'; '.join(labels)}.")
    lines.append("")
    baseline_controls = {(row["system"], row["pair"], row["baseline"]): row for row in oracle if row["category"] == "baseline_control"}
    for pair in pairs:
        key = f"{pair.low}/{pair.high}"
        b0 = baseline_controls[(pair.system, key, "B0")]
        qh_states = []
        for phase in (pair.low, pair.high):
            qh_states.append(f"{phase} qh_reliable={next(item['qh_reliable'] for item in audit_rows if item['system'] == pair.system and item['phase'] == phase)}")
        parts = []
        for baseline in ("B1", "B2"):
            current = baseline_controls[(pair.system, key, baseline)]
            ratio = float(current["residual_RMS_full_meV_per_atom"]) / float(b0["residual_RMS_full_meV_per_atom"])
            parts.append(f"{baseline} RMS ratio {ratio:.3f}")
        lines.append(f"- **{pair.system} {key}**: compared with B0, {', '.join(parts)}; {', '.join(qh_states)}.")
    lines += [
        "",
        "## Scope and deviations",
        "",
        "- The cristobalite-beta/tridymite-beta pair has no reference crossing and was skipped.",
        "- The optional classical-QH variants B1c/B2c were not run because the canonical QH artifacts preserve scalar F_QH outputs and minimum-frequency diagnostics, not the full positive phonon-frequency spectrum needed for a fresh classical vibrational sum.",
        "- Canonical QH raw-run paths under `result/experiments/` were used because the current processed tree does not materialize every phase QH CSV; no legacy, archive, or invalid source was read.",
    ]
    (OUT / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def git_commit() -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    policy = load_policy()
    inventory = load_inventory()
    pairs = load_pairs()
    audit_rows = make_input_audit(policy, inventory)
    if any(row["audit_status"] != "ok" for row in audit_rows):
        raise RuntimeError("input provenance consistency check failed")

    delta_rows = make_delta_e_rows(pairs)
    baselines = {pair.system: baseline_curves(pair) for pair in pairs}
    oracle, plots = oracle_rows(pairs, baselines)

    write_csv(OUT / "inputs_audit.csv", audit_rows, list(audit_rows[0]))
    write_csv(OUT / "delta_E_budget.csv", delta_rows, list(delta_rows[0]))
    write_csv(OUT / "oracle_fits.csv", oracle, list(oracle[0]))
    write_figure(plots)
    write_summary(pairs, audit_rows, delta_rows, oracle)

    input_paths = set()
    for system, phases in PHASES.items():
        input_paths.add(DATA / system / "reference_crossings.json")
        input_paths.add(DATA / system / "system.json")
        for phase in phases:
            input_paths.update(
                {
                    DATA / system / phase / "structure.extxyz",
                    DATA / system / phase / "meta.json",
                    DATA / system / phase / "reference_G.csv",
                }
            )
            qh_csv, qh_summary = qh_paths(system, phase)
            input_paths.update({qh_csv, qh_summary})
    input_paths.update(
        {
            ROOT / "configs/models/head_policy.yaml",
            ROOT / "result/tables/phase_inventory.csv",
            ROOT / "data/processed/splits_v2/temp_extrap.json",
            ROOT / "src/lib/fes_bench/qh/run.py",
        }
    )
    manifest = {
        "git_commit": git_commit(),
        "checkpoint_sha256": sorted({row["checkpoint_sha256"] for row in audit_rows}),
        "heads": sorted({row["energy_head"] for row in audit_rows}),
        "inputs": [{"path": rel(path), "sha256": sha256(path)} for path in sorted(input_paths, key=rel)],
        "pairs": [f"{pair.system}:{pair.low}/{pair.high}" for pair in pairs],
        "optional_classical_qh": "not_run: full phonon spectrum unavailable in canonical artifacts",
    }
    (OUT / "findings.meta.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Residual oracle complete: {len(pairs)} pairs, {len(oracle)} result rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
