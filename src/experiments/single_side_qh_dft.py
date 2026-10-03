"""DFT-static replacement and constant-sign control — plan experiment single_side_qh.

This script is a no-training extension of the single-sided QH experiment.  It
uses the canonical processed QH curves, replaces the DPA static term with the
adopted relaxed DFT static energy, and evaluates the constant-sign pair floor
on every frozen v2 split.  All outputs are provenance-bearing and no model is
run.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

import single_side_qh as base


SYSTEM_PHASES = {
    "hf": ("hcp", "bcc"),
    "ti": ("hcp", "bcc"),
    "zr": ("hcp", "bcc"),
    "sio2": ("quartz_beta", "cristobalite_beta", "tridymite_p63mmc"),
}
SEEDS = (11, 23, 37)
EPSILON_SIGN = -1.0e-9


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_commit(root: Path) -> str:
    try:
        return subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def load_split(root: Path, name: str) -> tuple[dict[str, Any], Path, str]:
    path = root / "data/processed/splits_v2" / f"{name}.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    return payload, path, sha256(path)


def reference_curve(root: Path, system: str, phase: str) -> tuple[np.ndarray, np.ndarray]:
    rows = read_csv(root / "data/processed" / system / phase / "reference_G.csv")
    return np.asarray([float(row["T_K"]) for row in rows]), np.asarray([float(row["G_eV_per_atom"]) for row in rows])


def fqh_curve(root: Path, system: str, phase: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    rows = read_csv(root / "data/processed" / system / phase / "fqh.csv")
    required = {"T_K", "F_QH_eV_per_atom"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"{system}/{phase} fqh.csv lacks {sorted(required)}")
    t = np.asarray([float(row["T_K"]) for row in rows])
    fqh = np.asarray([float(row["F_QH_eV_per_atom"]) for row in rows])
    if "E_static_eV_per_atom" in rows[0]:
        e_dpa_values = np.asarray([float(row["E_static_eV_per_atom"]) for row in rows])
    elif "F_vib_eV_per_atom" in rows[0]:
        e_dpa_values = fqh - np.asarray([float(row["F_vib_eV_per_atom"]) for row in rows])
    else:
        raise ValueError(f"{system}/{phase} fqh.csv has neither E_static nor F_vib")
    return t, fqh, e_dpa_values, float(np.ptp(e_dpa_values))


def load_dft_manifest(root: Path) -> dict[tuple[str, str], dict[str, str]]:
    path = root / "result/experiments/single_side_qh/dft_static_energy_manifest.csv"
    rows = read_csv(path)
    values = {(row["system"], row["phase"]): row for row in rows}
    expected = {(system, phase) for system, phases in SYSTEM_PHASES.items() for phase in phases}
    if set(values) != expected:
        raise AssertionError(f"DFT energy manifest keys differ: missing={sorted(expected - set(values))}")
    return values


def phase_curves(root: Path, manifest: dict[tuple[str, str], dict[str, str]]) -> tuple[dict[str, dict[float, float]], dict[str, dict[str, Any]]]:
    t1: dict[str, dict[float, float]] = {}
    dft: dict[str, dict[float, float]] = {}
    metadata: dict[str, dict[str, Any]] = {}
    for system, phases in SYSTEM_PHASES.items():
        for phase in phases:
            t_fqh, fqh, e_dpa, spread = fqh_curve(root, system, phase)
            t_ref, _ = reference_curve(root, system, phase)
            if not np.array_equal(t_fqh, t_ref):
                raise AssertionError(f"QH/reference temperature mismatch for {system}/{phase}")
            e_dft = float(manifest[(system, phase)]["E_DFT_eV_per_atom"])
            key = f"{system}:{phase}"
            t1[key] = {float(t): float(e + f) for t, e, f in zip(t_fqh, e_dpa, fqh)}
            dft[key] = {float(t): float(f - e + e_dft) for t, e, f in zip(t_fqh, e_dpa, fqh)}
            metadata[key] = {
                "E_DPA_eV_per_atom_from_F_QH_minus_F_vib_or_file": float(np.mean(e_dpa)),
                "F_QH_minus_F_vib_range_eV_per_atom": spread,
                "E_DFT_eV_per_atom": e_dft,
                "formula_T1": "E_DPA + F_QH",
                "formula_DFT_replaced": "F_QH - E_DPA + E_DFT = F_vib + E_DFT",
                "dft_source_path": manifest[(system, phase)]["source_path"],
                "dft_source_sha256": manifest[(system, phase)]["source_sha256"],
                "dft_convergence_status": manifest[(system, phase)]["convergence_status"],
            }
    return {"T1_E_plus_F_QH": t1, "DFT_static_replaced": dft}, metadata


def ordered_pairs(root: Path, system: str, phases: list[str]) -> list[tuple[str, str, str]]:
    payload = json.loads((root / "data/processed" / system / "reference_crossings.json").read_text(encoding="utf-8"))
    output = []
    for name, item in payload.get("pairs", {}).items():
        left, right = str(item["left"]), str(item["right"])
        if left in phases and right in phases:
            output.append((left, right, f"{system}:{name}"))
    return output


def roots(temperature: np.ndarray, values: np.ndarray) -> list[float]:
    found: list[float] = []
    for index in range(len(values) - 1):
        left, right = float(values[index]), float(values[index + 1])
        if left == 0.0:
            found.append(float(temperature[index]))
        if left * right < 0.0:
            found.append(float(temperature[index] - left * (temperature[index + 1] - temperature[index]) / (right - left)))
    if values[-1] == 0.0:
        found.append(float(temperature[-1]))
    unique: list[float] = []
    for value in found:
        if not unique or abs(value - unique[-1]) > 1.0e-8:
            unique.append(value)
    return unique


def balanced_sign(observed: np.ndarray, reference: np.ndarray) -> float:
    values: list[float] = []
    positive, negative = reference > 0.0, reference < 0.0
    if np.any(positive):
        values.append(float(np.mean(observed[positive] > 0.0)))
    if np.any(negative):
        values.append(float(np.mean(observed[negative] < 0.0)))
    return float(np.mean(values)) if values else 0.0


def pair_metric(root: Path, system: str, pair_name: str, temperature: np.ndarray, observed: np.ndarray) -> dict[str, Any]:
    left, right = pair_name.split("_minus_", 1)
    _, left_ref = reference_curve(root, system, left)
    _, right_ref = reference_curve(root, system, right)
    ref_t, _ = reference_curve(root, system, left)
    reference = np.interp(temperature, ref_t, left_ref - right_ref)
    crossing_payload = json.loads((root / "data/processed" / system / "reference_crossings.json").read_text(encoding="utf-8"))
    crossings = crossing_payload["pairs"].get(pair_name, {}).get("crossings", [])
    roots_ref = [float(item["T_K"]) for item in crossings]
    roots_pred = roots(temperature, observed)
    missed = max(0, len(roots_ref) - len(roots_pred))
    false = max(0, len(roots_pred) - len(roots_ref))
    if not roots_ref:
        tc_error: Any = "n/a:no_reference_crossing"
        tc_from_dg: Any = "n/a:no_reference_crossing"
    elif missed:
        tc_error = "n/a:missed_crossing"
        tc_from_dg = "n/a:missed_crossing"
    else:
        tc_error = float(roots_pred[0] - roots_ref[0])
        slope = float(crossings[0].get("slope_eV_per_atom_per_K", 0.0))
        tc_from_dg = float(abs(np.interp(roots_ref[0], temperature, observed - reference)) / slope) if slope else "n/a:zero_slope"
    return {
        "n_evaluation_points": len(temperature),
        "delta_G_MAE_eV_per_atom": float(np.mean(np.abs(observed - reference))),
        "delta_G_RMSE_eV_per_atom": float(np.sqrt(np.mean((observed - reference) ** 2))),
        "sign_accuracy": float(np.mean(np.sign(observed) == np.sign(reference))),
        "balanced_sign_accuracy": balanced_sign(observed, reference),
        "reference_Tc_K": roots_ref,
        "predicted_Tc_K": roots_pred,
        "Tc_error_K": tc_error,
        "Tc_err_from_dG_K": tc_from_dg,
        "false_crossings": false,
        "missed_crossings": missed,
    }


def split_rows(root: Path, split_name: str, split: dict[str, Any], curves: dict[str, dict[float, float]]) -> list[dict[str, Any]]:
    folds = split.get("folds", {"all": split})
    rows: list[dict[str, Any]] = []
    for fold_name, fold in folds.items():
        test: dict[tuple[str, str], set[int]] = defaultdict(set)
        train: dict[tuple[str, str], set[int]] = defaultdict(set)
        for item in fold.get("test", []):
            test[(str(item["system"]), str(item["phase"]))].add(int(item["T_index"]))
        for item in fold.get("train", []):
            train[(str(item["system"]), str(item["phase"]))].add(int(item["T_index"]))
        all_keys = set(test) | set(train)
        for system in sorted({key[0] for key in all_keys}):
            phases = [phase for candidate_system, phase in all_keys if candidate_system == system]
            for left, right, pair_key in ordered_pairs(root, system, phases):
                if (system, left) not in test and (system, right) not in test:
                    continue
                left_t, _ = reference_curve(root, system, left)
                right_t, _ = reference_curve(root, system, right)
                shared_test = test.get((system, left), set()) & test.get((system, right), set())
                shared = sorted(shared_test or ((train.get((system, left), set()) | test.get((system, left), set())) & (train.get((system, right), set()) | test.get((system, right), set()))))
                temperature = np.asarray([left_t[index] for index in shared], dtype=float)
                left_name, right_name = f"{system}:{left}", f"{system}:{right}"
                for predictor, phase_map in curves.items():
                    observed = np.asarray([phase_map[left_name][float(t)] - phase_map[right_name][float(t)] for t in temperature])
                    rows.append({"split": split_name, "fold": fold_name, "pair": pair_key, "predictor": predictor, "eval_scope": "test_shared" if shared_test else "test_plus_partner_train", **pair_metric(root, system, pair_key.split(":", 1)[1], temperature, observed)})
    return rows


def constant_sign_rows(root: Path, split_name: str, split: dict[str, Any]) -> list[dict[str, Any]]:
    curves = {f"{system}:{phase}": reference_curve(root, system, phase) for system, phases in SYSTEM_PHASES.items() for phase in phases}
    folds = split.get("folds", {"all": split})
    rows: list[dict[str, Any]] = []
    for fold_name, fold in folds.items():
        test: dict[tuple[str, str], set[int]] = defaultdict(set)
        train: dict[tuple[str, str], set[int]] = defaultdict(set)
        for item in fold.get("test", []):
            test[(str(item["system"]), str(item["phase"]))].add(int(item["T_index"]))
        for item in fold.get("train", []):
            train[(str(item["system"]), str(item["phase"]))].add(int(item["T_index"]))
        all_keys = set(test) | set(train)
        for system in sorted({key[0] for key in all_keys}):
            phases = [phase for candidate_system, phase in all_keys if candidate_system == system]
            for left, right, pair_key in ordered_pairs(root, system, phases):
                if (system, left) not in test and (system, right) not in test:
                    continue
                left_t, _ = curves[f"{system}:{left}"]
                shared_test = test.get((system, left), set()) & test.get((system, right), set())
                shared = sorted(shared_test or ((train.get((system, left), set()) | test.get((system, left), set())) & (train.get((system, right), set()) | test.get((system, right), set()))))
                temperature = np.asarray([left_t[index] for index in shared], dtype=float)
                observed = np.full(len(temperature), EPSILON_SIGN)
                metrics = pair_metric(root, system, pair_key.split(":", 1)[1], temperature, observed)
                rows.append({"split": split_name, "fold": fold_name, "pair": pair_key, "predictor": "constant_sign", "eval_scope": "test_shared" if shared_test else "test_plus_partner_train", **metrics})
    return rows


def summarize(rows: list[dict[str, Any]], *, exclude_no_reference: bool = False) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for split in sorted({str(row["split"]) for row in rows}):
        for predictor in sorted({str(row["predictor"]) for row in rows if row["split"] == split}):
            selected = [row for row in rows if row["split"] == split and row["predictor"] == predictor]
            if exclude_no_reference:
                selected = [row for row in selected if row.get("reference_Tc_K")]
            if not selected:
                continue
            weights = np.asarray([int(row["n_evaluation_points"]) for row in selected], dtype=float)
            def wmean(key: str) -> float:
                return float(np.average([float(row[key]) for row in selected], weights=weights))
            tc = [abs(float(row["Tc_error_K"])) for row in selected if isinstance(row["Tc_error_K"], (int, float))]
            output.append({"split": split, "predictor": predictor, "n_pair_rows": len(selected), "delta_G_MAE_eV_per_atom": wmean("delta_G_MAE_eV_per_atom"), "delta_G_RMSE_eV_per_atom": wmean("delta_G_RMSE_eV_per_atom"), "sign_accuracy": wmean("sign_accuracy"), "balanced_sign_accuracy": wmean("balanced_sign_accuracy"), "Tc_error_abs_mean_K": float(np.mean(tc)) if tc else None, "false_crossings": sum(int(row["false_crossings"]) for row in selected if isinstance(row["false_crossings"], (int, float))), "missed_crossings": sum(int(row["missed_crossings"]) for row in selected if isinstance(row["missed_crossings"], (int, float)))})
    return output


def update_readme(path: Path, dft_summary: list[dict[str, Any]], sign_summary: list[dict[str, Any]]) -> None:
    text = path.read_text(encoding="utf-8") if path.exists() else "# Single-sided QH diagnostic\n"
    marker = "\n## DFT static replacement and constant-sign control\n"
    if marker in text:
        text = text.split(marker, 1)[0]
    lines = [marker.rstrip(), "", "### DFT static replacement", "", "The canonical QH files contain the static term: for SiO2, `F_QH − F_vib` is constant; for metal QH files the row-wise `E_static_eV_per_atom` column is used because the minimized-volume static energy can vary with temperature. The replacement therefore uses `G_corr(T) = F_QH(T) − E_static,DPA(T) + E_DFT = F_vib(T) + E_DFT`, with both phases of each pair treated identically. `E_DFT` is the relaxed static `energy without entropy` from the `final-adopted-converged` OUTCAR manifest. The DFT output status is recorded per phase; the adopted production outputs are present, while the Hf neighboring-mesh deviation remains documented in the DFT convergence report.", "", "| Split | Predictor | ΔG MAE (meV/atom) | ΔG RMSE (meV/atom) | sign accuracy | balanced sign accuracy | mean |Tc error| (K) | false | missed |", "|---|---|---:|---:|---:|---:|---:|---:|---:|"]
    for row in dft_summary:
        lines.append(f"| {row['split']} | {row['predictor']} | {float(row['delta_G_MAE_eV_per_atom'])*1000:.3f} | {float(row['delta_G_RMSE_eV_per_atom'])*1000:.3f} | {float(row['sign_accuracy']):.4f} | {float(row['balanced_sign_accuracy']):.4f} | {row['Tc_error_abs_mean_K'] if row['Tc_error_abs_mean_K'] is not None else 'n/a'} | {row['false_crossings']} | {row['missed_crossings']} |")
    lines += ["", "### Constant-sign baseline", "", "`constant_sign` predicts the oriented pair ΔG as −1×10⁻⁹ eV/atom at every point, i.e. the high-temperature/right phase is always declared stable. The tiny nonzero magnitude avoids the exact-zero degenerate-prediction marker; this is a sign-only control, not an energy baseline. Balanced accuracy is the mean recall of the positive and negative reference sign classes present on the evaluated grid; for a one-sign pair the present class is reported.", "", "| Split | ΔG MAE (meV/atom) | ordinary sign accuracy | balanced sign accuracy |", "|---|---:|---:|---:|"]
    for row in sign_summary:
        lines.append(f"| {row['split']} | {float(row['delta_G_MAE_eV_per_atom'])*1000:.3f} | {float(row['sign_accuracy']):.4f} | {float(row['balanced_sign_accuracy']):.4f} |")
    lines += ["", "Pair-level rows, including fold scope and crossing markers, are in `dft_static_replacement_pair_metrics.csv` and `constant_sign_pair_metrics.csv`. Reproduce with:", "", "```bash", "PYTHONPATH=src python src/experiments/single_side_qh_dft.py --repo-root .", "```", ""]
    path.write_text(text.rstrip() + "\n" + "\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    args = parser.parse_args()
    root = args.repo_root.resolve()
    output = root / "result/experiments/single_side_qh"
    manifest = load_dft_manifest(root)
    curves, phase_meta = phase_curves(root, manifest)
    dft_rows: list[dict[str, Any]] = []
    sign_rows: list[dict[str, Any]] = []
    input_paths = [root / "result/experiments/single_side_qh/dft_static_energy_manifest.csv"]
    for split_name in ("temp_extrap", "phase_lopo", "system_loso"):
        split, split_path, _ = load_split(root, split_name)
        dft_rows.extend(split_rows(root, split_name, split, curves))
        sign_rows.extend(constant_sign_rows(root, split_name, split))
        if split_name == "system_loso":
            overlap_folds = {
                str(name): {"train": fold.get("train", []), "test": fold["overlap_T"]}
                for name, fold in split.get("folds", {}).items()
                if isinstance(fold, dict) and isinstance(fold.get("overlap_T"), list)
            }
            if overlap_folds:
                overlap_split = {"folds": overlap_folds}
                dft_rows.extend(split_rows(root, "system_loso_overlap_T", overlap_split, curves))
                sign_rows.extend(constant_sign_rows(root, "system_loso_overlap_T", overlap_split))
        input_paths.append(split_path)
    write_csv(output / "dft_static_replacement_pair_metrics.csv", dft_rows)
    write_csv(output / "constant_sign_pair_metrics.csv", sign_rows)
    # The DFT replacement report is defined on the five pairs with a frozen
    # reference crossing; the pair CSV still retains cristobalite--tridymite
    # with its explicit no-reference marker.
    dft_summary = summarize(dft_rows, exclude_no_reference=True)
    sign_summary = summarize(sign_rows)
    write_csv(output / "dft_static_replacement_summary.csv", dft_summary)
    write_csv(output / "constant_sign_summary.csv", sign_summary)
    (output / "dft_static_replaced_curves.json").write_text(json.dumps({"predictions": curves["DFT_static_replaced"]}, indent=2) + "\n", encoding="utf-8")
    (output / "dft_static_replacement.meta.json").write_text(json.dumps({
        "git_commit": git_commit(root),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "script": {"path": "src/experiments/single_side_qh_dft.py", "sha256": sha256(root / "src/experiments/single_side_qh_dft.py")},
        "inputs": [{"path": str(path.relative_to(root)), "sha256": sha256(path)} for path in input_paths],
        "phase_energy_audit": phase_meta,
        "dft_convergence_summary": "final-adopted-converged outputs are present for all nine phases; Hf neighboring-mesh criterion deviation is retained in result/experiments/dft_static/summary.md",
        "constant_sign_definition": "delta_G = -1e-9 eV/atom for every oriented pair and evaluation point",
        "no_training": True,
    }, indent=2) + "\n", encoding="utf-8")
    update_readme(output / "README.md", dft_summary, sign_summary)
    print(output / "dft_static_replacement_summary.csv")
    print(output / "constant_sign_summary.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
