"""Evaluate a static-energy replacement on top of the canonical QH curves.

This is a CPU-only post-processing experiment.  The canonical F_QH already
contains the DPA static energy, so the corrected curve is
F_QH - E_DPA + E_DFT.  No model is trained or evaluated here.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


SYSTEM_PHASES = {
    "hf": ("hcp", "bcc"),
    "ti": ("hcp", "bcc"),
    "zr": ("hcp", "bcc"),
    "sio2": ("quartz_beta", "cristobalite_beta", "tridymite_p63mmc"),
}
PAIRS = {
    "hf": (("hcp", "bcc"),),
    "ti": (("hcp", "bcc"),),
    "zr": (("hcp", "bcc"),),
    "sio2": (("quartz_beta", "cristobalite_beta"), ("quartz_beta", "tridymite_p63mmc")),
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


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


def dft_energy(work: Path, system: str, phase: str) -> tuple[float, Path, str]:
    key = "k555_encut1.000" if system == "sio2" else "k242424_encut2.000"
    path = work / "final-adopted-converged" / system / phase / "final_static" / key / "OUTCAR"
    if not path.is_file():
        raise FileNotFoundError(path)
    sys.path.insert(0, str(Path(__file__).resolve().parent / "dft_static"))
    from parse_results import parse_outcar

    row = parse_outcar(path)
    return float(row["energy_eV"]) / int(row["n_atoms"]), path, str(row["outcar_sha256"])


def curves(root: Path, work: Path) -> tuple[dict[str, object], dict[str, object], list[dict[str, object]], float]:
    dpa: dict[tuple[str, str], float] = {}
    fqh_curves: dict[str, dict[str, list[float]]] = {}
    corrected: dict[str, dict[str, list[float]]] = {}
    dft_rows: list[dict[str, object]] = []
    max_static_mismatch = 0.0
    for system, phases in SYSTEM_PHASES.items():
        for phase in phases:
            meta = json.loads((root / "data/processed" / system / phase / "meta.json").read_text(encoding="utf-8"))
            e_dpa = float(meta["representative"]["energy_eV_per_atom"]["after"])
            dpa[(system, phase)] = e_dpa
            rows = read_csv(root / "data/processed" / system / phase / "fqh.csv")
            e_dft, source, source_sha = dft_energy(work, system, phase)
            temperatures = [float(row["T_K"]) for row in rows]
            fqh = [float(row["F_QH_eV_per_atom"]) for row in rows]
            if "E_static_eV_per_atom" in rows[0]:
                static_values = np.asarray([float(row["E_static_eV_per_atom"]) for row in rows])
            else:
                fvib = [float(row["F_vib_eV_per_atom"]) for row in rows]
                static_values = np.asarray(fqh) - np.asarray(fvib)
            mismatch = float(np.max(np.abs(static_values - e_dpa)))
            max_static_mismatch = max(max_static_mismatch, mismatch)
            if float(np.ptp(static_values)) > 1.0e-8 and "E_static_eV_per_atom" not in rows[0]:
                raise AssertionError(f"{system}/{phase}: F_QH - F_vib is not a constant static term")
            key = f"{system}:{phase}"
            # Keep the repository's named T1 comparator exactly as defined:
            # E_DPA + F_QH.  The static-term audit below makes the resulting
            # duplicate-static convention explicit; DFT_corrected_QH removes
            # the term actually contained in F_QH before adding E_DFT.
            fqh_curves[key] = {"T_K": temperatures, "G_eV_per_atom": [e_dpa + value for value in fqh]}
            corrected[key] = {
                "T_K": temperatures,
                "G_eV_per_atom": [value - static_value + e_dft for value, static_value in zip(fqh, static_values)],
            }
            static_source = "E_static_eV_per_atom" if "E_static_eV_per_atom" in rows[0] else "F_QH_eV_per_atom - F_vib_eV_per_atom"
            dft_rows.append({
                "system": system,
                "phase": phase,
                "E_DPA_eV_per_atom": e_dpa,
                "E_DPA_in_F_QH_eV_per_atom": float(np.mean(static_values)),
                "E_DFT_eV_per_atom": e_dft,
                "E_DFT_minus_E_DPA_meV_per_atom": (e_dft - e_dpa) * 1000.0,
                "E_DPA_in_F_QH_min_eV_per_atom": float(np.min(static_values)),
                "E_DPA_in_F_QH_max_eV_per_atom": float(np.max(static_values)),
                "source": str(source.relative_to(work)),
                "source_sha256": source_sha,
                "static_term_source": static_source,
                "static_term_consistency": "match" if mismatch <= 1.0e-6 else "rowwise_qh_static_differs",
                "F_QH_minus_F_vib_max_abs_minus_E_DPA_eV_per_atom": mismatch,
            })
    return fqh_curves, corrected, dft_rows, max_static_mismatch


def numeric(value: object) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _fold_views(result: dict[str, object]) -> list[tuple[str, dict[str, object]]]:
    folds = result.get("folds")
    if isinstance(folds, dict):
        return [(str(name), view) for name, view in folds.items() if isinstance(view, dict)]
    return [("all", result)]


def _as_list(value: object) -> list[object]:
    return value if isinstance(value, list) else [value]


def pair_rows(result: dict[str, object], split: str, predictor: str) -> list[dict[str, object]]:
    """Return only the five reference-crossing pairs, retaining fold provenance."""
    rows: list[dict[str, object]] = []
    for fold_name, view in _fold_views(result):
        for pair_name, pair in view.get("pairs", {}).items():
            if not isinstance(pair, dict) or not pair.get("reference_Tc_K"):
                continue
            tc_error = _as_list(pair.get("Tc_error_K", "n/a:input_unavailable"))
            tc_from_dg = _as_list(pair.get("Tc_err_from_dG_K", "n/a:input_unavailable"))
            rows.append({
                "split": split,
                "predictor": predictor,
                "fold": fold_name,
                "pair": pair_name,
                "n_evaluation_points": int(pair.get("n_evaluation_points", view.get("n_test_frames", 0)) or 0),
                "delta_G_MAE_eV_per_atom": pair.get("delta_G_MAE_eV_per_atom", "n/a:input_unavailable"),
                "delta_G_RMSE_eV_per_atom": pair.get("delta_G_RMSE_eV_per_atom", "n/a:input_unavailable"),
                "sign_accuracy": pair.get("sign_accuracy", "n/a:input_unavailable"),
                "balanced_sign_accuracy": pair.get("balanced_sign_accuracy", "n/a:input_unavailable"),
                "reference_Tc_K": json.dumps(pair.get("reference_Tc_K", []), separators=(",", ":")),
                "predicted_Tc_K": json.dumps(pair.get("predicted_Tc_K", []), separators=(",", ":")),
                "Tc_error_K": json.dumps(tc_error, separators=(",", ":")),
                "Tc_err_from_dG_K": json.dumps(tc_from_dg, separators=(",", ":")),
                "false_crossings": pair.get("false_crossings", "n/a:input_unavailable"),
                "missed_crossings": pair.get("missed_crossings", "n/a:input_unavailable"),
            })
    return rows


def aggregate(rows: list[dict[str, object]]) -> dict[str, object]:
    weights = [int(row["n_evaluation_points"]) for row in rows]

    def weighted(field: str) -> float | str:
        values = [(numeric(row.get(field)), weight) for row, weight in zip(rows, weights)]
        values = [(value, weight) for value, weight in values if value is not None]
        return sum(value * weight for value, weight in values) / sum(weight for _, weight in values) if values else "n/a:input_unavailable"

    signed_tc: list[float] = []
    abs_tc: list[float] = []
    tc_markers: list[str] = []
    tc_dg_markers: list[str] = []
    for row in rows:
        for value in json.loads(str(row["Tc_error_K"])):
            if numeric(value) is not None:
                signed_tc.append(float(value)); abs_tc.append(abs(float(value)))
            elif isinstance(value, str) and value.startswith("n/a:"):
                tc_markers.append(value)
        for value in json.loads(str(row["Tc_err_from_dG_K"])):
            if isinstance(value, str) and value.startswith("n/a:"):
                tc_dg_markers.append(value)
    tc_marker = "n/a:missed_crossing" if "n/a:missed_crossing" in tc_markers else (tc_markers[0] if tc_markers else "n/a:input_unavailable")
    tc_dg_marker = "n/a:missed_crossing" if "n/a:missed_crossing" in tc_dg_markers else (tc_dg_markers[0] if tc_dg_markers else "n/a:input_unavailable")
    return {
        "n_evaluation_points": sum(weights),
        "n_pair_rows": len(rows),
        "delta_G_MAE_eV_per_atom": weighted("delta_G_MAE_eV_per_atom"),
        "delta_G_RMSE_eV_per_atom": weighted("delta_G_RMSE_eV_per_atom"),
        "sign_accuracy": weighted("sign_accuracy"),
        "balanced_sign_accuracy": weighted("balanced_sign_accuracy"),
        "Tc_error_K": float(np.mean(signed_tc)) if signed_tc else tc_marker,
        "Tc_error_K_abs_mean": float(np.mean(abs_tc)) if abs_tc else tc_marker,
        "Tc_err_from_dG_K": tc_dg_marker,
        "false_crossings": sum(int(row["false_crossings"]) for row in rows if numeric(row["false_crossings"]) is not None),
        "missed_crossings": sum(int(row["missed_crossings"]) for row in rows if numeric(row["missed_crossings"]) is not None),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--dft-work-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, default=Path("result/experiments/dft_corrected_qh"))
    args = parser.parse_args()
    root = args.repo_root.resolve()
    work = args.dft_work_root.resolve()
    output = (root / args.output_root).resolve()
    output.mkdir(parents=True, exist_ok=True)
    fqh, corrected, dft_rows, mismatch = curves(root, work)
    curves_payload = {
        "predictions": fqh,
        "definition": "T1 comparator: E_DPA + F_QH, with F_QH including the QH static term + F_vib",
    }
    (output / "t1_qh_curves.json").write_text(json.dumps(curves_payload, indent=2) + "\n", encoding="utf-8")
    corrected_payload = {
        "predictions": corrected,
        "definition": "DFT-corrected QH: F_QH - (static term contained in F_QH) + E_DFT",
    }
    (output / "dft_corrected_qh_curves.json").write_text(json.dumps(corrected_payload, indent=2) + "\n", encoding="utf-8")
    rows: list[dict[str, object]] = []
    pair_detail_rows: list[dict[str, object]] = []
    for split in ("temp_extrap", "phase_lopo"):
        split_path = root / "data/processed/splits_v2" / f"{split}.json"
        split_data = json.loads(split_path.read_text(encoding="utf-8"))
        for predictor, payload in (("T1_E_plus_F_QH", fqh), ("DFT_corrected_QH", corrected)):
            curve_path = output / ("t1_qh_curves.json" if predictor.startswith("T1") else "dft_corrected_qh_curves.json")
            from fes_bench.eval.run import evaluate

            result = evaluate(f"precomputed:{curve_path}", split_data, root / "data/processed", [0])
            details = pair_rows(result, split, predictor)
            pair_detail_rows.extend(details)
            summary = aggregate(details)
            rows.append({"split": split, "predictor": predictor, **summary})
    fields = ["split", "predictor", "n_evaluation_points", "n_pair_rows", "delta_G_MAE_eV_per_atom", "delta_G_RMSE_eV_per_atom", "sign_accuracy", "balanced_sign_accuracy", "Tc_error_K", "Tc_error_K_abs_mean", "Tc_err_from_dG_K", "false_crossings", "missed_crossings"]
    with (output / "metrics.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    pair_fields = list(pair_detail_rows[0]) if pair_detail_rows else []
    with (output / "pair_metrics.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=pair_fields)
        writer.writeheader()
        writer.writerows(pair_detail_rows)
    with (output / "dft_phase_energies.csv").open("w", encoding="utf-8", newline="") as handle:
        fields_dft = list(dft_rows[0])
        writer = csv.DictWriter(handle, fieldnames=fields_dft)
        writer.writeheader()
        writer.writerows(dft_rows)
    summary_lines = [
        "# DFT-corrected QH static-energy replacement",
        "",
        "The canonical `F_QH` static term was audited phase-by-phase before applying the correction. For SiO2 it is `F_QH - F_vib`; for metal QH files the row-wise `E_static_eV_per_atom` column is used. The correction therefore subtracts exactly the static term contained in `F_QH`, avoiding double counting.",
        "The T1 comparator is kept exactly as requested, `G_T1(T) = E_DPA + F_QH(T)`. Because the audit confirms that F_QH already contains its static term, this named comparator has the repository's legacy double-static convention; the corrected curve instead uses `G_corr(T) = F_QH(T) - static_term_contained_in_F_QH(T) + E_DFT`, with both sides of every pair treated identically. The five-pair summaries exclude cristobalite--tridymite because it has no reference crossing; pair-level output retains the five crossing pairs with fold provenance.",
        "",
        "| Split | Predictor | ΔG MAE (eV/atom) | balanced sign accuracy | Tc error | mean |Tc error| (K) | false | missed |",
        "|---|---|---:|---:|---|---:|---:|---:|",
    ]
    for row in rows:
        summary_lines.append("| {split} | {predictor} | {delta_G_MAE_eV_per_atom} | {balanced_sign_accuracy} | {Tc_error_K} | {Tc_error_K_abs_mean} | {false_crossings} | {missed_crossings} |".format(**row))
    summary_lines += ["", f"Maximum checked static-term difference relative to metadata E_DPA = {mismatch:.3e} eV/atom. The correction uses the exact QH-contained static term, not metadata E_DPA.", "", "DFT convergence status: all nine adopted static outputs are `final-adopted-converged`; the Hf neighboring-k-mesh deviation remains documented in the DFT static convergence report."]
    (output / "summary.md").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    input_paths = [root / "src/experiments/dft_corrected_qh.py"]
    for split in ("temp_extrap", "phase_lopo"):
        input_paths.append(root / "data/processed/splits_v2" / f"{split}.json")
    metadata = {
        "status": "complete",
        "training": "none",
        "device": "cpu",
        "formula": "G_corr(T) = F_QH(T) - static_term_contained_in_F_QH(T) + E_DFT",
        "t1_comparator_formula": "G_T1(T) = E_DPA + F_QH(T)",
        "dft_work_root": "external/dft_static_work/results",
        "max_F_QH_minus_F_vib_minus_E_DPA_eV_per_atom": mismatch,
        "splits": ["temp_extrap", "phase_lopo"],
        "dft_source": "final-adopted-converged",
        "reference_crossing_pair_count": 5,
        "pair_metrics": "pair_metrics.csv",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": git_commit(root),
        "script_sha256": sha256(root / "src/experiments/dft_corrected_qh.py"),
        "input_sha256": [{"path": str(path.relative_to(root)), "sha256": sha256(path)} for path in input_paths],
        "dft_source_sha256": {f"{row['system']}:{row['phase']}": row["source_sha256"] for row in dft_rows},
    }
    (output / "findings.meta.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
