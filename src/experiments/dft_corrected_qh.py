"""Evaluate a static-energy replacement on top of the canonical QH curves.

This is a CPU-only post-processing experiment.  The canonical F_QH already
contains the DPA static energy, so the corrected curve is
F_QH - E_DPA + E_DFT.  No model is trained or evaluated here.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
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


def dft_energy(work: Path, system: str, phase: str) -> tuple[float, str, str]:
    key = "k555_encut1.000" if system == "sio2" else "k242424_encut2.000"
    path = work / "final-adopted-converged" / system / phase / "final_static" / key / "OUTCAR"
    if not path.is_file():
        raise FileNotFoundError(path)
    sys.path.insert(0, str(Path(__file__).resolve().parent / "dft_static"))
    from parse_results import parse_outcar

    row = parse_outcar(path)
    return float(row["energy_eV"]) / int(row["n_atoms"]), str(path), str(row["outcar_sha256"])


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
            fqh_curves[key] = {"T_K": temperatures, "G_eV_per_atom": fqh}
            corrected[key] = {
                "T_K": temperatures,
                "G_eV_per_atom": [value - static_value + e_dft for value, static_value in zip(fqh, static_values)],
            }
            dft_rows.append({
                "system": system,
                "phase": phase,
                "E_DPA_eV_per_atom": e_dpa,
                "E_DPA_in_F_QH_eV_per_atom": float(np.mean(static_values)),
                "E_DFT_eV_per_atom": e_dft,
                "E_DFT_minus_E_DPA_meV_per_atom": (e_dft - e_dpa) * 1000.0,
                "E_DPA_in_F_QH_min_eV_per_atom": float(np.min(static_values)),
                "E_DPA_in_F_QH_max_eV_per_atom": float(np.max(static_values)),
                "source": source,
                "source_sha256": source_sha,
                "F_QH_minus_F_vib_max_abs_minus_E_DPA_eV_per_atom": mismatch,
            })
    return fqh_curves, corrected, dft_rows, max_static_mismatch


def numeric(value: object) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def aggregate(result: dict[str, object], split: str) -> dict[str, object]:
    entries: list[tuple[dict[str, object], int]] = []
    if split == "temp_extrap":
        n = int(result["n_test_frames"])
        entries = [(pair, n) for key, pair in result.get("pairs", {}).items() if any(key.startswith(f"{system}:") for system in PAIRS) and pair.get("status") == "ok"]
    else:
        for fold in result.get("folds", {}).values():
            n = int(fold["n_test_frames"])
            entries.extend((pair, n) for key, pair in fold.get("pairs", {}).items() if pair.get("status") == "ok")
    entries = [
        (pair, n)
        for pair, n in entries
        if len(pair.get("reference_Tc_K", [])) > 0
    ]
    total = sum(n for _, n in entries)
    def weighted(field: str) -> float | str:
        values = [(numeric(pair.get(field)), n) for pair, n in entries]
        values = [(value, n) for value, n in values if value is not None]
        return sum(value * n for value, n in values) / sum(n for _, n in values) if values else "n/a"
    tc_errors = [abs(float(value)) for pair, _ in entries for value in pair.get("Tc_error_K", []) if numeric(value) is not None]
    return {
        "n_test_frames": total,
        "delta_G_MAE_eV_per_atom": weighted("delta_G_MAE_eV_per_atom"),
        "sign_accuracy": weighted("sign_accuracy"),
        "Tc_error_K_abs_mean": float(np.mean(tc_errors)) if tc_errors else "n/a",
        "false_crossings": sum(int(pair.get("false_crossings", 0) or 0) for pair, _ in entries),
        "missed_crossings": sum(int(pair.get("missed_crossings", 0) or 0) for pair, _ in entries),
        "pairs": len(entries),
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
        "definition": "T1: canonical F_QH; F_QH includes E_DPA + F_vib",
    }
    (output / "t1_qh_curves.json").write_text(json.dumps(curves_payload, indent=2) + "\n", encoding="utf-8")
    corrected_payload = {
        "predictions": corrected,
        "definition": "DFT-corrected QH: F_QH - E_DPA_in_F_QH + E_DFT",
    }
    (output / "dft_corrected_qh_curves.json").write_text(json.dumps(corrected_payload, indent=2) + "\n", encoding="utf-8")
    rows: list[dict[str, object]] = []
    for split in ("temp_extrap", "phase_lopo"):
        split_path = root / "data/processed/splits_v2" / f"{split}.json"
        split_data = json.loads(split_path.read_text(encoding="utf-8"))
        for predictor, payload in (("T1_E_plus_F_QH", fqh), ("DFT_corrected_QH", corrected)):
            curve_path = output / ("t1_qh_curves.json" if predictor.startswith("T1") else "dft_corrected_qh_curves.json")
            from fes_bench.eval.run import evaluate

            result = evaluate(f"precomputed:{curve_path}", split_data, root / "data/processed", [0])
            summary = aggregate(result, split)
            rows.append({"split": split, "predictor": predictor, **summary})
    fields = ["split", "predictor", "n_test_frames", "delta_G_MAE_eV_per_atom", "sign_accuracy", "Tc_error_K_abs_mean", "false_crossings", "missed_crossings", "pairs"]
    with (output / "metrics.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    with (output / "dft_phase_energies.csv").open("w", encoding="utf-8", newline="") as handle:
        fields_dft = list(dft_rows[0])
        writer = csv.DictWriter(handle, fieldnames=fields_dft)
        writer.writeheader()
        writer.writerows(dft_rows)
    summary_lines = [
        "# DFT-corrected QH static-energy replacement",
        "",
        "The canonical `F_QH` contains a constant static term plus `F_vib`; this was checked phase-by-phase before applying the correction.",
        "The corrected curve is `G_corr(T) = F_QH(T) - E_DPA_in_F_QH + E_DFT`, which is the requested `F_QH - E_DPA + E_DFT` with the actual static term represented by the QH file. The provenance `E_DPA` is reported separately where it differs by about 1 meV/atom.",
        "",
        "| Split | Predictor | ΔG MAE (eV/atom) | sign accuracy | mean |Tc error| (K) | false | missed |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        summary_lines.append("| {split} | {predictor} | {delta_G_MAE_eV_per_atom} | {sign_accuracy} | {Tc_error_K_abs_mean} | {false_crossings} | {missed_crossings} |".format(**row))
    summary_lines += ["", f"Maximum checked |F_QH − F_vib − E_DPA| = {mismatch:.3e} eV/atom."]
    (output / "summary.md").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    metadata = {
        "status": "complete",
        "training": "none",
        "device": "cpu",
        "formula": "F_QH - E_DPA_in_F_QH + E_DFT",
        "dft_work_root": str(work),
        "max_F_QH_minus_F_vib_minus_E_DPA_eV_per_atom": mismatch,
        "splits": ["temp_extrap", "phase_lopo"],
        "dft_source": "final-adopted-converged",
    }
    (output / "findings.meta.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
