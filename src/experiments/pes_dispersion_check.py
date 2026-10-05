"""Check analytic PBE-D3(BJ) corrections on SiO2 structures — plan PES step 1.

This script deliberately reads structure files and static-energy metadata only;
it refuses to run if a free-energy table is present in its input list.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from ase.io import read
from dftd3.interface import DispersionModel, RationalDampingParam


PHASES = ("quartz_beta", "cristobalite_beta", "tridymite_p63mmc")
PAIR_LABELS = {
    "quartz_minus_cristobalite": ("quartz_beta", "cristobalite_beta"),
    "quartz_minus_tridymite": ("quartz_beta", "tridymite_p63mmc"),
}
S8 = 0.7875
A1 = 0.4289
A2 = 4.4407
S9 = 0.0  # Pairwise VASP IVDW=12 convention; ATM is not enabled.


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def d3_energy(path: Path) -> tuple[float, int]:
    atoms = read(path)
    model = DispersionModel(
        np.asarray(atoms.numbers, dtype=int),
        np.asarray(atoms.positions, dtype=float),
        np.asarray(atoms.cell.array, dtype=float),
        np.asarray(atoms.pbc, dtype=bool),
    )
    parameter = RationalDampingParam(s6=1.0, s8=S8, s9=S9, a1=A1, a2=A2)
    result = model.get_dispersion(parameter, grad=False)
    return float(result["energy"]), len(atoms)


def dft_energy(root: Path, path: Path) -> float:
    import sys

    sys.path.insert(0, str(root / "src/experiments/dft_static"))
    from parse_results import parse_outcar

    row = parse_outcar(path)
    return float(row["energy_eV"]) / int(row["n_atoms"])


def read_snapshot_energy(path: Path) -> float:
    atoms = read(path)
    value = atoms.info.get("energy", atoms.info.get("free_energy"))
    if value is None:
        header = path.read_text(encoding="utf-8").splitlines()[1]
        match = re.search(r"(?:energy|free_energy)=([-+0-9.eE]+)", header)
        value = float(match.group(1)) if match else None
    if value is None:
        raise ValueError(f"snapshot has no static energy: {path}")
    return float(value) / len(atoms)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--dft-work-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, default=Path("result/experiments/pes_dispersion_check"))
    args = parser.parse_args()
    root = args.repo_root.resolve()
    work = args.dft_work_root.resolve()
    output = (root / args.output_root).resolve()
    output.mkdir(parents=True, exist_ok=True)
    input_paths: list[Path] = []
    rows: list[dict[str, object]] = []
    energy_sets: dict[str, dict[str, float]] = {"dpa_plus_d3": {}, "dft_final": {}, "dft_snapshot": {}}
    for phase in PHASES:
        phase_dir = root / "data/processed/sio2" / phase
        meta_path = phase_dir / "meta.json"
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        dpa = float(meta["representative"]["energy_eV_per_atom"]["after"])
        current_path = phase_dir / "structure.extxyz"
        dft_dir = work / "final-adopted-converged/sio2" / phase / "final_static/k555_encut1.000"
        dft_structure = dft_dir / "CONTCAR"
        dft_outcar = dft_dir / "OUTCAR"
        snapshot_path = root / "data/processed/sio2/_archive_md_snapshot" / phase / "structure.source.extxyz"
        current_d3, n_current = d3_energy(current_path)
        relaxed_d3, n_relaxed = d3_energy(dft_structure)
        snapshot_d3, n_snapshot = d3_energy(snapshot_path)
        e_dft = dft_energy(root, dft_outcar)
        e_snapshot = read_snapshot_energy(snapshot_path)
        energy_sets["dpa_plus_d3"][phase] = dpa + current_d3 / n_current
        energy_sets["dft_final"][phase] = e_dft
        energy_sets["dft_snapshot"][phase] = e_snapshot
        rows.extend([
            {"phase": phase, "structure_kind": "current_dpa", "structure": str(current_path.relative_to(root)), "n_atoms": n_current, "E_DPA_or_DFT_eV_per_atom": dpa, "E_D3_eV_per_atom": current_d3 / n_current, "E_DPA_plus_D3_eV_per_atom": dpa + current_d3 / n_current, "structure_sha256": sha256(current_path)},
            {"phase": phase, "structure_kind": "dft_relaxed", "structure": str(dft_structure), "n_atoms": n_relaxed, "E_DPA_or_DFT_eV_per_atom": e_dft, "E_D3_eV_per_atom": relaxed_d3 / n_relaxed, "E_DPA_plus_D3_eV_per_atom": "n/a:DFT_energy_is_already_PBE-D3BJ", "structure_sha256": sha256(dft_structure)},
            {"phase": phase, "structure_kind": "reference_snapshot", "structure": str(snapshot_path.relative_to(root)), "n_atoms": n_snapshot, "E_DPA_or_DFT_eV_per_atom": e_snapshot, "E_D3_eV_per_atom": snapshot_d3 / n_snapshot, "E_DPA_plus_D3_eV_per_atom": "n/a:DFT_energy_is_already_PBE-D3BJ", "structure_sha256": sha256(snapshot_path)},
        ])
        input_paths.extend([meta_path, current_path, dft_structure, dft_outcar, snapshot_path])
    if any("reference_G.csv" in str(path) for path in input_paths):
        raise AssertionError("dispersion check must not read free-energy tables")
    summary_rows = []
    for label, (left, right) in PAIR_LABELS.items():
        model = energy_sets["dpa_plus_d3"][left] - energy_sets["dpa_plus_d3"][right]
        final = energy_sets["dft_final"][left] - energy_sets["dft_final"][right]
        snapshot = energy_sets["dft_snapshot"][left] - energy_sets["dft_snapshot"][right]
        summary_rows.append({"pair": label, "delta_E_DPA_plus_D3_eV_per_atom": model, "delta_E_DFT_final_eV_per_atom": final, "model_minus_DFT_meV_per_atom": (model - final) * 1000.0, "delta_E_DFT_snapshot_eV_per_atom": snapshot, "model_minus_snapshot_meV_per_atom": (model - snapshot) * 1000.0})
    fields = list(rows[0])
    with (output / "phase_results.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader(); writer.writerows(rows)
    with (output / "pair_summary.csv").open("w", encoding="utf-8", newline="") as handle:
        fields = list(summary_rows[0]); writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader(); writer.writerows(summary_rows)
    lines = [
        "# SiO2 analytic D3(BJ) dispersion check",
        "",
        f"D3(BJ) parameters: PBE s6=1, s8={S8}, s9={S9}, a1={A1}, a2={A2}; these match the pairwise VASP IVDW=12 PBE-D3(BJ) convention.",
        "The current-DPA estimate is E_DPA + analytic D3(BJ). DFT-relaxed and archived snapshot energies are already PBE-D3(BJ), so their D3 values are reported as corrections only and are not added again.",
        "",
        "| Pair | ΔE(DPA+D3) (eV/atom) | ΔE(DFT final) (eV/atom) | model−DFT (meV/atom) | ΔE(snapshot) (eV/atom) | model−snapshot (meV/atom) |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in summary_rows:
        lines.append("| {pair} | {delta_E_DPA_plus_D3_eV_per_atom:.8f} | {delta_E_DFT_final_eV_per_atom:.8f} | {model_minus_DFT_meV_per_atom:.3f} | {delta_E_DFT_snapshot_eV_per_atom:.8f} | {model_minus_snapshot_meV_per_atom:.3f} |".format(**row))
    lines += ["", "Expanse job 54600945 (PBE, IVDW=0) was not accessible from this environment; no pure-PBE VASP output was found locally, so the direct PBE versus PBE-D3(BJ) comparison remains pending."]
    (output / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    metadata = {"status": "complete_partial_expanse_pending", "functional": "PBE-D3(BJ)", "d3_parameters": {"s6": 1.0, "s8": S8, "s9": S9, "a1": A1, "a2": A2}, "git_commit": subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip(), "generated_at_utc": datetime.now(timezone.utc).isoformat(), "inputs": [{"path": str(path.relative_to(root)) if path.is_relative_to(root) else str(path), "sha256": sha256(path)} for path in input_paths], "free_energy_inputs": [], "expanse_job": {"job_id": "54600945", "status": "not_available_from_current_environment"}}
    (output / "findings.meta.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
