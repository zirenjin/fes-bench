"""Compare fixed-volume harmonic F_vib before/after Hf PES fine-tuning."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


EV_PER_KJMOL = 1.0 / 96.4853321233


def _phonopy_atoms(atoms):
    from phonopy.structure.atoms import PhonopyAtoms

    return PhonopyAtoms(symbols=atoms.get_chemical_symbols(), cell=atoms.cell.array, scaled_positions=atoms.get_scaled_positions())


def phase_free_energy(model, atoms, temperatures):
    from phonopy import Phonopy

    unit = _phonopy_atoms(atoms)
    phonon = Phonopy(unit, supercell_matrix=[2, 2, 2])
    phonon.generate_displacements(distance=0.01)
    displaced = phonon.supercells_with_displacements
    coords = np.asarray([s.positions.reshape(-1) for s in displaced], dtype=float)
    cells = np.asarray([np.asarray(s.cell).reshape(-1) for s in displaced], dtype=float)
    types = np.zeros(len(displaced[0]), dtype=np.int32)
    force_parts = []
    for start in range(0, len(displaced), 16):
        _, f, _ = model.eval(coords[start:start + 16], cells[start:start + 16], types)
        force_parts.append(np.asarray(f, dtype=float).reshape(-1, len(displaced[0]), 3))
    phonon.forces = np.concatenate(force_parts, axis=0)
    phonon.produce_force_constants()
    phonon.run_mesh([6, 6, 6], with_eigenvectors=False, is_mesh_symmetry=True, is_gamma_center=True)
    phonon.run_thermal_properties(t_min=min(temperatures), t_max=max(temperatures), t_step=1.0)
    thermal = phonon.get_thermal_properties_dict()
    lookup = {float(t): float(f) * EV_PER_KJMOL / len(atoms) for t, f in zip(thermal["temperatures"], thermal["free_energy"])}
    return {
        "n_atoms": len(atoms),
        "n_displacements": len(displaced),
        "min_frequency_THz": float(np.asarray(phonon.get_mesh_dict()["frequencies"]).min()),
        "static_energy_eV_per_atom": float(np.asarray(model.eval(atoms.positions.reshape(1, -1), atoms.cell.array.reshape(1, -1), np.zeros(len(atoms), dtype=np.int32))[0]).reshape(-1)[0] / len(atoms)),
        "F_vib_eV_per_atom": {str(t): lookup[float(t)] for t in temperatures},
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hcp", required=True)
    parser.add_argument("--bcc", required=True)
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--finetuned", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    import torch

    torch.serialization.add_safe_globals([slice])
    from ase.io import read
    from deepmd.infer import DeepPot

    temperatures = [1076.0, 1500.0, 1919.0, 2327.0]
    atoms = {"hcp": read(args.hcp), "bcc": read(args.bcc)}
    models = {
        "frozen_pretrained": DeepPot(args.baseline, head="Domains_Alloy"),
        "pes_finetuned": DeepPot(args.finetuned, head=None),
    }
    payload = {"system": "hf", "protocol": "fixed volume harmonic; supercell=[2,2,2], mesh=[6,6,6], displacement=0.01 A", "temperatures_K": temperatures, "training_data_used": "Hf DFT energy/force OUTCAR only; no TI/free-energy data", "models": {}}
    for label, model in models.items():
        payload["models"][label] = {phase: phase_free_energy(model, structure, temperatures) for phase, structure in atoms.items()}
        h = payload["models"][label]["hcp"]
        b = payload["models"][label]["bcc"]
        payload["models"][label]["deltaG_hcp_minus_bcc_eV_per_atom"] = {
            str(t): (h["static_energy_eV_per_atom"] + h["F_vib_eV_per_atom"][str(t)]) - (b["static_energy_eV_per_atom"] + b["F_vib_eV_per_atom"][str(t)]) for t in temperatures
        }
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    lines = ["# Hf hcp/bcc PES fine-tune QH comparison", "", "Fixed-volume harmonic comparison; no TI/free-energy labels used.", "", "| model | ΔG(hcp−bcc) 1076 K (meV/atom) | 1500 K | 1919 K | 2327 K |", "|---|---:|---:|---:|---:|"]
    for label, item in payload["models"].items():
        vals = [1000 * item["deltaG_hcp_minus_bcc_eV_per_atom"][str(t)] for t in temperatures]
        lines.append(f"| {label} | " + " | ".join(f"{v:.6g}" for v in vals) + " |")
    lines += ["", "F_QH here is the fixed-volume harmonic F_vib component; it is not a volume-minimized production QH run.", "", "## Deviations from design", "", "* The Hf archive has no usable virial/stress labels (`ISIF=0`), so the PES fine-tune and this comparison use energy+force only (`pref_v=0`).", "* A compact 6³ mesh is used for the before/after magnitude check; the protocol is recorded above."]
    out.with_suffix(".md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
