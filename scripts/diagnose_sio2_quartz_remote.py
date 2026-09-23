"""Standalone quartz head diagnostic for thu-GenSi (no repository dependency)."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import torch
torch.serialization.add_safe_globals([slice])
from ase.io import read
from deepmd.calculator import DP
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms


ROOT = Path(".")
PHASE = ROOT / "data/sio2/quartz_beta"
OUT = ROOT / "runs/phase2_sio2_quartz_head_20260922"
CHECKPOINT = "external/checkpoints/DPA-3.1-3M.pt"
HEADS = ["Domains_Alloy", "Domains_SSE_PBE"]


def to_phonopy(atoms):
    return PhonopyAtoms(
        symbols=atoms.get_chemical_symbols(),
        cell=atoms.cell.array,
        scaled_positions=atoms.get_scaled_positions(),
    )


def to_ase(cell):
    from ase import Atoms

    return Atoms(symbols=cell.symbols, cell=cell.cell, scaled_positions=cell.scaled_positions, pbc=True)


def one(head: str) -> dict[str, object]:
    atoms = read(PHASE / "structure.extxyz")
    calculator = DP(model=CHECKPOINT, head=head)
    atoms.calc = calculator
    energy = float(atoms.get_potential_energy() / len(atoms))
    fmax = float(np.abs(atoms.get_forces()).max())
    phonon = Phonopy(to_phonopy(atoms), supercell_matrix=[2, 2, 2])
    phonon.generate_displacements(distance=0.01)
    forces = []
    for supercell in phonon.supercells_with_displacements:
        displaced = to_ase(supercell)
        displaced.calc = calculator
        forces.append(displaced.get_forces())
    phonon.forces = forces
    phonon.produce_force_constants()
    phonon.run_mesh([12, 12, 12], with_eigenvectors=False, is_mesh_symmetry=True, is_gamma_center=True)
    mesh = phonon.get_mesh_dict()
    qpoints = np.asarray(mesh["qpoints"], dtype=float)
    frequencies = np.asarray(mesh["frequencies"], dtype=float)
    reciprocal = 2.0 * np.pi * np.linalg.inv(atoms.cell.array)
    q_cart = (qpoints - np.rint(qpoints)) @ reciprocal
    q_distance = np.linalg.norm(q_cart, axis=1)
    gamma_radius = 0.05 * 0.5 * float(np.min(np.linalg.norm(reciprocal, axis=1)))
    gamma = q_distance <= gamma_radius + 1e-12
    negative = frequencies < -0.05
    per_q = negative.sum(axis=1)
    negative_count = int(negative.sum())
    gamma_negative = int(negative[gamma].sum())
    return {
        "head": head,
        "single_point_energy_eV_per_atom": energy,
        "single_point_fmax_eV_per_A": fmax,
        "n_qpoints": int(len(qpoints)),
        "n_modes": int(frequencies.shape[1]),
        "negative_threshold_THz": -0.05,
        "gamma_radius_fraction": 0.05,
        "gamma_radius_Ainv": gamma_radius,
        "gamma_qpoints": int(gamma.sum()),
        "minimum_frequency_THz": float(frequencies.min()),
        "negative_mode_count": negative_count,
        "negative_mode_fraction": float(negative.mean()),
        "negative_qpoint_count": int((per_q > 0).sum()),
        "negative_qpoint_fraction": float(np.mean(per_q > 0)),
        "negative_modes_within_gamma": gamma_negative,
        "negative_mode_fraction_within_gamma": float(gamma_negative / negative_count) if negative_count else 0.0,
    }


OUT.mkdir(parents=True, exist_ok=True)
payload = {"system": "sio2", "phase": "quartz_beta", "heads": [one(head) for head in HEADS]}
(OUT / "head_diagnosis.json").write_text(json.dumps(payload, indent=2) + "\n")
lines = ["# SiO2 β-quartz head diagnosis", "", "- threshold: `-0.05 THz`", "- Γ neighborhood: `5%` of conservative BZ radius", "- ASR: `false`; force-constant symmetry: `false`", "", "| head | E (eV/atom) | fmax (eV/A) | min freq (THz) | negative modes | negative q-points | Γ negative / all negative |", "|---|---:|---:|---:|---:|---:|---:|"]
for item in payload["heads"]:
    lines.append(f"| {item['head']} | {item['single_point_energy_eV_per_atom']:.8f} | {item['single_point_fmax_eV_per_A']:.6g} | {item['minimum_frequency_THz']:.6g} | {item['negative_mode_count']}/{item['n_qpoints']*item['n_modes']} ({item['negative_mode_fraction']:.4g}) | {item['negative_qpoint_count']}/{item['n_qpoints']} ({item['negative_qpoint_fraction']:.4g}) | {item['negative_modes_within_gamma']}/{item['negative_mode_count'] if item['negative_mode_count'] else 0} |")
lines += ["", "Decision is made against the Hf rule: all-zone negative modes mean out-of-domain; zero or Γ-local negatives pass.", ""]
(OUT / "sio2_alloy_head_check.md").write_text("\n".join(lines))
print(OUT / "sio2_alloy_head_check.md")
