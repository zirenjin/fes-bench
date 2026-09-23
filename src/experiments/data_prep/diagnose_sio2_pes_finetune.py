"""Compare phonons from frozen DPA and the existing SiO2 PES-fine-tuned model.

This is a deliberately isolated, fixed-cell diagnostic.  It evaluates only
energy/forces from the two PES checkpoints and never reads FES/TI data.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def _to_phonopy(atoms):
    from phonopy.structure.atoms import PhonopyAtoms

    return PhonopyAtoms(
        symbols=atoms.get_chemical_symbols(),
        cell=atoms.cell.array,
        scaled_positions=atoms.get_scaled_positions(),
    )


def _to_ase(cell):
    from ase import Atoms

    return Atoms(
        symbols=cell.symbols,
        cell=cell.cell,
        scaled_positions=cell.scaled_positions,
        pbc=True,
    )


def evaluate(label: str, checkpoint: str, structure: Path, supercell: list[int], mesh: list[int], head: str | None = None) -> dict[str, object]:
    import torch

    torch.serialization.add_safe_globals([slice])
    from ase.io import read
    from deepmd.infer import DeepPot
    from phonopy import Phonopy

    atoms = read(structure)
    model = DeepPot(checkpoint, head=head)
    atom_types = np.asarray([0 if symbol == "Si" else 1 for symbol in atoms.get_chemical_symbols()], dtype=np.int32)
    base_energy, base_forces, _ = model.eval(
        atoms.positions.reshape(1, -1), atoms.cell.array.reshape(1, -1), atom_types
    )
    phonon = Phonopy(_to_phonopy(atoms), supercell_matrix=supercell)
    phonon.generate_displacements(distance=0.01)
    displaced = phonon.supercells_with_displacements
    coords = np.asarray([cell.positions.reshape(-1) for cell in displaced], dtype=np.float64)
    cells = np.asarray([np.asarray(cell.cell).reshape(-1) for cell in displaced], dtype=np.float64)
    repeats = len(displaced[0]) // len(atom_types)
    supercell_types = np.tile(atom_types, repeats)
    types = supercell_types
    force_parts = []
    for start in range(0, len(displaced), 16):
        _, chunk_force, _ = model.eval(coords[start : start + 16], cells[start : start + 16], types)
        force_parts.append(np.asarray(chunk_force, dtype=np.float64).reshape(-1, len(atoms), 3))
    forces = np.concatenate(force_parts, axis=0)
    phonon.forces = forces
    phonon.produce_force_constants()
    phonon.run_mesh(mesh, with_eigenvectors=False, is_mesh_symmetry=True, is_gamma_center=True)
    frequencies = np.asarray(phonon.get_mesh_dict()["frequencies"], dtype=float)
    negative = frequencies < -0.05
    per_q = negative.sum(axis=1)
    return {
        "label": label,
        "checkpoint": checkpoint,
        "head": head,
        "single_point_energy_eV_per_atom": float(np.asarray(base_energy).reshape(-1)[0] / len(atoms)),
        "single_point_fmax_eV_per_A": float(np.abs(np.asarray(base_forces)).max()),
        "n_qpoints": int(frequencies.shape[0]),
        "n_modes": int(frequencies.shape[1]),
        "negative_threshold_THz": -0.05,
        "minimum_frequency_THz": float(frequencies.min()),
        "negative_mode_count": int(negative.sum()),
        "negative_mode_fraction": float(negative.mean()),
        "negative_qpoint_count": int((per_q > 0).sum()),
        "negative_qpoint_fraction": float(np.mean(per_q > 0)),
        "n_displacements": len(forces),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--structure", required=True)
    parser.add_argument("--baseline-checkpoint", required=True)
    parser.add_argument("--finetuned-checkpoint", required=True)
    parser.add_argument("--baseline-head", default="Domains_Alloy")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    structure = Path(args.structure)
    items = [
        evaluate("frozen_pretrained", args.baseline_checkpoint, structure, [1, 1, 1], [12, 12, 12], args.baseline_head),
        evaluate("pes_finetuned", args.finetuned_checkpoint, structure, [1, 1, 1], [12, 12, 12], None),
    ]
    payload = {
        "system": "sio2",
        "phase": "quartz_beta",
        "protocol": "fixed-cell finite displacement; supercell=[1,1,1] (unit-cell force constants); mesh=[12,12,12]; displacement=0.01 A",
        "training_data_used": "DFT PES energy/force/virial dataset only for the fine-tuned checkpoint provenance; no FES/TI data",
        "items": items,
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    lines = [
        "# SiO2 beta-quartz PES fine-tune phonon comparison",
        "",
        "Fixed-cell finite displacement; no FES/TI data used.",
        "",
        "| model | min freq (THz) | negative modes | negative q-points | fmax (eV/A) |",
        "|---|---:|---:|---:|---:|",
    ]
    for item in items:
        lines.append(
            f"| {item['label']} | {item['minimum_frequency_THz']:.8g} | "
            f"{item['negative_mode_count']}/{item['n_qpoints'] * item['n_modes']} "
            f"({item['negative_mode_fraction']:.6g}) | {item['negative_qpoint_count']}/{item['n_qpoints']} "
            f"({item['negative_qpoint_fraction']:.6g}) | {item['single_point_fmax_eV_per_A']:.6g} |"
        )
    lines += [
        "",
        "Interpretation: a large reduction/disappearance of all-zone negatives after PES fine-tuning would identify a model/PES issue; persistent all-zone negatives support a physical/structure instability.",
        "",
        "## Deviations from design",
        "",
        "* This is the required sanity comparison, not the full production QH grid; the cell is fixed and the unit-cell force-constant protocol is recorded above.",
        "* The existing PES checkpoint has `pref_v=0`; no virial label is silently substituted.",
    ]
    output.with_suffix(".md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
