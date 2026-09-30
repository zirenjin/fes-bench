"""Configuration-driven finite-displacement quasi-harmonic free energies.

The calculation is deliberately small and transparent: every displaced
supercell is evaluated with the configured DeepMD calculator, force constants
are fitted by phonopy, and the minimum of E(V)+F_vib(V,T) is reported on the
reference temperature grid.  It is intended for benchmark provenance, not as
an opaque phonopy-wrapper shortcut.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from fes_bench.config import ConfigError, load_mapping


EV_PER_KJMOL = 1.0 / 96.4853321233


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    return parser


def _resolve(config_path: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else (config_path.parent / path).resolve()


def _reference_temperatures(path: Path) -> list[float]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    temperatures = sorted({float(row["T_K"]) for row in rows})
    if not temperatures:
        raise ConfigError(f"{path} has no temperatures")
    return temperatures


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


def _thermal_for_volume(
    phonon, calculator, temperatures: list[float], mesh: list[int], negative_threshold: float
) -> tuple[dict[float, float], float, float, float]:
    displaced = phonon.supercells_with_displacements
    forces = []
    for supercell in displaced:
        atoms = _to_ase(supercell)
        atoms.calc = calculator
        forces.append(atoms.get_forces())
    phonon.forces = forces
    phonon.produce_force_constants()
    phonon.run_mesh(mesh, with_eigenvectors=False, is_mesh_symmetry=True)
    phonon.run_thermal_properties(
        t_min=min(temperatures), t_max=max(temperatures), t_step=1.0
    )
    properties = phonon.get_thermal_properties_dict()
    by_t = {
        float(temp): float(free_energy) * EV_PER_KJMOL / len(phonon.unitcell)
        for temp, free_energy in zip(properties["temperatures"], properties["free_energy"])
    }
    missing = [temperature for temperature in temperatures if temperature not in by_t]
    if missing:
        raise ConfigError(f"phonopy did not return requested temperatures: {missing[:5]}")
    phonon.run_mesh(mesh, with_eigenvectors=False, is_mesh_symmetry=True)
    frequencies = phonon.get_mesh_dict()["frequencies"]
    negative = frequencies < negative_threshold
    per_q = negative.sum(axis=1)
    return (
        {temperature: by_t[temperature] for temperature in temperatures},
        float(frequencies.min()),
        float(negative.mean()),
        float((per_q > 0).mean()),
    )


def _write_phase(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["T_K", "F_QH_eV_per_atom", "E_static_eV_per_atom", "volume_scale", "min_frequency_THz", "checkpoint_sha256", "head"],
        )
        writer.writeheader()
        writer.writerows(rows)


def run(config_path: str | Path) -> int:
    try:
        import numpy as np
        from ase.io import read
        from deepmd.calculator import DP
        import torch
        torch.serialization.add_safe_globals([slice])
        from phonopy import Phonopy
    except ImportError as exc:
        raise ConfigError("ASE, phonopy, numpy, and DeepMD are required for QH") from exc
    config_file = Path(config_path).expanduser().resolve()
    config = load_mapping(config_file)
    required = ("data_root", "system", "checkpoint", "phases")
    if any(not isinstance(config.get(key), str) for key in required[:3]) or not isinstance(config.get("phases"), list):
        raise ConfigError("data_root, system, checkpoint, and phases are required")
    phases = config["phases"]
    if not phases or not all(isinstance(phase, str) for phase in phases):
        raise ConfigError("phases must be a non-empty list of strings")
    supercell_matrix = config.get("supercell_matrix", [2, 2, 2])
    supercell_by_phase = config.get("supercell_matrix_by_phase", {})
    mesh = config.get("mesh", [16, 16, 16])
    volume_scales = config.get("volume_scales", [0.98, 1.0, 1.02])
    if not all(isinstance(value, (int, float)) and value > 0 for value in volume_scales):
        raise ConfigError("volume_scales must contain positive numeric values")
    if len(supercell_matrix) != 3 or len(mesh) != 3:
        raise ConfigError("supercell_matrix and mesh must each have length 3")
    root = _resolve(config_file, config["data_root"]) / config["system"]
    output_root = _resolve(config_file, str(config.get("output_root", "../../result/experiments/quasi_harmonic/raw_runs"))) / config["system"]
    output_root.mkdir(parents=True, exist_ok=True)
    calculator = DP(model=config["checkpoint"], head=config.get("head"))
    summary: dict[str, object] = {
        "system": config["system"],
        "checkpoint": config["checkpoint"],
        "head": config.get("head"),
        "checkpoint_sha256": config.get("checkpoint_sha256"),
        "supercell_matrix": supercell_matrix,
        "mesh": mesh,
        "volume_scales": volume_scales,
        "phases": {},
    }
    for phase in phases:
        phase_dir = root / phase
        atoms0 = read(phase_dir / "structure.extxyz")
        temperatures = _reference_temperatures(phase_dir / "reference_G.csv")
        candidates: list[tuple[float, float, dict[float, float], float, float, float]] = []
        for scale in volume_scales:
            atoms = atoms0.copy()
            atoms.set_cell(atoms0.cell.array * float(scale) ** (1.0 / 3.0), scale_atoms=True)
            atoms.calc = calculator
            static_energy = atoms.get_potential_energy() / len(atoms)
            phase_supercell = supercell_by_phase.get(phase, supercell_matrix)
            phonon = Phonopy(_to_phonopy(atoms), supercell_matrix=phase_supercell)
            phonon.generate_displacements(distance=float(config.get("displacement_distance_A", 0.01)))
            thermal, min_frequency, negative_fraction, negative_qpoint_fraction = _thermal_for_volume(
                phonon,
                calculator,
                temperatures,
                mesh,
                float(config.get("imaginary_frequency_cutoff_THz", -0.05)),
            )
            candidates.append((float(scale), float(static_energy), thermal, min_frequency, negative_fraction, negative_qpoint_fraction))
        rows: list[dict[str, object]] = []
        for temperature in temperatures:
            scale, energy, thermal, min_frequency, _, _ = min(candidates, key=lambda item: item[1] + item[2][temperature])
            rows.append(
                {
                    "T_K": f"{temperature:.12g}",
                    "F_QH_eV_per_atom": f"{energy + thermal[temperature]:.16g}",
                    "E_static_eV_per_atom": f"{energy:.16g}",
                    "volume_scale": f"{scale:.12g}",
                    "min_frequency_THz": f"{min_frequency:.12g}",
                    "checkpoint_sha256": config.get("checkpoint_sha256", ""),
                    "head": config.get("head", ""),
                }
            )
        _write_phase(output_root / f"{phase}_fqh.csv", rows)
        summary["phases"][phase] = {
            "supercell_matrix": supercell_by_phase.get(phase, supercell_matrix),
            "n_temperatures": len(rows),
            "minimum_frequency_THz_by_volume": {str(scale): min_frequency for scale, _, _, min_frequency, _, _ in candidates},
            "negative_mode_fraction_by_volume": {str(scale): negative_fraction for scale, _, _, _, negative_fraction, _ in candidates},
            "negative_qpoint_fraction_by_volume": {str(scale): negative_qpoint_fraction for scale, _, _, _, _, negative_qpoint_fraction in candidates},
            "output": f"{phase}_fqh.csv",
        }
        print(f"QH complete: {phase} temperatures={len(rows)}")
    (output_root / "qh_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        return run(args.config)
    except (ConfigError, OSError, ValueError, KeyError) as exc:
        print(f"QH calculation failed: {exc}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
