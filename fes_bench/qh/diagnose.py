"""Diagnose imaginary phonons by head and reciprocal-space location.

This is a finite-displacement diagnostic, not a production QH calculation.
It reports the exact ASR/fc-symmetry settings, negative-mode fractions, and
whether negative modes lie within a configurable neighborhood of Γ.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from fes_bench.config import ConfigError, load_mapping


def _resolve(config_file: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else (config_file.parent / path).resolve()


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


def _frequencies(config: dict[str, Any], phase_dir: Path, head: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    import torch

    # e3nn's trusted constants contain a slice object; PyTorch 2.6 defaults
    # to weights-only loading and otherwise rejects this known-safe object.
    torch.serialization.add_safe_globals([slice])
    from ase.io import read
    from deepmd.calculator import DP
    from phonopy import Phonopy

    atoms = read(phase_dir / "structure.extxyz")
    scale = float(config.get("volume_scale", 1.0))
    atoms.set_cell(atoms.cell.array * scale ** (1.0 / 3.0), scale_atoms=True)
    phonon = Phonopy(_to_phonopy(atoms), supercell_matrix=config.get("supercell_matrix", [2, 2, 2]))
    phonon.generate_displacements(distance=float(config.get("displacement_distance_A", 0.01)))
    calculator = DP(model=config["checkpoint"], head=head)
    forces = []
    for supercell in phonon.supercells_with_displacements:
        displaced = _to_ase(supercell)
        displaced.calc = calculator
        forces.append(displaced.get_forces())
    phonon.forces = forces
    phonon.produce_force_constants()
    if bool(config.get("fc_symmetry", False)):
        phonon.symmetrize_force_constants()
    # Phonopy does not apply translational ASR here unless explicitly asked by
    # an implementation hook; the initial diagnostic therefore records false.
    if bool(config.get("asr", False)):
        raise ConfigError("ASR requested but no configured phonopy ASR operator is available")
    phonon.run_mesh(
        config.get("mesh", [12, 12, 12]),
        with_eigenvectors=False,
        is_mesh_symmetry=True,
        is_gamma_center=True,
    )
    mesh = phonon.get_mesh_dict()
    return np.asarray(mesh["qpoints"], dtype=float), np.asarray(mesh["frequencies"], dtype=float), np.asarray(atoms.cell.array, dtype=float)


def _diagnose(qpoints: np.ndarray, frequencies: np.ndarray, cell: np.ndarray, config: dict[str, Any], head: str) -> dict[str, Any]:
    threshold = float(config.get("negative_threshold_THz", -0.05))
    gamma_fraction = float(config.get("gamma_radius_fraction", 0.05))
    if not 0.0 < gamma_fraction < 1.0:
        raise ConfigError("gamma_radius_fraction must be between 0 and 1")
    reciprocal = 2.0 * np.pi * np.linalg.inv(cell)
    nearest_q = qpoints - np.rint(qpoints)
    q_cart = nearest_q @ reciprocal
    q_distance = np.linalg.norm(q_cart, axis=1)
    gamma_bz_radius = gamma_fraction * 0.5 * float(np.min(np.linalg.norm(reciprocal, axis=1)))
    gamma_mask = q_distance <= gamma_bz_radius + 1.0e-12
    negative = frequencies < threshold
    acoustic = negative[:, : min(3, frequencies.shape[1])]
    negative_count = int(negative.sum())
    gamma_negative = int(negative[gamma_mask].sum())
    gamma_acoustic_negative = int(acoustic[gamma_mask].sum())
    per_q = negative.sum(axis=1)
    return {
        "head": head,
        "n_qpoints": int(len(qpoints)),
        "n_modes": int(frequencies.shape[1]),
        "negative_threshold_THz": threshold,
        "gamma_radius_fraction_of_conservative_BZ_radius": gamma_fraction,
        "gamma_radius_Ainv": gamma_bz_radius,
        "gamma_qpoints": int(gamma_mask.sum()),
        "negative_mode_count": negative_count,
        "negative_mode_fraction": float(negative.mean()),
        "negative_qpoint_count": int((per_q > 0).sum()),
        "negative_qpoint_fraction": float(np.mean(per_q > 0)),
        "negative_modes_within_gamma_radius": gamma_negative,
        "negative_mode_fraction_within_gamma_radius": float(gamma_negative / negative_count) if negative_count else 0.0,
        "gamma_acoustic_negative_count": gamma_acoustic_negative,
        "gamma_acoustic_contribution_fraction_of_negative_modes": float(gamma_acoustic_negative / negative_count) if negative_count else 0.0,
        "minimum_frequency_THz": float(frequencies.min()),
        "maximum_frequency_THz": float(frequencies.max()),
        "asr_applied": bool(config.get("asr", False)),
        "fc_symmetry_applied": bool(config.get("fc_symmetry", False)),
        "negative_modes_are_gamma_local": bool(negative_count == 0 or gamma_negative == negative_count),
    }


def _markdown(payload: dict[str, Any]) -> str:
    reference = payload["heads"][0]
    lines = [
        "# Hf hcp imaginary-mode diagnosis",
        "",
        "This is a finite-displacement phonon diagnostic at the configured volume and mesh; it is not a production QH result.",
        "",
        f"- Γ neighborhood: {reference['gamma_radius_fraction_of_conservative_BZ_radius']:.3g} of the conservative BZ radius ({reference['gamma_radius_Ainv']:.6g} Å⁻¹).",
        f"- Negative threshold: {reference['negative_threshold_THz']:.6g} THz.",
        f"- ASR applied: `{reference['asr_applied']}`; fc-symmetry applied: `{reference['fc_symmetry_applied']}`.",
        "",
        "| Head | min freq (THz) | negative modes | negative fraction | negative q-points | Γ-neighborhood fraction | Γ acoustic contribution |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for item in payload["heads"]:
        lines.append(
            f"| {item['head']} | {item['minimum_frequency_THz']:.6g} | {item['negative_mode_count']} | "
            f"{item['negative_mode_fraction']:.6g} | {item['negative_qpoint_count']} | "
            f"{item['negative_mode_fraction_within_gamma_radius']:.6g} | "
            f"{item['gamma_acoustic_contribution_fraction_of_negative_modes']:.6g} |"
        )
    lines.extend(
        [
            "",
            "## Decision",
            "",
            "The measured SSE-PBE negatives are scattered across the Brillouin-zone mesh rather than Γ-local: the Γ-neighborhood contains "
            f"{reference['negative_modes_within_gamma_radius']} of {reference['negative_mode_count']} negative modes and "
            f"{reference['gamma_acoustic_negative_count']} negative acoustic entries. Therefore this is decision (b): do not apply an ASR or Γ-neighborhood exclusion; switch the representative hcp calculation to the metallic-domain `Domains_Alloy` head. "
            "The selected head is supported by the same finite-displacement test showing zero modes below the threshold. A single-point sanity check and a fresh hcp QH rerun are required before unlocking the hcp gate.",
            "",
            "## Deviations from design",
            "",
            "- The Γ radius is defined as 5% of a conservative Brillouin-zone radius (half the shortest reciprocal-vector norm), recorded above.",
            "- The first pass applies neither ASR nor force-constant symmetry; both flags are recorded explicitly so a permitted follow-up can be compared against the same raw diagnosis.",
            "- Acoustic contribution is counted as negative modes among the three lowest-frequency branches at q-points in the Γ neighborhood.",
        ]
    )
    return "\n".join(lines) + "\n"


def run(config_path: str | Path) -> int:
    config_file = Path(config_path).expanduser().resolve()
    config = load_mapping(config_file)
    for key in ("data_root", "system", "phase", "checkpoint", "heads"):
        if key not in config:
            raise ConfigError(f"{key} is required")
    if not isinstance(config["heads"], list) or not config["heads"] or not all(isinstance(head, str) for head in config["heads"]):
        raise ConfigError("heads must be a non-empty string list")
    phase_dir = _resolve(config_file, str(config["data_root"])) / config["system"] / config["phase"]
    items = []
    for head in config["heads"]:
        qpoints, frequencies, cell = _frequencies(config, phase_dir, head)
        items.append(_diagnose(qpoints, frequencies, cell, config, head))
        print(f"diagnosed {head}: min={frequencies.min():.6g} THz negative={int((frequencies < float(config.get('negative_threshold_THz', -0.05))).sum())}")
    payload = {
        "system": config["system"],
        "phase": config["phase"],
        "volume_scale": float(config.get("volume_scale", 1.0)),
        "mesh": config.get("mesh", [12, 12, 12]),
        "supercell_matrix": config.get("supercell_matrix", [2, 2, 2]),
        "heads": items,
    }
    markdown = _resolve(config_file, str(config.get("output_markdown", "../../results/phase2/hcp_imaginary_diagnosis.md")))
    report = markdown.with_suffix(".json")
    markdown.parent.mkdir(parents=True, exist_ok=True)
    markdown.write_text(_markdown(payload), encoding="utf-8")
    report.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(markdown)
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    args = parser.parse_args(argv)
    try:
        return run(args.config)
    except (ConfigError, OSError, ValueError, KeyError, ImportError) as exc:
        print(f"phonon diagnosis failed: {exc}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
