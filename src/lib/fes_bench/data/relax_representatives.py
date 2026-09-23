"""Relax exported representative structures with DPA and record an audit trail."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from fes_bench.config import ConfigError, load_mapping


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="YAML/JSON relaxation configuration")
    return parser


def _resolve(config_path: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else (config_path.parent / path).resolve()


def _space_group(atoms) -> str | None:
    try:
        import spglib
        dataset = spglib.get_symmetry_dataset(
            (atoms.cell.array, atoms.get_scaled_positions(), atoms.numbers), symprec=0.01
        )
        if dataset is None:
            return None
        return str(getattr(dataset, "international", dataset["international"]))
    except Exception:
        return None


def run(config_path: str | Path) -> int:
    try:
        from ase.filters import UnitCellFilter
        from ase.io import read, write
        from ase.optimize import FIRE
        from deepmd.calculator import DP
    except ImportError as exc:
        raise ConfigError("ASE and deepmd are required for representative relaxation") from exc
    config_file = Path(config_path).expanduser().resolve()
    config = load_mapping(config_file)
    if not isinstance(config.get("data_root"), str) or not isinstance(config.get("system"), str):
        raise ConfigError("data_root and system are required strings")
    checkpoint = config.get("checkpoint")
    head = config.get("head")
    phases = config.get("phases")
    if not isinstance(checkpoint, str) or not isinstance(phases, list) or not all(isinstance(item, str) for item in phases):
        raise ConfigError("checkpoint and string phase list are required")
    if head is not None and not isinstance(head, str):
        raise ConfigError("head must be a string when supplied")
    fmax = float(config.get("fmax_eV_per_A", 0.01))
    max_steps = int(config.get("max_steps", 1000))
    cell_relaxed = bool(config.get("cell_relaxed", True))
    system_dir = _resolve(config_file, config["data_root"]) / config["system"]
    calculator = DP(model=checkpoint, head=head)
    for phase in phases:
        phase_dir = system_dir / phase
        meta_path = phase_dir / "meta.json"
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        source_path = phase_dir / "structure.source.extxyz"
        if not source_path.is_file():
            print(f"skipped: {phase}: missing {source_path.name}")
            continue
        atoms = read(source_path)
        atoms.calc = calculator
        before_energy = atoms.get_potential_energy() / len(atoms)
        before_volume = atoms.get_volume() / len(atoms)
        before_positions = atoms.get_positions().copy()
        before_group = _space_group(atoms)
        logfile = phase_dir / "representative_relaxation.log"
        optimizer = FIRE(UnitCellFilter(atoms) if cell_relaxed else atoms, logfile=str(logfile))
        optimizer.run(fmax=fmax, steps=max_steps)
        after_energy = atoms.get_potential_energy() / len(atoms)
        after_volume = atoms.get_volume() / len(atoms)
        after_group = _space_group(atoms)
        displacement = ((atoms.get_positions() - before_positions) ** 2).sum(axis=1) ** 0.5
        write(phase_dir / "structure.extxyz", atoms, format="extxyz")
        representative = dict(meta.get("representative", {}))
        representative.update(
            {
                "status": "relaxed",
                "calculator": "deepmd.calculator.DP",
                "checkpoint": checkpoint,
                "head": head,
                "optimizer": "FIRE",
                "cell_relaxed": cell_relaxed,
                "fmax_eV_per_A": fmax,
                "max_steps": max_steps,
                "converged": bool(optimizer.converged()),
                "energy_eV_per_atom": {"before": before_energy, "after": after_energy},
                "volume_A3_per_atom": {"before": before_volume, "after": after_volume},
                "space_group": {"before": before_group, "after": after_group},
                "max_displacement_A": float(displacement.max()),
                "relaxation_log": logfile.name,
            }
        )
        meta["representative"] = representative
        meta_path.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
        print(f"relaxed: {phase} converged={representative['converged']} max_disp_A={representative['max_displacement_A']:.6g}")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        return run(args.config)
    except (ConfigError, OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"representative relaxation failed: {exc}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
