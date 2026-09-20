"""Export a spglib primitive conventional source as a benchmark structure."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from fes_bench.config import ConfigError, load_mapping


def _resolve(config_path: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else (config_path.parent / path).resolve()


def run(config_path: str | Path) -> int:
    try:
        import spglib
        from ase import Atoms
        from ase.io import read, write
    except ImportError as exc:
        raise ConfigError("ASE and spglib are required") from exc
    config_file = Path(config_path).expanduser().resolve()
    config = load_mapping(config_file)
    if not all(isinstance(config.get(key), str) for key in ("data_root", "system")):
        raise ConfigError("data_root and system are required strings")
    phases = config.get("phases")
    if not isinstance(phases, list):
        raise ConfigError("phases must be a list")
    system_dir = _resolve(config_file, config["data_root"]) / config["system"]
    for spec in phases:
        if not isinstance(spec, dict) or not isinstance(spec.get("name"), str) or not isinstance(spec.get("source"), str):
            raise ConfigError("each phase needs name and source strings")
        phase_dir = system_dir / spec["name"]
        source = _resolve(config_file, spec["source"])
        atoms = read(source)
        primitive = spglib.standardize_cell(
            (atoms.cell.array, atoms.get_scaled_positions(), atoms.numbers),
            to_primitive=True,
            no_idealize=False,
            symprec=float(spec.get("symprec_A", 0.01)),
        )
        if primitive is None:
            raise ConfigError(f"spglib could not standardize {source}")
        cell, scaled_positions, numbers = primitive
        result = Atoms(numbers=numbers, cell=cell, scaled_positions=scaled_positions, pbc=True)
        phase_dir.mkdir(parents=True, exist_ok=True)
        write(phase_dir / "structure.extxyz", result, format="extxyz")
        meta_path = phase_dir / "meta.json"
        meta = json.loads(meta_path.read_text(encoding="utf-8")) if meta_path.exists() else {}
        meta["representative"] = {
            "status": "source_standardized",
            "source": str(source),
            "method": "spglib.standardize_cell(to_primitive=True, no_idealize=False)",
            "source_atoms": len(atoms),
            "primitive_atoms": len(result),
            "fixed_for_all_temperatures": True,
        }
        meta_path.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
        print(f"standardized: {spec['name']} atoms={len(atoms)}->{len(result)}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    args = parser.parse_args(argv)
    try:
        return run(args.config)
    except (ConfigError, OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"structure standardization failed: {exc}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
