"""Run a read-only single-point force/energy sanity check for QH heads."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from fes_bench.config import ConfigError, load_mapping


def _resolve(config_file: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else (config_file.parent / path).resolve()


def run(config_path: str | Path) -> int:
    config_file = Path(config_path).expanduser().resolve()
    config = load_mapping(config_file)
    for key in ("structure", "checkpoint", "heads", "output"):
        if key not in config:
            raise ConfigError(f"{key} is required")
    if not isinstance(config["heads"], list) or not config["heads"]:
        raise ConfigError("heads must be a non-empty list")
    try:
        import torch
        torch.serialization.add_safe_globals([slice])
        import numpy as np
        from ase.io import read
        from deepmd.calculator import DP
    except ImportError as exc:
        raise ConfigError("ASE, NumPy, and DeepMD are required") from exc
    atoms = read(_resolve(config_file, str(config["structure"])))
    rows = []
    for head in config["heads"]:
        atoms.calc = DP(model=config["checkpoint"], head=head)
        energy = float(atoms.get_potential_energy() / len(atoms))
        forces = np.asarray(atoms.get_forces(), dtype=float)
        rows.append({
            "head": str(head),
            "n_atoms": int(len(atoms)),
            "energy_eV_per_atom": energy,
            "max_force_eV_per_A": float(np.abs(forces).max()),
            "rms_force_eV_per_A": float(np.sqrt(np.mean(forces ** 2))),
        })
        print(f"{head}: E={energy:.12g} eV/atom max|F|={rows[-1]['max_force_eV_per_A']:.12g} eV/A")
    output = _resolve(config_file, str(config["output"]))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps({"structure": str(config["structure"]), "checkpoint": str(config["checkpoint"]), "heads": rows}, indent=2) + "\n", encoding="utf-8")
    print(output)
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    args = parser.parse_args(argv)
    try:
        return run(args.config)
    except (ConfigError, OSError, ValueError, KeyError, ImportError) as exc:
        print(f"single-point sanity failed: {exc}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
