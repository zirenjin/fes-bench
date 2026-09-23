"""Export provenance-preserving representative structures from source trajectories."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from fes_bench.config import ConfigError, load_mapping


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="YAML/JSON representative-source configuration")
    return parser


def _resolve(config_path: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else (config_path.parent / path).resolve()


def run(config_path: str | Path) -> int:
    try:
        from ase.io import read, write
    except ImportError as exc:
        raise ConfigError("ASE is required to export representative structures") from exc
    config_file = Path(config_path).expanduser().resolve()
    config = load_mapping(config_file)
    if not isinstance(config.get("data_root"), str) or not isinstance(config.get("system"), str):
        raise ConfigError("data_root and system are required strings")
    sources = config.get("sources")
    if not isinstance(sources, list) or not sources:
        raise ConfigError("sources must be a non-empty list")
    system_dir = _resolve(config_file, config["data_root"]) / config["system"]
    for source_spec in sources:
        if not isinstance(source_spec, dict) or not isinstance(source_spec.get("phase"), str):
            raise ConfigError("each source needs a phase")
        phase = source_spec["phase"]
        phase_dir = system_dir / phase
        meta_path = phase_dir / "meta.json"
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        source = source_spec.get("source")
        representative = dict(meta.get("representative", {}))
        if source is None:
            representative.update({"status": "source_missing", "source": None, "reason": source_spec.get("reason", "not supplied")})
            meta["representative"] = representative
            meta_path.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
            print(f"missing: {phase}: {representative['reason']}")
            continue
        if not isinstance(source, str):
            raise ConfigError(f"source for {phase} must be a string or null")
        source_path = _resolve(config_file, source)
        if not source_path.is_file():
            raise ConfigError(f"source for {phase} is absent: {source_path}")
        atoms = read(source_path, index=-1)
        destination = phase_dir / "structure.source.extxyz"
        write(destination, atoms, format="extxyz")
        representative.update(
            {
                "status": "source_exported_pending_relaxation",
                "source": str(source_path),
                "source_format": "VASP OUTCAR",
                "source_selection": "last ionic image in the configured source OUTCAR",
                "source_temperature_K": source_spec.get("source_temperature_K"),
                "atom_count": len(atoms),
                "unrelaxed_structure": destination.name,
            }
        )
        meta["representative"] = representative
        meta_path.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
        print(f"exported: {phase} atoms={len(atoms)} source={source_path}")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        return run(args.config)
    except (ConfigError, OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"representative export failed: {exc}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
