"""Configuration-driven audit for canonical benchmark datasets."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from fes_bench.config import ConfigError, load_mapping

from .load import default_data_root, load
from .schema import SchemaError


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="YAML/JSON configuration path")
    return parser


def run(config_path: str | Path) -> int:
    resolved_config = Path(config_path).expanduser().resolve()
    config = load_mapping(resolved_config)
    raw_root = config.get("data_root")
    if isinstance(raw_root, str):
        configured_root = Path(raw_root).expanduser()
        root = (
            configured_root.resolve()
            if configured_root.is_absolute()
            else (resolved_config.parent / configured_root).resolve()
        )
    else:
        root = default_data_root()
    requested = config.get("systems", [])
    if not isinstance(requested, list) or not all(isinstance(item, str) for item in requested):
        raise ConfigError("systems must be a list of system directory names")
    systems = requested or sorted(
        entry.name for entry in root.iterdir() if entry.is_dir() and (entry / "system.json").is_file()
    )
    print(f"data_root={root}")
    print(f"systems={len(systems)}")
    failed = False
    for system in systems:
        system_json = root / system / "system.json"
        try:
            declared = json.loads(system_json.read_text(encoding="utf-8"))["phases"]
            phases = [load(system, phase, root) for phase in declared]
        except (OSError, KeyError, TypeError, json.JSONDecodeError, SchemaError) as exc:
            print(f"{system}: INVALID: {exc}")
            failed = True
            continue
        temperatures = [point.T_K for phase in phases for point in phase.G_table]
        pressures = [point.P_GPa for phase in phases for point in phase.G_table]
        print(
            f"{system}: phases={len(phases)} T_K=[{min(temperatures):g}, {max(temperatures):g}] "
            f"P_GPa=[{min(pressures):g}, {max(pressures):g}]"
        )
        crossings_path = root / system / "reference_crossings.json"
        try:
            crossing_data = json.loads(crossings_path.read_text(encoding="utf-8"))
            pairs = crossing_data["pairs"]
            summary = ", ".join(
                f"{name}:{len(value.get('crossings', []))}"
                for name, value in pairs.items()
            )
            print(f"{system}: crossings={summary or 'none'}")
        except (OSError, KeyError, TypeError, json.JSONDecodeError) as exc:
            print(f"{system}: crossings=INVALID ({exc})")
            failed = True
        sources = ", ".join(
            f"{phase.name}:{phase.meta.get('representative', {}).get('source', phase.meta.get('representative', {}).get('checkpoint', 'not-recorded'))}"
            for phase in phases
        )
        print(f"{system}: representative_sources={sources}")
    return 2 if failed else 0


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        return run(args.config)
    except (ConfigError, OSError, SchemaError) as exc:
        print(f"audit failed: {exc}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
