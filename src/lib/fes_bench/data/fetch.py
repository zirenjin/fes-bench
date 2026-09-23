"""Fetch and checksum raw benchmark archives from explicit configuration."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any
from urllib.request import urlopen

from fes_bench.config import ConfigError, load_mapping


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="YAML/JSON source configuration")
    return parser


def _resolve(config_path: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else (config_path.parent / path).resolve()


def _download(url: str, destination: Path) -> str:
    temporary = destination.with_suffix(destination.suffix + ".part")
    digest = hashlib.sha256()
    with urlopen(url, timeout=120) as response, temporary.open("wb") as output:
        while chunk := response.read(1024 * 1024):
            output.write(chunk)
            digest.update(chunk)
    temporary.replace(destination)
    return digest.hexdigest()


def run(config_path: str | Path) -> int:
    config_file = Path(config_path).expanduser().resolve()
    config = load_mapping(config_file)
    raw_root_value = config.get("raw_root")
    if not isinstance(raw_root_value, str):
        raise ConfigError("raw_root must be a path string")
    raw_root = _resolve(config_file, raw_root_value)
    entries = config.get("downloads")
    if not isinstance(entries, list) or not entries:
        raise ConfigError("downloads must be a non-empty list")

    manifests: dict[str, list[dict[str, Any]]] = {}
    for entry in entries:
        if not isinstance(entry, dict):
            raise ConfigError("each downloads entry must be a mapping")
        system = entry.get("system")
        url = entry.get("url")
        filename = entry.get("filename")
        doi = entry.get("doi")
        if not all(isinstance(value, str) and value for value in (system, url, filename, doi)):
            raise ConfigError("each download needs non-empty system, url, filename, and doi")
        destination = raw_root / system / filename
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.is_file():
            digest = hashlib.sha256(destination.read_bytes()).hexdigest()
            action = "reused"
        else:
            digest = _download(url, destination)
            action = "downloaded"
        manifests.setdefault(system, []).append(
            {
                "url": url,
                "doi": doi,
                "filename": str(destination.relative_to(raw_root.parent)),
                "sha256": digest,
                "downloaded_utc": datetime.now(UTC).isoformat(),
                "action": action,
            }
        )
        print(f"{action}: {system}/{filename} sha256={digest}")

    for system, records in manifests.items():
        manifest_path = raw_root.parent / system / "download.json"
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        manifest_path.write_text(
            json.dumps({"schema_version": 1, "files": records}, indent=2) + "\n",
            encoding="utf-8",
        )
        print(f"manifest: {manifest_path}")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        return run(args.config)
    except (ConfigError, OSError, ValueError) as exc:
        print(f"fetch failed: {exc}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
