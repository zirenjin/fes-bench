"""Safely extract configured ZIP or tar archives into the benchmark data area."""

from __future__ import annotations

import argparse
import tarfile
import zipfile
from pathlib import Path

from fes_bench.config import ConfigError, load_mapping


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="YAML/JSON extraction configuration")
    return parser


def _resolve(config_path: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else (config_path.parent / path).resolve()


def _safe_target(destination: Path, member_name: str) -> Path:
    candidate = (destination / member_name).resolve()
    if destination not in candidate.parents and candidate != destination:
        raise ConfigError(f"archive member escapes destination: {member_name!r}")
    return candidate


def _extract_zip(archive: Path, destination: Path) -> int:
    with zipfile.ZipFile(archive) as source:
        members = source.infolist()
        for member in members:
            _safe_target(destination, member.filename)
        source.extractall(destination)
        return len(members)


def _extract_tar(archive: Path, destination: Path) -> int:
    with tarfile.open(archive) as source:
        members = source.getmembers()
        for member in members:
            if member.issym() or member.islnk():
                raise ConfigError(f"links are disallowed in raw archives: {member.name!r}")
            _safe_target(destination, member.name)
        source.extractall(destination, members=members, filter="data")
        return len(members)


def run(config_path: str | Path) -> int:
    config_file = Path(config_path).expanduser().resolve()
    config = load_mapping(config_file)
    archives = config.get("archives")
    if not isinstance(archives, list) or not archives:
        raise ConfigError("archives must be a non-empty list")
    for item in archives:
        if not isinstance(item, dict) or not isinstance(item.get("archive"), str) or not isinstance(item.get("destination"), str):
            raise ConfigError("each archive needs string archive and destination fields")
        archive = _resolve(config_file, item["archive"])
        destination = _resolve(config_file, item["destination"])
        if not archive.is_file():
            raise ConfigError(f"archive does not exist: {archive}")
        destination.mkdir(parents=True, exist_ok=True)
        suffixes = archive.suffixes
        if archive.suffix == ".zip":
            count = _extract_zip(archive, destination)
        elif suffixes[-2:] == [".tar", ".gz"] or archive.suffix in {".tgz", ".tar"}:
            count = _extract_tar(archive, destination)
        else:
            raise ConfigError(f"unsupported archive format: {archive}")
        print(f"extracted: {archive} -> {destination} members={count}")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        return run(args.config)
    except (ConfigError, OSError, tarfile.TarError, zipfile.BadZipFile) as exc:
        print(f"extract failed: {exc}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
