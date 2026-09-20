"""Small dependency-free configuration loader used by command modules."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


class ConfigError(ValueError):
    """Raised when a command configuration cannot be decoded safely."""


def load_mapping(path: str | Path) -> dict[str, Any]:
    """Load a mapping from JSON or JSON-compatible YAML.

    Phase 0 configurations are JSON documents with a ``.yaml`` extension;
    JSON is valid YAML and needs no optional parser. For conventional YAML,
    install PyYAML in the execution environment.
    """

    config_path = Path(path).expanduser().resolve()
    try:
        raw = config_path.read_text(encoding="utf-8")
    except OSError as exc:
        raise ConfigError(f"cannot read config {config_path}: {exc}") from exc

    try:
        value = json.loads(raw)
    except json.JSONDecodeError:
        try:
            import yaml  # type: ignore[import-not-found]
        except ImportError as exc:
            raise ConfigError(
                f"{config_path} is not JSON-compatible YAML; install PyYAML "
                "to use conventional YAML syntax"
            ) from exc
        value = yaml.safe_load(raw)

    if not isinstance(value, dict):
        raise ConfigError(f"config {config_path} must contain a mapping")
    return value
