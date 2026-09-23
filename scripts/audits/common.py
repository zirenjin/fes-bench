"""Shared output helpers for audit scripts."""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def pair_records(metrics: dict[str, Any]) -> list[tuple[str, dict[str, Any]]]:
    if isinstance(metrics.get("pairs"), dict):
        return [(str(key), value) for key, value in metrics["pairs"].items() if isinstance(value, dict)]
    result: list[tuple[str, dict[str, Any]]] = []
    for fold_name, fold in metrics.get("folds", {}).items() if isinstance(metrics.get("folds"), dict) else []:
        if not isinstance(fold, dict):
            continue
        for pair, record in pair_records(fold):
            item = dict(record); item["fold"] = fold_name; result.append((pair, item))
    return result


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def commit(root: Path) -> str:
    try:
        return subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def write_audit(root: Path, name: str, fields: list[str], rows: list[dict[str, Any]], inputs: list[Path], missing: list[dict[str, str]], extra: dict[str, Any] | None = None) -> None:
    out = root / "results" / "audits" / name
    out.mkdir(parents=True, exist_ok=True)
    with (out / "findings.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader(); writer.writerows({key: ("" if value is None else value) for key, value in row.items()} for row in rows)
    meta = {"generated_utc": datetime.now(timezone.utc).isoformat(), "git_commit": commit(root), "inputs": [{"path": str(path.relative_to(root)), "sha256": sha256(path)} for path in sorted(set(inputs)) if path.exists()], "missing": missing}
    if extra:
        meta.update(extra)
    (out / "findings.meta.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
