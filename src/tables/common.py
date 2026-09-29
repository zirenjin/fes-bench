"""Shared readers and provenance writers for table scripts."""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


def root_path(value: str | Path = ".") -> Path:
    return Path(value).resolve()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def git_commit(root: Path) -> str:
    try:
        return subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def raw_runs(root: Path, split: str | None = None) -> Iterable[tuple[Path, dict[str, Any]]]:
    pattern = f"result/experiments/external_baselines/raw_{split or '*'}/*/*/metrics.json"
    for path in sorted(root.glob(pattern)):
        yield path, json.loads(path.read_text(encoding="utf-8"))


def inventory(root: Path) -> dict[str, Any]:
    path = root / "result/experiments/data_prep/raw_inventory/inventory/seed_none/metrics.json"
    return json.loads(path.read_text(encoding="utf-8"))["metrics"]


def csv_write(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows({key: ("" if value is None else value) for key, value in row.items()} for row in rows)


def meta_write(root: Path, csv_path: Path, inputs: list[Path], missing: list[dict[str, str]], extra: dict[str, Any] | None = None) -> None:
    payload = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": git_commit(root),
        "inputs": [{"path": str(path.relative_to(root)), "sha256": sha256(path)} for path in sorted(set(inputs))],
        "missing": missing,
    }
    if extra:
        payload.update(extra)
    csv_path.with_suffix(".meta.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def numeric(value: Any) -> float | None:
    try:
        return None if value is None else float(value)
    except (TypeError, ValueError):
        return None


def pair_records(metrics: dict[str, Any]) -> list[tuple[str, dict[str, Any]]]:
    if isinstance(metrics.get("pairs"), dict):
        return [(str(key), value) for key, value in metrics["pairs"].items() if isinstance(value, dict)]
    records: list[tuple[str, dict[str, Any]]] = []
    for fold_name, fold in metrics.get("folds", {}).items() if isinstance(metrics.get("folds"), dict) else []:
        if not isinstance(fold, dict):
            continue
        for pair, record in pair_records(fold):
            record = dict(record)
            record["fold"] = fold_name
            records.append((pair, record))
    return records


def mean(values: list[float]) -> float | None:
    return sum(values) / len(values) if values else None


def aggregate(metrics: dict[str, Any]) -> dict[str, Any]:
    pairs = [record for _, record in pair_records(metrics)]
    return {
        "G_MAE_eV_per_atom": numeric(metrics.get("G_MAE_eV_per_atom")),
        "delta_G_MAE_eV_per_atom": mean([float(record["delta_G_MAE_eV_per_atom"]) for record in pairs if numeric(record.get("delta_G_MAE_eV_per_atom")) is not None]),
        "delta_G_RMSE_eV_per_atom": mean([float(record["delta_G_RMSE_eV_per_atom"]) for record in pairs if numeric(record.get("delta_G_RMSE_eV_per_atom")) is not None]),
        "sign_accuracy": mean([float(record["sign_accuracy"]) for record in pairs if numeric(record.get("sign_accuracy")) is not None]),
        "Tc_error_K": mean([abs(float(error)) for record in pairs for error in record.get("Tc_error_K", []) if numeric(error) is not None]),
        "false_crossings": sum(int(record.get("false_crossings", 0) or 0) for record in pairs),
        "missed_crossings": sum(int(record.get("missed_crossings", 0) or 0) for record in pairs),
    }
