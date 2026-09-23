"""Verify frozen split content hashes and train/test disjointness."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from fes_bench.config import load_mapping


def _payload_hash(payload: dict[str, Any]) -> str:
    unsigned = {key: value for key, value in payload.items() if key != "sha256"}
    canonical = json.dumps(unsigned, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def _triples(rows: object, label: str) -> set[tuple[str, str, int]]:
    if not isinstance(rows, list):
        raise ValueError(f"{label} must be a list")
    triples: set[tuple[str, str, int]] = set()
    for row in rows:
        if not isinstance(row, dict):
            raise ValueError(f"{label} has a non-object entry")
        try:
            triple = (str(row["system"]), str(row["phase"]), int(row["T_index"]))
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"{label} has an invalid frame triple") from exc
        if triple in triples:
            raise ValueError(f"{label} repeats frame triple {triple}")
        triples.add(triple)
    return triples


def verify(splits_dir: Path) -> dict[str, object]:
    """Return an integrity report or raise on a modified frozen split."""

    reports: dict[str, dict[str, object]] = {}
    for path in sorted(splits_dir.glob("*.json")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(payload, dict) or "sha256" not in payload:
            raise ValueError(f"{path} is not a frozen split JSON")
        stored = str(payload["sha256"])
        computed = _payload_hash(payload)
        if stored != computed:
            raise ValueError(f"{path} SHA256 mismatch: stored {stored}, computed {computed}")
        if "folds" in payload:
            folds = payload["folds"]
            if not isinstance(folds, dict):
                raise ValueError(f"{path} folds must be an object")
            fold_reports = {}
            for name, fold in folds.items():
                if not isinstance(fold, dict):
                    raise ValueError(f"{path} fold {name} is not an object")
                train = _triples(fold.get("train"), f"{path}:{name}:train")
                test = _triples(fold.get("test"), f"{path}:{name}:test")
                if train & test:
                    raise ValueError(f"{path} fold {name} overlaps train and test")
                fold_reports[str(name)] = {"train_frames": len(train), "test_frames": len(test)}
            reports[path.stem] = {"sha256": stored, "folds": fold_reports, "git_commit": payload.get("git_commit")}
        else:
            train = _triples(payload.get("train"), f"{path}:train")
            test = _triples(payload.get("test"), f"{path}:test")
            if train & test:
                raise ValueError(f"{path} overlaps train and test")
            reports[path.stem] = {
                "sha256": stored,
                "train_frames": len(train),
                "test_frames": len(test),
                "git_commit": payload.get("git_commit"),
            }
    if not reports:
        raise ValueError(f"no frozen split JSON files in {splits_dir}")
    return {"status": "pass", "splits": reports}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    args = parser.parse_args(argv)
    config_path = Path(args.config).resolve()
    config = load_mapping(config_path)
    root = Path(str(config.get("splits_dir", "../splits")))
    splits_dir = root if root.is_absolute() else (config_path.parent / root).resolve()
    report = verify(splits_dir)
    output = Path(str(config.get("output", "../result/experiments/leakage/split_integrity.json")))
    output = output if output.is_absolute() else (config_path.parent / output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
