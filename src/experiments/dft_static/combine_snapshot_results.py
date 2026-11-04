"""Combine per-system snapshot-bias tables without recomputing completed systems."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


PHASES = {
    "hf": ("hcp", "bcc"),
    "ti": ("hcp", "bcc"),
    "zr": ("hcp", "bcc"),
    "sio2": ("quartz_beta", "cristobalite_beta", "tridymite_p63mmc"),
}
PAIRS = {
    "hf": (("hcp", "bcc"),),
    "ti": (("hcp", "bcc"),),
    "zr": (("hcp", "bcc"),),
    "sio2": (("quartz_beta", "cristobalite_beta"), ("quartz_beta", "tridymite_p63mmc")),
}
FIELDS = [
    "row_type", "system", "phase", "phase_low", "phase_high",
    "bias_meV_per_atom", "std_meV_per_atom", "n_snapshots", "status", "head", "source",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def canonical_rows(source_by_system: dict[str, Path]) -> list[dict[str, str]]:
    rows_by_system: dict[str, list[dict[str, str]]] = {}
    for system, path in source_by_system.items():
        rows_by_system[system] = read_rows(path)
    rows: list[dict[str, str]] = []
    for system, phases in PHASES.items():
        phase_rows = {
            row["phase"]: row
            for row in rows_by_system[system]
            if row["row_type"] == "phase"
        }
        pair_rows = {
            (row["phase_low"], row["phase_high"]): row
            for row in rows_by_system[system]
            if row["row_type"] == "pair"
        }
        rows.extend(phase_rows[phase] for phase in phases)
        rows.extend(pair_rows[pair] for pair in PAIRS[system])
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--ti", type=Path, required=True)
    parser.add_argument("--zr", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    args = parser.parse_args()

    source_by_system = {"hf": args.base, "sio2": args.base, "ti": args.ti, "zr": args.zr}
    rows = canonical_rows(source_by_system)
    bad = [row for row in rows if row["status"] != "ok"]
    if bad:
        raise SystemExit("snapshot table still contains unavailable rows")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    base_meta_path = args.base.with_name("dpa_vs_dft_snapshots.meta.json")
    base_meta = json.loads(base_meta_path.read_text(encoding="utf-8")) if base_meta_path.is_file() else {}
    metadata = {
        "status": "complete",
        "checkpoint": str(args.checkpoint),
        "checkpoint_sha256": sha256(args.checkpoint),
        "head_policy": "configs/models/head_policy.yaml",
        "systems": list(PHASES),
        "inference_mode": "Hf and SiO2 completed in the base table; Ti and Zr completed with isolated DeepMD workers.",
        "source_csv_sha256": {system: sha256(path) for system, path in source_by_system.items()},
        "base_metadata": base_meta,
    }
    (args.output.parent / "dpa_vs_dft_snapshots.meta.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({"output": str(args.output), "rows": len(rows), "status": "complete"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
