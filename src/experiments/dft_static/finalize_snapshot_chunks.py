"""Convert per-phase resume-safe snapshot chunks into a system CSV."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


FIELDS = [
    "row_type", "system", "phase", "phase_low", "phase_high",
    "bias_meV_per_atom", "std_meV_per_atom", "n_snapshots", "status",
    "head", "source",
]


def read_phase(path: Path) -> dict[str, object]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if value.get("status") != "ok":
        raise SystemExit(f"phase result is not complete: {path}")
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--system", choices=("ti", "zr"), required=True)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    phases = {phase: read_phase(args.root / phase / "phase_result.json") for phase in ("hcp", "bcc")}
    rows: list[dict[str, object]] = []
    for phase in ("hcp", "bcc"):
        value = phases[phase]
        rows.append({
            "row_type": "phase", "system": args.system, "phase": phase,
            "phase_low": "", "phase_high": "",
            "bias_meV_per_atom": value["bias_meV_per_atom"],
            "std_meV_per_atom": value["std_meV_per_atom"],
            "n_snapshots": value["n_snapshots"], "status": "ok",
            "head": value["head"], "source": value["source"],
        })
    low = phases["hcp"]
    high = phases["bcc"]
    rows.append({
        "row_type": "pair", "system": args.system, "phase": "",
        "phase_low": "hcp", "phase_high": "bcc",
        "bias_meV_per_atom": float(high["bias_meV_per_atom"]) - float(low["bias_meV_per_atom"]),
        "std_meV_per_atom": "", "n_snapshots": "", "status": "ok",
        "head": "Domains_Alloy", "source": "phase bias difference",
    })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({"system": args.system, "rows": len(rows), "status": "complete"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
