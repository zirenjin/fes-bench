"""Audit reference ΔG standard deviations and extrema — plan experiment E9."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from _common import write_audit


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--repo-root", type=Path, default=Path(".")); args = parser.parse_args(); root = args.repo_root.resolve()
    rows, inputs = [], []
    for path in sorted((root / "result/experiments/reference_statistics/raw_runs/reference_statistics/seed_none").glob("metrics.json")):
        inputs.append(path); payload = json.loads(path.read_text(encoding="utf-8"))["metrics"]
        for record in payload.get("records", []): rows.append(record)
    fields = ["split", "fold", "subset", "system", "pair", "n_frames", "temperature_range_K", "mean_eV_per_atom", "std_eV_per_atom", "min_eV_per_atom", "max_eV_per_atom"]
    write_audit(root, "reference_statistics", fields, rows, inputs, [])
    return 0


if __name__ == "__main__": raise SystemExit(main())
