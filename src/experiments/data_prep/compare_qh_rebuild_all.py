"""Build the combined SiO2 and metal QH rebuild comparison — plan Table 2b."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, default=Path("result/experiments/quasi_harmonic_rebuild"))
    args = parser.parse_args()
    root = args.repo_root.resolve(); output = args.output if args.output.is_absolute() else root / args.output
    output.mkdir(parents=True, exist_ok=True)
    jobs = [("sio2", "result/experiments/quasi_harmonic_ideal_rebuild", "sio2"), ("hf", "result/experiments/quasi_harmonic_10a", "hf"), ("ti", "result/experiments/quasi_harmonic_10a", "ti"), ("zr", "result/experiments/quasi_harmonic_10a", "zr")]
    rows = []
    for system, new_root, subdir in jobs:
        target = output / subdir
        subprocess.run([sys.executable, "src/experiments/data_prep/compare_qh_rebuild.py", "--repo-root", str(root), "--system", system, "--new-root", new_root, "--output", str(target)], check=True)
        with (target / "findings.csv").open(encoding="utf-8", newline="") as handle:
            rows.extend(csv.DictReader(handle))
    fields = list(rows[0]) if rows else ["system", "phase"]
    with (output / "findings.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader(); writer.writerows(rows)
    (output / "findings.meta.json").write_text(json.dumps({"git_commit": subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip(), "rows": len(rows), "sources": ["result/experiments/quasi_harmonic_ideal_rebuild", "result/experiments/quasi_harmonic_10a"]}, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
