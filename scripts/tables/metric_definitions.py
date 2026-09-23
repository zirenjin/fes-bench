"""Metric names and definitions extracted from the evaluator contract — plan Table 4."""

from __future__ import annotations

import argparse
import ast
from pathlib import Path

from common import csv_write, inventory, meta_write


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--output-root", type=Path, default=Path("results/tables"))
    args = parser.parse_args(); root = args.repo_root.resolve(); output = args.output_root if args.output_root.is_absolute() else root / args.output_root
    evaluator = root / "fes_bench/eval/run.py"
    docstring = ast.get_docstring(ast.parse(evaluator.read_text(encoding="utf-8"))) or ""
    rows = []
    in_section = False
    for line in docstring.splitlines():
        if line.strip() == "Metric definitions:":
            in_section = True
            continue
        if not in_section or not line.strip().startswith("- "):
            continue
        metric, definition, unit, direction = [item.strip() for item in line.strip()[2:].split("|")]
        rows.append({"metric": metric, "definition": definition, "unit": unit, "direction": direction})
    csv_path = output / "metric_definitions.csv"; csv_write(csv_path, ["metric", "definition", "unit", "direction"], rows)
    meta_write(root, csv_path, [root / "results/raw_runs/_inventory/inventory/seed_none/metrics.json", evaluator], [])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
