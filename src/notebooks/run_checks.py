"""Execute the table notebook checks locally and write result/tables/_checks.txt."""

from __future__ import annotations

import builtins
import json
import os
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    notebook = root / "src/notebooks/fes_bench_tables.ipynb"
    cells = json.loads(notebook.read_text(encoding="utf-8"))["cells"]
    namespace = {"__name__": "__main__", "display": print, "json": json, "Path": Path, "TABLES": root / "result/tables", "DATA": root / "data/processed"}
    os.environ.setdefault("MPLBACKEND", "Agg")
    for index, cell in enumerate(cells):
        if cell.get("cell_type") != "code" or index in {1, 2, 22}:
            continue
        source = "".join(cell.get("source", []))
        if not source.strip():
            continue
        if index == 3:
            source = source.replace('TABLES = root / "result" / "tables"', f'TABLES = pathlib.Path({str(root / "result/tables")!r})')
            source = source.replace('DATA = root / "data" / "processed"', f'DATA = pathlib.Path({str(root / "data/processed")!r})')
        exec(compile(source, f"{notebook}:cell-{index}", "exec"), namespace)
        if index == 5:
            checks = "\n".join(namespace["CHECKS"]) + "\n"
            # Provenance check is kept outside the notebook so it cannot be
            # accidentally bypassed by display/export code.
            audit = root / "result/experiments/checkpoint_consistency/findings.csv"
            if audit.exists():
                import csv
                rows = list(csv.DictReader(audit.open(encoding="utf-8")))
                bad = [r for r in rows if r.get("status") == "mismatch"]
                checks += ("⚠ checkpoint/head mismatch: %d\n" % len(bad)) if bad else "✓ checkpoint/head provenance consistent for audited canonical rows\n"
            (root / "result/tables/_checks.txt").write_text(checks, encoding="utf-8")
            return 1 if any(line.startswith("⚠") for line in namespace["CHECKS"]) else 0
    raise RuntimeError("notebook consistency-check cell was not executed")


if __name__ == "__main__":
    raise SystemExit(main())
