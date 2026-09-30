"""Execute the table notebook checks locally and write result/tables/_checks.txt."""

from __future__ import annotations

import builtins
import csv
import json
import os
import re
from pathlib import Path


def _fallback_checks(root: Path) -> list[str]:
    """Run the notebook's consistency assertions without optional pandas/matplotlib."""
    out: list[str] = []
    def rows(name):
        with (root / "result/tables" / f"{name}.csv").open(encoding="utf-8", newline="") as fh:
            return list(csv.DictReader(fh))
    for split in ("temp_extrap", "phase_lopo", "system_loso"):
        rs = rows(f"predictor_comparison_{split}")
        c = next((r for r in rs if r.get("predictor") == "constant_delta_g"), None)
        floor_fields = {"delta_G_MAE_eV_per_atom", "delta_G_RMSE_eV_per_atom", "sign_accuracy", "Tc_error_K", "Tc_err_from_dG_K", "false_crossings", "missed_crossings", "skill_score"}
        if split in {"phase_lopo", "system_loso"} and c and any(c.get(k) not in {"n/a:no_training_phase", ""} for k in floor_fields):
            out.append(f"⚠ {split}: constant ΔG is not marked n/a:no_training_phase")
        else: out.append(f"✓ {split}: constant ΔG floor labels are consistent")
        z = next((r for r in rs if r.get("predictor") == "zero"), None)
        if z and any("Tc" in k and v not in {"", "n/a:degenerate_prediction"} for k,v in z.items()):
            out.append(f"⚠ {split}: zero crossing fields are not degenerate_prediction")
        else: out.append(f"✓ {split}: zero crossing fields use n/a:degenerate_prediction")
        required={"predictor","delta_G_MAE_eV_per_atom","sign_accuracy","skill_score","floor_predictor"}
        missing=sorted(required-set(rs[0]) if rs else required)
        out.append(f"⚠ {split}: main table missing columns {missing}" if missing else f"✓ {split}: main table required columns present")
        empties=[r.get("predictor") for r in rs if r.get("predictor") not in {"", "self_check_reference"} and any(v=="" for k,v in r.items() if k not in {"notes"})]
        out.append(f"⚠ {split}: unmarked empty cells {empties}" if empties else f"✓ {split}: non-self rows have no unmarked empty cells")
    inv={r["system"]:r for r in rows("system_inventory")}; phases=rows("phase_inventory")
    for system, row in inv.items():
        pr=[p for p in phases if p.get("system")==system]
        true=sum(p.get("qh_reliable")=="True" for p in pr)
        out.append(f"✓ {system}: qh_reliable_count={row.get('qh_reliable_count')} agrees with phase_inventory ({true})" if str(true)==str(row.get("qh_reliable_count")) else f"⚠ {system}: qh_reliable_count mismatch")
    splitrows=rows("split_definitions")
    bad=[r for r in splitrows if r.get("split")=="temp_extrap" and r.get("crossing_in_test") and not all(x.endswith("=test") for x in r["crossing_in_test"].split(";"))]
    out.append("⚠ temp_extrap has crossings outside test" if bad else "✓ temp_extrap crossings all in test")
    bad=[p.get("phase") for p in phases if p.get("qh_reliable") and not p.get("imaginary_fraction")]
    out.append(f"⚠ phase_inventory missing imaginary_fraction {bad}" if bad else "✓ phase_inventory reliable rows have imaginary_fraction")
    bad=[p.get("phase") for p in phases if p.get("qh_supercell_shortest_edge_A") and float(p["qh_supercell_shortest_edge_A"])<10]
    out.append(f"⚠ QH supercell edge <10 A {bad}" if bad else "✓ QH supercell shortest edges meet 10 A")
    expected_folds = {"hf", "ti", "zr"}
    # Main system_loso CSVs must not retain the obsolete v1 SiO2 fold.  The
    # separately named overlap_T support CSV is intentionally Hf-only.
    for table_name in ("predictor_comparison_system_loso_folds", "crossing_errors_system_loso"):
        fold_table = root / "result/tables" / f"{table_name}.csv"
        if not fold_table.exists():
            out.append(f"⚠ {table_name} is missing")
            continue
        fold_rows = list(csv.DictReader(fold_table.open(encoding="utf-8", newline="")))
        observed = {r.get("fold", "") for r in fold_rows if r.get("eval_subset", "full") == "full" and r.get("fold", "")}
        out.append(f"⚠ {table_name} fold set is {sorted(observed)}, expected ['hf', 'ti', 'zr']" if observed != expected_folds else f"✓ {table_name} fold set is exactly {{hf, ti, zr}}")
    return out


def _non_english_checks(root: Path) -> list[str]:
    """The checked-in executable notebook must contain no CJK text."""
    pattern = re.compile(r"[\u3000-\u303f\u3400-\u4dbf\u4e00-\u9fff\uf900-\ufaff]")
    warnings = []
    for path in sorted((root / "src/notebooks").glob("*.ipynb")):
        if pattern.search(path.read_text(encoding="utf-8")):
            warnings.append(f"⚠ non-English text in {path.relative_to(root)}")
    return warnings


def _metric_csv_empty_checks(root: Path) -> list[str]:
    """Reject unmarked blanks in fold-level and crossing-level metric CSVs."""
    out: list[str] = []
    paths = sorted((root / "result/tables").glob("predictor_comparison_*_folds.csv"))
    paths += sorted((root / "result/tables").glob("crossing_errors_*.csv"))
    for path in paths:
        with path.open(encoding="utf-8", newline="") as fh:
            for number, row in enumerate(csv.DictReader(fh), start=2):
                empty = [key for key, value in row.items() if value == ""]
                if empty:
                    ident = row.get("predictor", row.get("pair", "row"))
                    out.append(f"⚠ {path.relative_to(root)}:{number} ({ident}) has unmarked empty cells: {','.join(empty)}")
    return out


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
        try:
            exec(compile(source, f"{notebook}:cell-{index}", "exec"), namespace)
        except ModuleNotFoundError as exc:
            if exc.name not in {"pandas", "matplotlib", "numpy"}:
                raise
            checks = _fallback_checks(root)
            checks.extend(_metric_csv_empty_checks(root))
            checks.extend(_non_english_checks(root))
            (root / "result/tables/_checks.txt").write_text("\n".join(checks) + "\n", encoding="utf-8")
            return 1 if any(line.startswith("⚠") for line in checks) else 0
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
                try:
                    from experiments.checkpoint_consistency import policy_warnings
                except ModuleNotFoundError:
                    from src.experiments.checkpoint_consistency import policy_warnings
                policy_bad = policy_warnings(root)
                checks += "".join(f"⚠ {item}\n" for item in policy_bad)
            checks += "".join(f"{line}\n" for line in _non_english_checks(root))
            checks += "".join(f"{line}\n" for line in _metric_csv_empty_checks(root))
            (root / "result/tables/_checks.txt").write_text(checks, encoding="utf-8")
            return 1 if any(line.startswith("⚠") for line in checks.splitlines()) else 0
    raise RuntimeError("notebook consistency-check cell was not executed")


if __name__ == "__main__":
    raise SystemExit(main())
