"""Execute the table notebook checks locally and write result/tables/_checks.txt."""

from __future__ import annotations

import re


def _canonical_pair_name(value: str) -> str:
    value = value.strip()
    return value.replace("bcc_minus_hcp", "hcp_minus_bcc").replace("bcc-hcp", "hcp_minus_bcc")


def _pair_names(path: Path, *, canonicalize: bool = True) -> set[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        field = next((name for name in ("pair", "pair_name", "pair_id", "phase_pair") if name in (reader.fieldnames or [])), None)
        if field is None:
            return set()
        values = {row[field].strip() for row in reader if row.get(field, "").strip()}
        return {_canonical_pair_name(value) for value in values} if canonicalize else values


def _additional_checks(root: Path) -> list[str]:
    checks: list[str] = []
    for relative in ("README.md", "result/README.md"):
        path = root / relative
        bad = [index for index, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1) if re.search(r"[\u3400-\u9fff]", line)]
        checks.append(f"⚠ {relative}: non-English Han characters at lines {bad}" if bad else f"✓ {relative}: non-English check passed")

    detail_files = sorted((root / "result/experiments/t1_qh").glob("**/pair_details_*.csv"))
    for path in detail_files:
        raw_names = _pair_names(path, canonicalize=False)
        noncanonical = sorted(name for name in raw_names if "bcc_minus_hcp" in name or name == "bcc-hcp")
        checks.append(f"⚠ {path.relative_to(root)}: metal pair names are not hcp_minus_bcc: {noncanonical}" if noncanonical else f"✓ {path.relative_to(root)}: pair names use canonical direction")

    for split in ("temp_extrap", "phase_lopo", "system_loso"):
        predictor_files = sorted((root / "result/tables").glob(f"predictor_comparison_{split}*.csv"))
        predictor_files += [path for path in sorted((root / "result/tables").glob(f"crossing_errors_{split}*.csv")) if "_overlap_T" not in path.name]
        t1_files = [path for path in detail_files if split in path.name or split in str(path.parent)]
        sets: list[tuple[str, set[str]]] = []
        for path in predictor_files + t1_files:
            names = _pair_names(path)
            if names:
                sets.append((str(path.relative_to(root)), names))
        if not sets:
            checks.append(f"ℹ {split}: no predictor/t1 pair-detail files found")
            continue
        # Crossing-error tables intentionally list only pairs with a reference
        # crossing, while T1 detail tables also retain no-reference pairs.
        expected = next((values for name, values in sets if "crossing_errors_" in name), sets[0][1])
        mismatches = [name for name, values in sets[1:] if not expected.issubset(values) and not values.issubset(expected)]
        checks.append(f"⚠ {split}: predictor and t1 pair-name sets differ: {mismatches}" if mismatches else f"✓ {split}: predictor and t1 pair-name sets agree")
    return checks

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
    bad_slopes = []
    for crossing_path in sorted((root / "data/processed").glob("*/reference_crossings.json")):
        payload = json.loads(crossing_path.read_text(encoding="utf-8"))
        for pair_name, pair in payload.get("pairs", {}).items():
            for index, crossing in enumerate(pair.get("crossings", []), start=1):
                slope = crossing.get("slope_eV_per_atom_per_K")
                r_squared = crossing.get("r_squared")
                if not isinstance(slope, (int, float)) or slope <= 0 or not isinstance(r_squared, (int, float)) or r_squared < 0.95:
                    bad_slopes.append(f"{crossing_path.parent.name}:{pair_name}#{index}")
    out.append(f"⚠ reference crossing slope audit failed: {bad_slopes}" if bad_slopes else "✓ reference crossing slopes are positive with R²≥0.95")
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


def _crossing_slope_consistency_checks(root: Path) -> list[str]:
    """Ensure figure-1 fits, reference crossings, and Table 6 use one slope."""
    tolerance = 1e-10
    reference: dict[tuple[str, str, int], float] = {}
    for crossing_path in sorted((root / "data/processed").glob("*/reference_crossings.json")):
        payload = json.loads(crossing_path.read_text(encoding="utf-8"))
        system = str(payload.get("system", crossing_path.parent.name))
        for pair_name, pair in payload.get("pairs", {}).items():
            for index, crossing in enumerate(pair.get("crossings", []), start=1):
                slope = crossing.get("slope_eV_per_atom_per_K")
                if isinstance(slope, (int, float)):
                    reference[(system, str(pair_name), index)] = abs(float(slope))

    def lookup(system: str, pair_name: str, index: int) -> float | None:
        value = reference.get((system, pair_name, index))
        if value is not None:
            return value
        if "_minus_" in pair_name:
            left, right = pair_name.split("_minus_", 1)
            return reference.get((system, f"{right}_minus_{left}", index))
        return None

    warnings: list[str] = []
    table_failures: list[str] = []
    for path in sorted((root / "result/tables").glob("crossing_errors_*.csv")):
        with path.open(encoding="utf-8", newline="") as handle:
            for number, row in enumerate(csv.DictReader(handle), start=2):
                pair = row.get("pair", "")
                if ":" not in pair:
                    continue
                system, pair_name = pair.split(":", 1)
                try:
                    index = int(row.get("crossing_index", "1"))
                    observed = float(row.get("crossing_slope_eV_per_atom_per_K", ""))
                except (TypeError, ValueError):
                    continue
                expected = lookup(system, pair_name, index)
                if expected is not None and abs(abs(observed) - expected) > tolerance:
                    table_failures.append(f"{path.name}:{number}")
    if table_failures:
        warnings.append(f"⚠ crossing_errors slope mismatch with reference_crossings: {table_failures}")
    else:
        warnings.append("✓ crossing_errors slopes match reference_crossings (absolute sign convention)")

    # Reconstruct the ΔG-at-crossing inputs from the raw runs and verify the
    # table's derived Tc_err_from_dG values use the same fitted slope.
    raw_dg: dict[tuple[str, str, str, str, str, int], float] = {}

    def add_raw(split: str, predictor: str, subset: str, fold: str, metrics: dict) -> None:
        pairs = metrics.get("pairs", {}) if isinstance(metrics, dict) else {}
        if not isinstance(pairs, dict):
            return
        for pair_name, record in pairs.items():
            values = record.get("delta_G_MAE_at_crossing_eV_per_atom", [])
            if not isinstance(values, list):
                continue
            for index, value in enumerate(values, start=1):
                if isinstance(value, (int, float)):
                    raw_dg[(split, predictor, subset, fold, str(pair_name), index)] = float(value)

    for split in ("temp_extrap", "phase_lopo", "system_loso"):
        patterns = [f"result/experiments/external_baselines/raw_{split}/*/*/metrics.json"]
        if split == "temp_extrap":
            patterns.append("result/experiments/t3_temp_extrap/raw_runs/*/*/metrics.json")
        for pattern in patterns:
            for raw_path in sorted(root.glob(pattern)):
                try:
                    body = json.loads(raw_path.read_text(encoding="utf-8"))
                    predictor = raw_path.parts[-3]
                    metrics = body.get("metrics", {})
                    folds = metrics.get("folds", {}) if isinstance(metrics, dict) else {}
                    if isinstance(folds, dict) and folds:
                        for fold_name, fold_metrics in folds.items():
                            add_raw(split, predictor, "main", str(fold_name), fold_metrics)
                            if isinstance(fold_metrics, dict) and isinstance(fold_metrics.get("overlap_T"), dict):
                                add_raw(split, predictor, "overlap_T", str(fold_name), fold_metrics["overlap_T"])
                    else:
                        add_raw(split, predictor, "main", "all", metrics)
                except (OSError, json.JSONDecodeError):
                    continue

    tc_failures: list[str] = []
    for path in sorted((root / "result/tables").glob("crossing_errors_*.csv")):
        split_name = path.name[len("crossing_errors_") : -len(".csv")]
        subset_default = "overlap_T" if split_name.endswith("_overlap_T") else "main"
        split = split_name.removesuffix("_overlap_T")
        with path.open(encoding="utf-8", newline="") as handle:
            for number, row in enumerate(csv.DictReader(handle), start=2):
                tc_value = row.get("Tc_err_from_dG_K", "")
                try:
                    observed_tc = float(tc_value)
                    index = int(row.get("crossing_index", "1"))
                except (TypeError, ValueError):
                    continue
                pair = row.get("pair", "")
                if ":" not in pair:
                    continue
                system, pair_name = pair.split(":", 1)
                slope = lookup(system, pair_name, index)
                if slope is None:
                    continue
                subset = row.get("subset") or subset_default
                fold = row.get("fold") or "all"
                dg = raw_dg.get((split, row.get("predictor", ""), subset, fold, pair_name, index))
                if dg is None and "_minus_" in pair_name:
                    left, right = pair_name.split("_minus_", 1)
                    dg = raw_dg.get((split, row.get("predictor", ""), subset, fold, f"{right}_minus_{left}", index))
                if dg is not None and abs(observed_tc - abs(dg) / slope) > tolerance:
                    tc_failures.append(f"{path.name}:{number}")
    if tc_failures:
        warnings.append(f"⚠ Tc_err_from_dG mismatch with fitted slopes: {tc_failures}")
    else:
        warnings.append("✓ Tc_err_from_dG values match fitted reference slopes")

    figure_path = root / "result/figures/figure1/panel_b.csv"
    figure_failures: list[str] = []
    if figure_path.exists():
        with figure_path.open(encoding="utf-8", newline="") as handle:
            for number, row in enumerate(csv.DictReader(handle), start=2):
                system = row.get("system", "")
                low = row.get("low", "")
                high = row.get("high", "")
                pair_name = None
                for (candidate_system, candidate_pair, _), _value in reference.items():
                    if candidate_system != system:
                        continue
                    if set(candidate_pair.split("_minus_")) == {low, high}:
                        pair_name = candidate_pair
                        break
                try:
                    observed = abs(float(row.get("slope_meV_per_atom_per_K", ""))) * 1e-3
                except (TypeError, ValueError):
                    continue
                expected = lookup(system, pair_name or "", 1)
                if expected is not None and abs(observed - expected) > tolerance:
                    figure_failures.append(f"panel_b.csv:{number}")
    if figure_failures:
        warnings.append(f"⚠ figure 1 slope mismatch with reference_crossings: {figure_failures}")
    elif figure_path.exists():
        warnings.append("✓ figure 1 panel_b slopes match reference_crossings")
    return warnings


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
            checks = _fallback_checks(root) + _additional_checks(root)
            checks.extend(_metric_csv_empty_checks(root))
            checks.extend(_crossing_slope_consistency_checks(root))
            checks.extend(_non_english_checks(root))
            (root / "result/tables/_checks.txt").write_text("\n".join(checks) + "\n", encoding="utf-8")
            return 1 if any(line.startswith("⚠") for line in checks) else 0
        if index == 5:
            checks = "\n".join(namespace["CHECKS"]) + "\n"
            checks += "\n".join(_additional_checks(root)) + "\n"
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
            checks += "".join(f"{line}\n" for line in _crossing_slope_consistency_checks(root))
            (root / "result/tables/_checks.txt").write_text(checks, encoding="utf-8")
            return 1 if any(line.startswith("⚠") for line in checks.splitlines()) else 0
    raise RuntimeError("notebook consistency-check cell was not executed")


if __name__ == "__main__":
    raise SystemExit(main())
