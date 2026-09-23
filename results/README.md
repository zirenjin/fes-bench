# Results index

`results/raw_runs/` is the only machine-readable input area for table scripts;
each `metrics.json` has a top-level `provenance` block. `results/tables/` holds
CSV summaries and matching `.meta.json` files. `results/audits/` holds
structured findings and provenance for audit experiments. `results/_legacy/`
contains the six historical `phaseX_delivery.md` records and is ignored by
Git; they are preserved locally but deliberately excluded from the remote
repository.

## Reproduce the current tables

From the repository root:

```bash
python3 scripts/reference_delta_g_stats.py --data-root data --splits-root splits \
  --output results/reference_delta_g_stats --include-system-full-grids
PYTHONPATH=. python3 scripts/data/normalize_raw_runs.py --manifest configs/raw_runs.json
PYTHONPATH=scripts/tables python3 scripts/tables/system_inventory.py
PYTHONPATH=scripts/tables python3 scripts/tables/phase_inventory.py
PYTHONPATH=scripts/tables python3 scripts/tables/split_definitions.py
PYTHONPATH=scripts/tables python3 scripts/tables/metric_definitions.py
for split in temp_extrap phase_lopo system_loso; do
  PYTHONPATH=scripts/tables python3 scripts/tables/predictor_comparison.py --split "$split"
  PYTHONPATH=scripts/tables python3 scripts/tables/crossing_errors.py --split "$split"
done
```

The `interp_const` Tc-error cells in
`results/tables/crossing_errors_temp_extrap.csv` come from the corresponding
`results/raw_runs/temp_extrap/interp_const/seed_none/metrics.json`; its
`provenance.input_metrics_sha256` identifies the original
`results/external_baselines/metrics.json`, and the split/config hashes identify
the frozen evaluator inputs.

## Reproduce the audits

```bash
for audit in crossing_reevaluation calibration_window synthetic_recovery skill_floor mae_convention split_sign_structure imaginary_modes representative_structures reference_statistics leakage reproduction_compare; do
  PYTHONPATH=scripts/audits python3 "scripts/audits/${audit}.py"
done
```

`skill_floor/findings.csv` retains the denominator used by each skill score in
the `floor_predictor` column; `constant_delta_g` and `global_mean_delta_g` are
not silently unified. E1 rows additionally carry `aggregation` values
`include_all_pairs` and `exclude_no_reference_crossing`; negative skills are
explicitly marked by `skill_status=NEGATIVE`.

## CSV missing-value vocabulary

CSV cells use explicit markers while the matching `.meta.json` keeps the
detailed field-level reason:

- `n/a:no_reference_crossing` — the pair has no sign change on the frozen grid.
- `n/a:pair_only_predictor` — the predictor reports pair-level ΔG only and has
  no phase-level scalar G output.
- `n/a:no_training_phase` — the held-out fold lacks a training phase/system
  required for that value.

The current generated tables contain 33, 15, and 32 cells of these three
types, respectively. `split_definitions.csv` derives `crossing_in_test` from
the frozen split index and reference crossing JSON without modifying either;
all `temp_extrap` reference crossings are asserted to be in the test region.

## Known issues

`constant_delta_g` is structurally unavailable on LOPO/LOSO because the held
out phase/system has no training pair labels. Consequently the three main
comparison tables currently use different skill denominators: `constant_delta_g`
for `temp_extrap`, and `global_mean_delta_g` for LOPO/LOSO. This is recorded,
not unified in this round; unifying it requires a new recomputation.

The historical 2.56 meV/atom source is closed in
`results/audits/mae_convention/findings.meta.json`: it used four-phase gauge-C,
one quartz–C2221 pair, the full grid, and a different checkpoint, so it does
not affect current conclusions. The SiO₂ `qh_comparison.json` copied from the
isolated run is recorded with source path and SHA-256 in the QH raw-run
provenance.

## Current status and missing cells

The generated `.meta.json` files are authoritative for missing cells and state
whether a value is unavailable, structurally N/A, or not implemented. Large
raw data and external checkpoints are intentionally not included; their
download/provenance contracts are described in `data/README.md`.
