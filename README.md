# fes-bench

## Reorganized table + audit layout

The repository is organized around configuration-driven raw runs, tables, and
audits. Paths in the reproducibility commands below are relative to the
repository root; public datasets and external checkpoints are not committed.

```text
configs/systems/       system contracts and reference levels
configs/splits/        split-generation parameters
configs/predictors/    predictor contracts
data/raw/              ignored public archives
data/processed/        ignored derived data
scripts/data/           normalization and download entry points
scripts/qh/             configured QH wrapper
scripts/eval/           configured evaluator wrapper
scripts/tables/         CSV + meta.json table builders
scripts/audits/         findings.csv + findings.meta.json audits
results/raw_runs/       normalized JSON runs with provenance
results/tables/         reproducible CSV summaries
results/audits/         structured audit findings
results/_legacy/        local-only archived phase records (ignored by Git)
```

## Reproduce every current conclusion

```bash
python3 scripts/reference_delta_g_stats.py --data-root data --splits-root splits \
  --output results/reference_delta_g_stats --include-system-full-grids
PYTHONPATH=. python3 scripts/data/normalize_raw_runs.py --manifest configs/raw_runs.json
for table in system_inventory phase_inventory split_definitions metric_definitions; do
  PYTHONPATH=scripts/tables python3 "scripts/tables/${table}.py"
done
for split in temp_extrap phase_lopo system_loso; do
  PYTHONPATH=scripts/tables python3 scripts/tables/predictor_comparison.py --split "$split"
  PYTHONPATH=scripts/tables python3 scripts/tables/crossing_errors.py --split "$split"
done
for audit in crossing_reevaluation calibration_window synthetic_recovery skill_floor mae_convention split_sign_structure imaginary_modes representative_structures reference_statistics leakage reproduction_compare; do
  PYTHONPATH=scripts/audits python3 "scripts/audits/${audit}.py"
done
```

The output index and provenance route for each conclusion is
`results/README.md`; public data URLs, SHA-256 values, and landing paths are in
`data/README.md` and `data/*/download.json`.

Reproducible infrastructure for benchmarking polymorph Gibbs free-energy
predictions. This is an independent repository: `deepmd-kit` is an external
runtime dependency used only by later training and QH phases.

## Current status

The repository now contains a reproducible infrastructure layer through the
data audit, QH diagnostic, frozen splitting, unified evaluation, external
baseline fixtures, and synthetic recovery stages. It does **not** claim a
complete production benchmark: Ti/Zr representative structures, CaSiO3
absolute per-phase `G(T,P)`, valid QH curves, and checkpoint re-evaluations
remain explicit holds in their respective delivery records.

The active minimal reference domain is Hf plus the three DaRUS-consistent
SiO2 phases (β-quartz, β-cristobalite, and P6₃/mmc tridymite). C222₁ is
excluded rather than replaced with the incompatible COD representative.

## Data contract

Each phase is stored at `data/<system>/<phase>/`:

```
structure.extxyz
reference_G.csv
meta.json
```

`reference_G.csv` has exactly these required columns:

```
T_K,P_GPa,G_eV_per_atom,level,source
```

`meta.json` must provide the declared electronic-structure provenance and
method fields. `data/<system>/system.json` records the system-wide phase list,
type map, reference grid, and truth level. The loader validates that all
temperature/pressure/free-energy values are finite and that every tabulated
phase is named by `system.json`.

The public loading API is:

```python
from fes_bench.data import load
phase = load("sio2", "beta_quartz")
```

By default it reads this repository's `data/` directory. Callers may pass a
different `data_root`, or set `FES_BENCH_DATA_ROOT`, for a copied dataset.

## Commands

All command-line modules take a configuration file:

```bash
python -m fes_bench.data.audit --config configs/audit.yaml
python -m fes_bench.splits.make --config configs/splits.yaml
python -m fes_bench.splits.verify --config configs/split_verify.yaml
python -m fes_bench.eval.run --predictor reference --split splits/temp_extrap.json \
  --data-root data --output-root results
python -m fes_bench.eval.plot --predictor reference --split splits/phase_lopo.json \
  --data-root data --output results/reference/phase_lopo/delta_g_curves.png
python -m fes_bench.synthetic.e3_recovery --config configs/synthetic/e3_recovery.yaml
pytest
```

## Development environment

The test dependency is the optional `test` extra declared in `pyproject.toml`:

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -e '.[test]'
python -m pytest
```

The current local interpreter cannot run pytest because that package is not
installed (`No module named pytest`). Install the optional extra in a virtual
environment as above; no dependency installation was performed during the
legacy-record migration.

The supplied YAML files use JSON syntax, which is valid YAML and keeps Phase 0
free of a parser dependency. Later phases may use normal YAML when their
environment explicitly provides a YAML parser.

## External-baseline note

The Bartel 2018 baseline uses only atomic volume as its structural input once
composition is fixed: reduced mass is then identical for every polymorph.
Consequently it is deliberately a coarse temperature-dependent comparator,
not a phase-sensitive free-energy model.

## Repository layout

```
fes_bench/data/       data download, parsing, schema, representative structures
fes_bench/qh/         phonopy + DPA quasi-harmonic implementation (Phase 2)
fes_bench/splits/     frozen split generation (Phase 3)
fes_bench/eval/       evaluation (Phase 3)
fes_bench/baselines/  external baselines (Phase 5)
fes_bench/models/     deepmd FES predictor wrappers (Phase 4)
fes_bench/synthetic/  synthetic recovery tests (Phase 5)
configs/              configuration files
data/                 generated datasets and tracked manifests
splits/               frozen split JSON (Phase 3)
results/              generated results and tracked manifests
tests/                executable tests
```

## Delivery records and deviations

Historical phase delivery records are preserved locally under
`results/_legacy/` and deliberately excluded from Git. The tracked
reproduction entry point is `results/README.md`; it indexes the structured
raw runs, tables, audits, and their missing-cell metadata.
- `results/completion_audit_20260921.json` — machine-readable artifact audit;
  `pass_with_holds` means infrastructure checks pass while the declared source
  and physical holds remain active.

The independent repository is initialized. Frozen split JSON files record the
real provenance commit `b5b0b2e679b2da66a593d713afc677f1adeab61e`, and carry
embedded SHA-256 values that are revalidated by `fes_bench.splits.verify`; the
passing audit is tracked at `results/phase3/split_integrity.json`.
