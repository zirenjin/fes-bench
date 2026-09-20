# fes-bench

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

Each phase record names its artifacts, actual command output, acceptance
evidence, and deviations:

- `results/phase1_delivery.md` — source/representative holds and data audit.
- `results/phase2_delivery.md` — QH diagnostic gate; no invalid curve is used
  as a production physical baseline.
- `results/phase3_delivery.md` — frozen splits, reference/noise acceptance,
  V100 evaluation, tables, and curves.
- `results/phase4_delivery.md` — scoped `feat/fes-head` physics-baseline
  interface and focused frozen-export evidence.
- `results/phase5_delivery.md` — external baseline fixtures, E3 recovery, and
  the non-fabricated E1/E2 availability hold.
- `results/completion_audit_20260921.json` — machine-readable artifact audit;
  `pass_with_holds` means infrastructure checks pass while the declared source
  and physical holds remain active.

The independent repository currently has no initial Git commit. Frozen split
JSON files therefore record `git_commit: unavailable`, but carry embedded
SHA-256 values that are revalidated by `fes_bench.splits.verify`; the passing
audit is tracked at `results/phase3/split_integrity.json`.
