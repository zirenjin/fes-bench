# fes-bench

`fes-bench` is a configuration-driven benchmark for polymorph Gibbs free-energy
predictions. It separates immutable input data, executable code, and generated
evidence so each reported number has a repository-relative route back to its
script and provenance-bearing raw output.

## Quick start

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -e '.[test]'
bash src/reproduce_all.sh
```

The reproduction command re-normalizes imported historical runs, regenerates
audits and CSV tables, and does not train a model. It needs only the tracked
processed inputs; `data/raw/` is needed only to repeat public-data import.

## Three-layer layout

```text
data/                 public-data contract, ignored archives, tracked processed inputs
src/lib/fes_bench/    reusable data, QH, splitting, evaluation, baseline, and model code
src/experiments/      one experiment-oriented producer per result/experiments/ directory
src/tables/           read-only CSV builders
configs/              system, split, predictor, QH, and audit contracts
result/               generated experiment evidence, tables, and cleaned legacy records
```

`data/processed/` contains the versioned `reference_G.csv`, representative
structures, phase metadata, system contracts, and frozen split JSON files.
`result/experiments/` contains raw `metrics.json` files with provenance,
structured findings, and a directory-level metrics index. Table scripts read
only that result area; they never evaluate a model or accept hand-entered data.

## Reproduce

Run the one-line command in Quick start. The corresponding producer and raw
evidence for each current output are:

| Output | Producer | Raw evidence |
| --- | --- | --- |
| `result/tables/system_inventory.csv` | `src/tables/system_inventory.py` | `result/experiments/data_prep/raw_inventory/`, `reference_statistics/` |
| `result/tables/phase_inventory.csv` | `src/tables/phase_inventory.py` | `result/experiments/data_prep/raw_inventory/`, `quasi_harmonic/` |
| `result/tables/split_definitions.csv` | `src/tables/split_definitions.py` | `result/experiments/data_prep/raw_inventory/` |
| `result/tables/metric_definitions.csv` | `src/tables/metric_definitions.py` | evaluator docstrings under `src/lib/fes_bench/eval/` |
| `result/tables/predictor_comparison_<split>.csv` | `src/tables/predictor_comparison.py --split <split>` | `result/experiments/external_baselines/raw_<split>/` |
| `result/tables/crossing_errors_<split>.csv` | `src/tables/crossing_errors.py --split <split>` | `result/experiments/external_baselines/raw_<split>/` |
| `result/experiments/crossing_reevaluation/` | `src/experiments/crossing_reevaluation.py` | `crossing_reevaluation/raw_runs/` |
| `result/experiments/calibration_window/` | `src/experiments/calibration_window.py` | `calibration_window/raw_runs/` |
| `result/experiments/synthetic_recovery/` | `src/experiments/synthetic_recovery.py` | `synthetic_recovery/raw_runs/` |
| remaining audit directories | matching `src/experiments/<name>.py` | each directory's `metrics.json` provenance index and `findings.meta.json` |

For example, the `interp_const` Tc-error row in
`result/tables/crossing_errors_temp_extrap.csv` is read from
`result/experiments/external_baselines/raw_temp_extrap/interp_const/seed_none/metrics.json`.
Its provenance records the frozen split and source metric hashes.

## Data

The data landing paths, public DOI records, download commands, and SHA-256
checksums are documented in [data/README.md](data/README.md). Downloaded
archives go in ignored `data/raw/`; the checked-in, analysis-ready files go in
`data/processed/`. No public archive is committed.

## Current status

The active reference domain has 4 systems, 9 phases, and 6 within-system phase
pairs: Hf (2/1), SiO2 (3/3), Ti (2/1), and Zr (2/1), where each parenthesis is
phases/pairs. The tabulated reference grids span 693–1,649 points per phase.

Canonical SiO₂ QH production uses Domains_SSE_PBE for all 3 active phases;
1 of 3 (β-quartz, 0.161% imaginary modes at equilibrium volume) passes the 1%
physical reliability gate, while β-cristobalite (4.745%) and β-tridymite
(5.565%) do not. Hf, Ti, and Zr use Domains_Alloy; their current 10 Å QH
inventory has hcp reliable and bcc held for each metal. Ti and Zr
representative structures are complete.

Explicit holds are: CaSiO3 lacks a usable absolute per-phase `G(T,P)` source;
SiO2 C2221 is excluded because no compatible DaRUS reference is available;
and the current SiO2 QH curves fail the physical reliability gate.

## Notes

### Controlled ablation status in the main tables

`E + F_QH` is the no-training controlled baseline in the three v2 main tables:
ΔG-MAE is 215.182 meV/atom (temp_extrap), 169.440 meV/atom (phase_lopo),
and 333.809 meV/atom (system_loso); Hf overlap-T is 412.214 meV/atom.
The E1 polynomial and tlog checkpoints were trained on all four phases and the
full 851–2499 K grid, so they are in-sample and cannot populate learned rows
for the other splits. The T3 SiO2 seed-11 diagnostic was corrected for the
type map and FES label scale and reached 1.861 meV/atom in the training region.
The v2 temp-extrapolation T3 rows are now canonical for three seeds: polynomial
has ΔG-MAE 86.409 meV/atom (skill −18.114), and tlog has 31.557 meV/atom
(skill −5.981), both on the test region. The QH-residual T2 row is also
complete for three seeds: ΔG-MAE 31.239 meV/atom (skill −5.910), sign
accuracy 0.5581, ranking accuracy 0.1045, and mean T_c error 102.36 K;
all five crossing pairs have missed predicted crossings in this test region.

The no-training predictor is expected to have identical metrics on phase_lopo and system_loso: both test-point unions cover the full relative-free-energy grid with the same weighting proportions.

CSV N/A values are explicit: `n/a:no_reference_crossing`,
`n/a:pair_only_predictor`, `n/a:no_training_phase`,
`n/a:degenerate_prediction`, `n/a:missed_crossing`, `n/a:undefined_crossing_slope`, and `n/a:input_unavailable`. Detailed reasons remain
in the matching `.meta.json`. Skill scores now use the zero floor uniformly:
`1 − MAE / MAE_zero`, where `MAE_zero` is the mean reference `|ΔG|` on the
evaluated points. This removes split-dependent denominators and makes every
skill score directly comparable.

Bartel now covers Hf using the Domains_Alloy-computed Hf E0. Its Hf error is
substantially larger than its SiO₂ error; the SiO₂-only table reports 28.2
meV/atom for Bartel ΔG MAE.

Historical delivery records are cleaned and tracked under `result/_legacy/`.
They are not separate generated evidence and are not used by table scripts.
Install the optional test extra shown above before running `python -m pytest`.
The local lightweight run passes 19 tests and skips 2 torch-dependent modules
because the local interpreter has no usable torch. On thu-GenSi, the independent
`/share/jzr/conda-envs/fes-bench-torch` environment runs the full suite
(`21 passed in 43.65s`), including both torch-dependent modules. The CUDA
mat-agent environment is reserved for T3 training.

The β-quartz equilibrium-volume soft-mode diagnostic gives the lowest mode
−1.15 THz at Γ; negative modes span 19 q-points, so they are not Γ-localized.
It is annotated as a possible physical soft mode; the 1%
imaginary-mode threshold may be insensitive to soft-mode-driven transitions.
