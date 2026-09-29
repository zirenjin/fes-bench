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

SiO2 QH production has diagnostic outputs for 3 of 3 active SiO2 phases; 0 of
3 pass the physical reliability gate because of imaginary-mode diagnostics.
Hf, Ti, and Zr have no current QH production `metrics.json` in this result
layer, which is data/provenance unavailable here rather than a statement that
their representative structures are missing. Ti and Zr representative
structures are complete.

Explicit holds are: CaSiO3 lacks a usable absolute per-phase `G(T,P)` source;
SiO2 C2221 is excluded because no compatible DaRUS reference is available;
and the current SiO2 QH curves fail the physical reliability gate.

## Notes

### 主表受控消融块现状

三个冻结切分上目前没有合法的 learned predictor 结果：E1 的 polynomial 与 tlog checkpoint 均为全四相、851–2499 K 全温区训练，属于 in-sample 评测，不能接入 temp_extrap、phase_lopo 或 system_loso 主表。要填入这些行，必须分别在每个冻结切分的训练集上训练对应的 polynomial 与 tlog 模型，并只在该切分测试区重评；本轮不训练。

无训练 predictor 在 phase_lopo 与 system_loso 上指标相同是预期的：两者的测试点并集都是全部相对的全温区，权重比例一致。

CSV N/A values are explicit: `n/a:no_reference_crossing`,
`n/a:pair_only_predictor`, `n/a:no_training_phase`,
`n/a:degenerate_prediction`, and `n/a:input_unavailable`. Detailed reasons remain
in the matching `.meta.json`. `constant_delta_g` is structurally unavailable
on LOPO/LOSO, so skill denominators remain `constant_delta_g` for temp-extrap
and `global_mean_delta_g` for LOPO/LOSO; this migration records rather than
recomputes that difference.

Historical delivery records are cleaned and tracked under `result/_legacy/`.
They are not separate generated evidence and are not used by table scripts.
The local interpreter currently has no `pytest`; install the optional test
extra shown above before running `python -m pytest`. The local lightweight run
passes 19 tests and skips 2 torch-dependent modules because torch is unavailable;
the isolated thu-GenSi run passes 19 and skips the same 2 modules because its
torch installation lacks `libtorch_global_deps.so`.
