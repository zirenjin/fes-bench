# FES-Bench v0.1 (2026-11-02)

**A benchmark for phase-relative harmonic free energies predicted by machine-learning interatomic potentials (MLIPs).** It measures how accurate those free energies must be to place a polymorph transition, which capabilities ordinary accuracy metrics miss, and where the practical workflow fails.

## Contents

| Path | What it is |
|---|---|
| `manifest/ENTRIES.csv` | 521 entries in 239 formulas. Each has its reference (Alexandria PBE phonon benchmark, CC BY 4.0) and tier, plus validity flags: the harmonic gate at −0.05 THz for DFT, F0 and MatterSim, the U-scheme / magnetic flags, and the F0 relaxed-pipeline status. Also stored: structure metadata and the stored per-entry wall time. |
| `manifest/PAIRS.csv` | 351 same-formula pairs: the reference ΔG at 300/600/1000/1500 K, the first reference crossing in 50–1500 K, ΔS at the crossing, the ±50 K budget ε₅₀, and common-validity. |
| `manifest/SPLITS.csv` | The formula-grouped outer folds (3 repeats) used by the intervention study. |
| `predictions/*.npz` | Stored curves on 0–1500 K (50 K grid): `E0` + `F`, or `G`. Includes the reference, F0 (DPA-3.1-3M), MatterSim-v1.0.0-5M and the frozen equal-weight average, each under S0 (no reference information) and S1 (static ΔE₀ given). Also the out-of-fold intervention routes A–C. |
| `score.py` | The single scoring entry point: `python score.py predictions/X.npz`. It writes `scores/X__{all,common_valid}.json`, every metric with its weighting and a 95% formula-cluster bootstrap. |
| `make_main.py` | The main table (`tables/main_table.{csv,tex}`, plus the common-valid version) and the four main figures (`figures/fig1–4`). |
| `audit_numerics.py`, `AUDIT.md` | Numerical robustness checks and the claim–evidence register. |
| `build.py` | Rebuilds `manifest/` and `predictions/` from the frozen source outputs listed in it. |

## Tasks

- **Main: fixed geometry, harmonic.** Fixed DFT geometry, harmonic, the same stored displacements and phonopy path for every model.
- **Pipeline diagnostic: relaxed geometry.** The MLIP-relaxed pipeline, available for F0 only. Its failures are recorded in `ENTRIES.csv`.
- **QHA diagnostic: HfO₂, lenient tier.** Reported separately and never used for model ranking. The data stay local until the authors give permission.

## Status

All results are development evidence on data analysed before the protocols were frozen. MatterSim is a model replication. **There is no independent confirmation set.** That is the entry condition for any next-round method claim.

## Scoring rules (after external review)

- **Reference-invalid pairs are not scored.** 31 pairs touch an entry whose reference minimum frequency is below -0.05 THz; they have no trustworthy target.
- **Accuracy uses model-valid pairs.** A model-invalid pair counts as a failure instead.
- **Success is over all 320 reference-valid pairs.** It requires the pair to be model-valid and correct.
- **Paired comparisons use the 219 common-valid pairs.** These are valid for DFT and both models.
- **Truncated-harmonic curves are diagnostic only.** They appear in `tables/main_table_all_diagnostic.tex`.
- **The energy budget and the temperature outcome are separate metrics.**
  - In-eps50 is |e| <= |dS_ref| x 50 K, an energy criterion.
  - Hit +-50 K means a model crossing lies within 50 K of the reference crossing.
- **Crossings are pairwise.** `crossing_on_lower_envelope` marks the ones that are also candidate-set phase boundaries (76 of 93).
- **Abstention is a leave-one-formula-out empirical-quantile rule.** It reports a decided fraction, not a coverage guarantee.

## Reproduction levels

- **Rescoring** (`score.py`, `make_main.py`, `audit_numerics.py`) uses only `manifest/` and `predictions/`. Predictions cover all three intervention repeats.
- **Rebuilding** (`build.py`) needs the source analysis outputs listed in it.
- **Not reproducible from the stored G(T) curves alone:** the tail audit, the intervention training diagnostics and the curvature-correction diagnostic. They come from those source analyses.

## Reproduce

```
python build.py
python score.py predictions/*.npz     # writes scores/<name>__{model_valid,common_valid,all}.json
python make_main.py
python audit_numerics.py
```

None of these steps needs model inference or DFT.
