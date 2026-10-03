# Single-sided QH diagnostic

This diagnostic compares anchored thermal shapes on the QH-reliable side only. It does not evaluate phase pairs, ΔG, ranking, or T_c.

## Results

| Phase | Absolute-G RMSE | One-sided QH RMSE | R = QH / absolute | Improvement | Absolute slope | QH slope |
|---|---:|---:|---:|---:|---:|---:|
| Hf hcp | 421.351 ± 0.671 | 28.069 ± 0.001 | 0.067 ± 0.000 | +93.3% | -1.53214 | 0.10651 |
| Ti hcp | 214.456 ± 0.131 | 52.499 ± 0.000 | 0.245 ± 0.000 | +75.5% | -0.69241 | 0.15826 |
| Zr hcp | 118.946 ± 0.144 | 58.413 ± 0.000 | 0.491 ± 0.001 | +50.9% | -0.41045 | 0.16830 |
| SiO2 beta-quartz | 1382.801 ± 0.192 | 59.610 ± 0.003 | 0.043 ± 0.000 | +95.7% | -2.48975 | 0.10201 |

Each model has one train row and one test row per phase and seed in `per_phase_seed_metrics.csv`; the summary uses the common test seeds only.

The anchor is the final temperature in the frozen temp_extrap v2 training window for each phase. All curves are converted to anchored shape relative to that temperature before the primary RMSE is computed.

QH is quantum phonopy while the canonical DaRUS/TI reference is classical MD. Therefore the anchored thermal-shape error, not raw absolute-G MAE, is the primary comparison.

## Interpretation

- **Hf hcp:** QH reduces the anchored thermal residual (R = 0.067); the three seeds are consistent (R std = 0.0001). The improvement is mainly slope only (mean absolute offset: 8526.3 → 14170.4 meV; slope magnitude: 1.53214 → 0.10651 meV/atom/K).
- **Ti hcp:** QH reduces the anchored thermal residual (R = 0.245); the three seeds are consistent (R std = 0.0001). The improvement is mainly slope only (mean absolute offset: 7524.1 → 13285.1 meV; slope magnitude: 0.69241 → 0.15826 meV/atom/K).
- **Zr hcp:** QH reduces the anchored thermal residual (R = 0.491); the three seeds are consistent (R std = 0.0006). The improvement is mainly slope only (mean absolute offset: 8094.4 → 14981.5 meV; slope magnitude: 0.41045 → 0.16830 meV/atom/K).
- **SiO2 beta-quartz:** QH reduces the anchored thermal residual (R = 0.043); the three seeds are consistent (R std = 0.0000). The improvement is mainly offset and slope (mean absolute offset: 1830.3 → 1149.8 meV; slope magnitude: 2.48975 → 0.10201 meV/atom/K).

Across all four QH-reliable phases, QH gives R < 1, but this is a single-phase thermal-shape result only. It does not establish pairwise phase-stability improvement, does not validate half-QH on a dynamically unstable high-temperature partner, and does not identify a true entropy error.

## Reproduce

```bash
python src/experiments/single_side_qh.py --repo-root .
```

The script verifies QH reliability, the frozen split hash, exact temperature grids, common seed sets, eV/atom input units, and zero anchored residual at T_ref. Prediction CSVs are the canonical T3 tlog and T2 QH-residual outputs exported from the isolated runs; their source paths, checkpoint placeholders, heads, hashes, and commits are recorded in `per_phase_seed_metrics.csv` and `findings.meta.json`.

## Extensions: controls and pair diagnostic

This extension keeps the original temp_extrap v2, seed 11/23/37 convention and writes all new controls under this directory. No model was retrained.

### 1. Training-only thermo-form upper bounds

For each QH-reliable phase, `thermo_abs_G` fits `G_ref(T)` and `thermo_QH_residual` fits `G_ref(T) − F_QH(T)` with `a + bT + cT ln(T)` on train points only. The reported error is the same anchored test RMSE used above.

| Phase | thermo abs-G RMSE (meV/atom) | thermo QH-residual RMSE (meV/atom) | R = QH/abs | learned T3 abs-G | learned T2 QH |
|---|---:|---:|---:|---:|---:|
| Hf hcp | 1.841 | 1.835 | 0.997 | 421.351 | 28.069 |
| Ti hcp | 0.991 | 0.947 | 0.956 | 214.456 | 52.499 |
| Zr hcp | 2.210 | 2.144 | 0.970 | 118.946 | 58.413 |
| SiO2 beta-quartz | 1.496 | 0.016 | 0.011 | 1382.801 | 59.610 |

These controls separate the representation model from the physics baseline: the thermo-form rows are non-learning ceilings for the same basis and train/test window. Here the learned T3 errors are much larger than the thermo-form ceiling, so its dominant limitation is the learned absolute-G representation/extrapolation. T2 is closer than T3 but still above the QH thermo ceiling, especially for SiO2; that residual gap is not explained by the smooth baseline alone.

### 2. Einstein smooth baseline

`einstein_parameters.csv` fits the single positive parameter θ_E to the anchored training shape of F_QH. `thermo_Einstein_residual` fits `G_ref − F_E` with the same three-term thermo form. `learned_Einstein_transplant` is explicitly a no-retraining diagnostic: the existing T2 learned residual `(T2 − F_QH)` is added to F_E. It is not presented as a newly trained model.

| Phase | theta_E (K) | thermo Einstein RMSE (meV/atom) | learned Einstein transplant (meV/atom) | QH learned RMSE (meV/atom) |
|---|---:|---:|---:|---:|
| Hf hcp | 130.172 | 1.836 | 28.067 | 28.069 |
| Ti hcp | 227.847 | 0.955 | 52.478 | 52.499 |
| Zr hcp | 162.522 | 2.157 | 58.391 | 58.413 |
| SiO2 beta-quartz | 577.863 | 0.896 | 58.112 | 59.610 |

The Einstein comparison is a shape control, not a claim that one oscillator is a physical phonon spectrum. Similar QH and Einstein errors would support a smooth-shape interpretation; a large difference would indicate phase-specific QH information matters.

### 3. QH-invalid negative controls

The invalid phases use the same non-learning comparison after the adopted imaginary-mode exclusion. They are not used in the valid-side conclusions.

| Phase | QH-reliable | thermo abs-G RMSE | thermo QH-residual RMSE | R = QH/abs |
|---|---|---:|---:|---:|
| Hf bcc | false | 1.216 | 1.213 | 0.997 |
| Ti bcc | false | 0.440 | 2.593 | 5.893 |
| Zr bcc | false | 1.328 | 1.287 | 0.969 |
| SiO2 beta-cristobalite | false | 1.309 | 0.205 | 0.156 |
| SiO2 beta-tridymite | false | 1.439 | 0.067 | 0.046 |

The expected R ≥ 1 is a diagnostic expectation, not a forced acceptance criterion; the measured values are reported as-is.

### 4. Pair-level mixed predictor

`single_side_QH_combo` uses the learned T2 QH curve for the low-temperature/QH-reliable phase and the learned T3 absolute-G curve for the partner. It is compared with `T3_absolute_both` and the existing canonical `thermo_form_fit`. Pair metrics retain the frozen crossing labels, including N/A reasons.

| Predictor | Pair rows | ΔG MAE (meV/atom) | sign accuracy | mean |Tc error| (K) |
|---|---:|---:|---:|---:|
| T3_absolute_both | 18 | 31.764 ± 34.762 | 0.4845 | 36.405 |
| single_side_QH_combo | 15 | 644.868 ± 465.591 | 0.8615 | 100.407 |
| thermo_form_fit | 6 | 0.425 ± 0.394 | 0.9948 | 2.796 |

For the mixed pair, the T2 low-side curve is shifted to the T3 low-side value at the frozen train anchor; the shift is recorded as `gauge_alignment_meV` and uses no test labels. The third SiO2 pair has no QH-reliable low side and is explicitly N/A. The pair diagnostic is the direct test of whether single-sided QH helps a phase transition. It does not change the main benchmark tables or imply that QH is valid on the high-temperature unstable side.

## Reproduce extensions

```bash
python src/experiments/single_side_qh_extensions.py --repo-root .
```

The script validates exact frozen train/test grids, eV/atom prediction files, split hashes, QH source files, common seeds, and the no-retraining provenance of all controls.

## DFT static replacement and constant-sign control

### DFT static replacement

The canonical QH files contain the static term: for SiO2, `F_QH − F_vib` is constant; for metal QH files the row-wise `E_static_eV_per_atom` column is used because the minimized-volume static energy can vary with temperature. The replacement therefore uses `G_corr(T) = F_QH(T) − E_static,DPA(T) + E_DFT = F_vib(T) + E_DFT`, with both phases of each pair treated identically. `E_DFT` is the relaxed static `energy without entropy` from the `final-adopted-converged` OUTCAR manifest. The DFT output status is recorded per phase; the adopted production outputs are present, while the Hf neighboring-mesh deviation remains documented in the DFT convergence report.

| Split | Predictor | ΔG MAE (meV/atom) | ΔG RMSE (meV/atom) | sign accuracy | balanced sign accuracy | mean |Tc error| (K) | false | missed |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| phase_lopo | DFT_static_replaced | 122.352 | 124.594 | 0.5293 | 0.5000 | n/a | 0 | 10 |
| phase_lopo | T1_E_plus_F_QH | 167.628 | 169.316 | 0.5293 | 0.5000 | n/a | 0 | 10 |
| system_loso | DFT_static_replaced | 246.990 | 251.775 | 0.4844 | 0.5000 | n/a | 0 | 3 |
| system_loso | T1_E_plus_F_QH | 333.575 | 337.181 | 0.4844 | 0.5000 | n/a | 0 | 3 |
| system_loso_overlap_T | DFT_static_replaced | 290.342 | 291.527 | 1.0000 | 1.0000 | n/a | 0 | 1 |
| system_loso_overlap_T | T1_E_plus_F_QH | 412.096 | 412.932 | 1.0000 | 1.0000 | n/a | 0 | 1 |
| temp_extrap | DFT_static_replaced | 119.614 | 120.278 | 0.5589 | 0.5000 | n/a | 0 | 5 |
| temp_extrap | T1_E_plus_F_QH | 160.504 | 161.011 | 0.5589 | 0.5000 | n/a | 0 | 5 |

### Constant-sign baseline

`constant_sign` predicts the oriented pair ΔG as −1×10⁻⁹ eV/atom at every point, i.e. the high-temperature/right phase is always declared stable. The tiny nonzero magnitude avoids the exact-zero degenerate-prediction marker; this is a sign-only control, not an energy baseline. Balanced accuracy is the mean recall of the positive and negative reference sign classes present on the evaluated grid; for a one-sign pair the present class is reported.

| Split | ΔG MAE (meV/atom) | ordinary sign accuracy | balanced sign accuracy |
|---|---:|---:|---:|
| phase_lopo | 4.706 | 0.5711 | 0.6057 |
| system_loso | 9.047 | 0.4844 | 0.5000 |
| system_loso_overlap_T | 21.931 | 1.0000 | 1.0000 |
| temp_extrap | 3.614 | 0.3224 | 0.6070 |

Pair-level rows, including fold scope and crossing markers, are in `dft_static_replacement_pair_metrics.csv` and `constant_sign_pair_metrics.csv`. Reproduce with:

```bash
PYTHONPATH=src python src/experiments/single_side_qh_dft.py --repo-root .
```
