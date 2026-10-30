# Carbon fixed-volume harmonic entropy transfer: results

## Scope and status

This is a development analysis of the existing independent PBE-D3(BJ), quantum-harmonic fixed-volume reference at 300 K. It contains one carbon system and two phases (graphite and diamond), with seven nested volume structures per phase. It is not 14 independent materials and is not a blind confirmation set. The primary target is `Svib(300 K)`; QHA `G(T)` data were not used for model selection or phase-boundary claims.

Frozen DPA-3.1-3M / `Domains_Drug` mean-pooled 128D representations were extracted on CPU without reading thermodynamic labels. The selected branch contains C. Checkpoint SHA-256: `86dd3a804d78ca5d203ebf98747e8f16dff9713ba8950097ceb760b161e19907`. Descriptor NPZ SHA-256: `021ea40901c288243f25e456b5f3c6983a510a05575de548762a84a41dde5e50`. The run protocol hash is `e211b4ee9aa2586aba57a1b06e64b1b34c234e1174b54f88ce30d57dc0f9ca07` and was frozen before descriptor extraction and target scoring.

## Main results

MAE is in µeV/(atom·K). The primary scaling is coordinate-wise standardization, matching the predeclared SiO₂ A/B/C setting. Within-phase values are leave-one-volume-out across all seven held volumes. Endpoint extrapolation and five-volume interpolation are reported separately in `METRICS.csv`.

| Evaluation | Arm A: volume | Arm B: DPA | Arm C: DPA + volume |
|---|---:|---:|---:|
| Graphite, within-phase LOOV MAE | 0.0893 | **0.00743** | 0.00748 |
| Diamond, within-phase LOOV MAE | 0.1129 | **0.00920** | 0.00917 |
| Graphite→diamond whole-phase MAE | 17.01 | **15.27** | 15.30 |
| Diamond→graphite whole-phase MAE | **5.12** | 22.53 | 22.36 |

The within-phase experiment is strongly favorable to DPA over the scalar-volume linear baseline, but adding volume to DPA barely changes the primary results: the DPA-conditioned volume error reduction `E(DPA)-E(DPA+V)` is −0.000048 µeV/(atom·K) for graphite and +0.000034 for diamond. The tiny differences are descriptive and do not establish useful incremental volume information.

Whole-phase transfer is poor and asymmetric. With graphite as the sole source, all three arms overpredict diamond entropy by about 15–17 µeV/(atom·K). With diamond as the sole source, DPA-based arms underpredict graphite by about 22.4–22.5 µeV/(atom·K); volume alone has a smaller 5.12 µeV/(atom·K) MAE in that direction. These are seven correlated volume states in each direction, not independent phase or material replicates. The result is consistent with substantial phase-level entropy offsets that cannot be calibrated from a single source phase.

## Scaling sensitivity and attribution

Under the predeclared equal-total-group-weight sensitivity, DPA+volume MAE is 0.01073 for graphite and 0.00729 for diamond, compared with DPA-only 0.00743 and 0.00920. Thus weighting changes the direction of the combined-arm difference by phase. In whole-phase transfer, equal weighting changes some outcomes substantially (for diamond→graphite, C improves to 12.32 while B remains 22.53), illustrating that the readout's regularization geometry matters. This is a sensitivity result, not evidence that one scaling is preferred; the primary remains coordinate-standardized.

Within-phase, DPA’s improvement over volume is large. Yet these are closely related structures within each phase, and DPA pretraining structure exposure is undocumented. This supports within-phase volume-series prediction in this benchmark, not structure-zero-shot or material-level transfer. DPA+V versus DPA differences are the predeclared paired contrasts in `INCREMENTS.csv`; positive values mean the first arm's MAE minus combined-arm MAE, so positive favors DPA+V for `V_given_DPA`.

## Fitting, computation, and limits

The run produced 96 outer fits and 168 held-out predictions across A/B/C, two scaling rules, 14 within-phase folds, and two directed whole-phase folds. Alpha was chosen by nested leave-one-volume-out MAE within source data only; all standardization was fit-only. No target label was used in scaling, feature weighting, or alpha selection. The frozen feature extraction took 11.47 CPU seconds on the configured remote runtime. No DFT, MD, VASP, or Bohrium job was submitted.

Reference labels are the fixed-volume harmonic vibrational entropy, not the full finite-temperature entropy including anharmonic/electronic terms. The four-temperature QHA `G` table has only 16 points and was not used here; no Tc or phase-event score is appropriate in this experiment. Carbon results remain a development probe because this reference set was previously inspected in the benchmark/website workflow.

## Files

- `PROTOCOL.json`, `FOLDS.csv`, `SELECTION_FREEZE.json`: frozen design and folds.
- `PREDICTIONS.csv`: every held-out prediction with reference and error.
- `METRICS.csv`, `INCREMENTS.csv`: primary and sensitivity metrics/paired contrasts.
- `FIT_AND_SELECTION.csv`, `INNER_SELECTION_SCORES.csv`: source fit and nested model-selection records.
- `MODEL_PROVENANCE.json`, `EXTRACTION_PROVENANCE.json`, `INPUT_STRUCTURE_MAP.csv`: checkpoint and input provenance.
- `run_carbon_entropy_transfer.py`: local fitting/scoring implementation.

## Follow-up mechanism analysis

A predeclared follow-up tested nonlinear volume-only controls and decomposed the whole-phase errors. The quadratic `V/V0` baseline outperformed DPA within both phases (MAE 0.00467 for graphite and 0.00372 for diamond), so the original advantage over linear volume does not establish DPA-specific incremental information. Cross-phase DPA predictions also have both large mean bias and incorrect `Svib`-versus-`V/V0` slopes. A fixed one-point target-label offset correction removes most of the mean error, and a fixed two-point affine residual correction further reduces it; these are supervised calibration diagnostics, not zero-shot scores. See [FOLLOWUP_RESULTS.md](FOLLOWUP_RESULTS.md) and its CSV products for details.

The existing convergence materials provide q-mesh sensitivity for `Fvib`, but not a direct uncertainty estimate for `Svib(300 K)`. The very small within-phase model errors therefore describe agreement with the saved 16³ harmonic labels, not demonstrated accuracy below DFT phonon uncertainty.
