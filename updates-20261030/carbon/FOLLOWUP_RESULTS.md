# Carbon entropy-transfer mechanism follow-up

## Frozen scope

This is a development-only post-analysis on the same 14 graphite/diamond fixed-volume `Svib(300 K)` labels. The follow-up design was written before the new baseline, decomposition, and calibration scoring in `FOLLOWUP_PROTOCOL.json` (SHA-256 `e860b01835f038130ce859fe18ce6ac60db46a24c604f9720e7076c5c15ec73a`). No additional calculations or model queries were run.

## A. Same-phase nonlinear volume baselines

MAE is in µeV/(atom·K); each number pools seven leave-one-volume-out predictions within that phase. `V/V0` was fit-only standardized where relevant. The quadratic model is fixed OLS on `[1, V/V0, (V/V0)^2]`; RBF kernel ridge uses nested source-only CV; nearest-volume and mean are fixed baselines.

| Phase | DPA (parent run) | Linear volume (parent) | Quadratic `V/V0` | RBF KRR | Nearest volume | Training mean |
|---|---:|---:|---:|---:|---:|---:|
| Graphite | 0.00743 | 0.0893 | **0.00467** | 0.1925 | 0.6014 | 1.1841 |
| Diamond | 0.00920 | 0.1129 | **0.00372** | 0.1946 | 0.9265 | 1.9042 |

The quadratic volume-only baseline beats DPA by 37% in graphite and 60% in diamond on these small within-phase folds. At interior points its MAE is 0.00335 (graphite) and 0.00211 (diamond); at the two endpoint extrapolation points it is 0.00799 and 0.00775. DPA's corresponding interior/endpoint MAEs are 0.00478/0.01406 and 0.00419/0.02173. Thus the earlier “DPA is about 12× better than linear volume” comparison was baseline-dependent: a simple quadratic volume trend explains the same-phase curves at least as well. The result does not show that DPA contributes information beyond volume for this task.

## B. Cross-phase error shape

The target phases use the same `V/V0` interval, but their absolute volumes do not overlap: graphite 7.544–8.507 and diamond 5.229–5.897 Å³/atom, with a 1.647 Å³/atom interval gap. Source-standardized DPA descriptor distances are also very large: target-centroid RMS z-distance is 84.0 for graphite→diamond and 39.0 for diamond→graphite; mean nearest-source RMS z-distances are 84.1 and 39.1. With only seven source structures and 128 descriptor dimensions, these are descriptive OOD flags, not calibrated probabilities.

| Transfer / model | Reference mean | Predicted mean | Signed bias | Reference slope | Predicted slope | Slope error |
|---|---:|---:|---:|---:|---:|---:|
| Graphite→diamond, DPA | 24.64 | 39.91 | +15.27 | +47.19 | +25.14 | −22.06 |
| Graphite→diamond, DPA+V | 24.64 | 39.94 | +15.30 | +47.19 | +24.95 | −22.24 |
| Diamond→graphite, DPA | 50.62 | 28.09 | −22.53 | +29.28 | −8.76 | −38.04 |
| Diamond→graphite, DPA+V | 50.62 | 28.27 | −22.36 | +29.28 | −8.14 | −37.42 |

Slopes are linear fits of entropy against target-phase `V/V0`, in µeV/(atom·K) per unit `V/V0`; they summarize only these seven points. Cross-phase error is not just a constant anchor shift: the DPA slope is too shallow in graphite→diamond and has the wrong sign in diamond→graphite. DPA+V barely changes either transfer curve under the primary coordinate-standardized rule.

## C. Explicit target-label calibration diagnostic

This diagnostic deliberately uses target entropy labels and is **not zero-shot transfer**. The one-point correction uses target `v1000` to shift the curve by a constant, then scores the other six volumes. The two-point correction uses `v0980` and `v1020` to fit an affine residual in `V/V0`, then scores the remaining five.

| Transfer / model | Calibration | MAE before | MAE after | Scored target volumes |
|---|---|---:|---:|---:|
| Graphite→diamond, DPA | one point | 15.25 | 0.888 | 6 |
| Graphite→diamond, DPA | two points | 15.23 | 0.167 | 5 |
| Diamond→graphite, DPA | one point | 22.52 | 1.522 | 6 |
| Diamond→graphite, DPA | two points | 22.52 | 0.032 | 5 |

A target-phase anchor removes most of the large mean bias, while a second fixed anchor corrects much of the remaining volume-response error. This supports an anchor-plus-response decomposition as a useful diagnostic for this dataset, but the tiny correlated sample and predetermined calibration points do not establish a generally sufficient two-label calibration strategy.

## Reference precision

The labels exactly reproduce the existing independent `SVIB_CV_CANDIDATE.csv` values at 300 K to floating-point precision. The source workflow checks 16³ q-mesh spectra and 12³/16³/20³ q-mesh effects on `Fvib`, but the archived convergence report does not provide an uncertainty estimate specifically for `Svib(300 K)`. Therefore the 0.003–0.009 µeV/(atom·K) model MAEs cannot be interpreted as accuracy relative to the physical or numerical uncertainty of DFT phonon entropy. They measure agreement with this saved 16³ harmonic label set.

## Interpretation

For these two phases, same-phase `Svib(V)` is captured very well by a quadratic scalar-volume curve; the DPA linear readout is not the strongest baseline once that control is included. Whole-phase transfer has large phase-level bias and also incorrect volume response, with both volume range shift and strong descriptor-space shift. A small amount of target entropy calibration largely repairs the prediction, but that is a distinct supervised task. These observations motivate broader tests; they do not establish a general structural-manifold law from one system.

## Artifacts

- `VOLUME_BASELINE_METRICS.csv`, `VOLUME_BASELINE_PREDICTIONS.csv`, `VOLUME_BASELINE_SELECTION.csv`
- `CROSS_PHASE_DECOMPOSITION.csv`
- `TARGET_CALIBRATION_DIAGNOSTIC.csv`
- `run_carbon_entropy_mechanism_followup.py`
