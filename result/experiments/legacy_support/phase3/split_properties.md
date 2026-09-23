# Frozen-split properties relevant to SiO2 phase ordering

## `temp_extrap`: cristobalite_beta − quartz_beta

The frozen SiO2 temperature split uses train indices 0–591 (851–1442 K) and test indices 592–1648 (1443–2499 K).  For the reference difference ΔG = G(cristobalite) − G(quartz):

| region | points | positive | negative | zero | range (meV/atom) |
|---|---:|---:|---:|---:|---:|
| train | 592 | 592 | 0 | 0 | +0.756 to +4.884 |
| test | 1,057 | 100 | 957 | 0 | −7.494 to +0.748 |

Thus the training set contains only the positive ordering.  The reference crossing is 1542.71 K, so 957/1,057 test points (90.5%) have the opposite, negative ordering; the first 100 test points remain positive because the frozen threshold is 100 K below the crossing.  This is the intended behavior of the temperature-extrapolation split: the sign reversal and its crossing are held out, not evidence by itself that a predictor has failed.

Paper-method wording: “For `temp_extrap`, the cutoff is selected so every reference crossing lies in the test region with a 100 K buffer below it. Consequently, a pair may have one sign only in training and the opposite sign predominantly in test; phase-order accuracy must be interpreted as extrapolation across an unseen sign reversal.”

## Deviations from design

No deviation. The counts were computed directly from the frozen `splits/temp_extrap.json` and the canonical SiO2 reference table; no model was trained or evaluated to establish this property.
