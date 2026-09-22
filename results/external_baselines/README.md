# External-baseline results

Bartel coefficients are audited in `fes_bench/baselines/bartel2018.py`; phase-level baselines are unavailable for held-out phases without training rows and are marked explicitly. `global_mean_delta_g` is pair-only and uses all available training pair labels in the fold.

## Deviations from design

No new model was trained. Bartel rows lacking a recorded representative DPA 0-K energy are reported unavailable rather than substituting a reference free energy.
