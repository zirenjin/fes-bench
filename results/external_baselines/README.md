# External-baseline results

Bartel coefficients are audited in `fes_bench/baselines/bartel2018.py`; phase-level baselines are unavailable for held-out phases without training rows and are marked explicitly. `global_mean_delta_g` is pair-only and uses all available training pair labels in the fold. The mat-agent phase-ID MLP rerun and its negative temp-extrap skill are summarized in `phase_id_mlp_matagent.md`.

## Deviations from design

No new model was trained. Bartel rows lacking a recorded representative DPA 0-K energy are reported unavailable rather than substituting a reference free energy. The phase-ID MLP is structurally unavailable on LOPO/LOSO held-out phases/systems because no training rows exist for those targets.
