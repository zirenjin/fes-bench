# Phase 6 readiness gate

| gate | result | evidence |
|---|---|---|
| Environment sync | partial | isolated thu-GenSi `/share/jzr/codex_fes_overnight_20260922_1790052610710037915`; V100 probe timed out (`dcwq1547908.bohrium.tech:22`, retried 2026-09-22 05:54 UTC) and relay could not resolve `bohrium-v100` |
| Python syntax | pass | `compile.log`, `compile_exit:0` |
| Full Phase 3/4/5 pytest | running | vendored pytest under isolated thu-GenSi workspace; PID 2148593, log `results/phase6_readiness_pytest_vendor_full.log` |
| Frozen splits | pass | `results/phase3/split_integrity.json` |
| Trivial floors | pass | `results/trivial_floor/metrics.json` |
| External baselines | partial | `results/external_baselines/`; mat-agent rerun is live for phase-ID MLP; Hf DPA E0 helper hit a PyTorch/e3nn compatibility load failure and was terminated |
| SiO₂ QH | running | remote PID 1996888, log under the isolated workspace |

Conclusion: **not yet Phase 6 ready**. The remaining blockers are concrete environment/runtime gates and the still-running QH calculation; no benchmark-model training was started.

## Deviations from design

The first full pytest attempt was blocked by a missing package; a vendored pytest install was then made only inside the isolated remote workspace and the corrected full run is now active. No frozen split, metric definition, or head policy was changed.
