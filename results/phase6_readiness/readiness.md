# Phase 6 readiness gate

| gate | result | evidence |
|---|---|---|
| Environment sync | partial | isolated thu-GenSi `/share/jzr/codex_fes_overnight_20260922_1790052610710037915`; V100 probe timed out (`dcwq1547908.bohrium.tech:22`, retried 2026-09-22 07:06 UTC). Authenticated Bohr resource query lists V100 SKUs, but account balance is 0 CNY, so no paid replacement job was submitted. |
| Python syntax | pass | `compile.log`, `compile_exit:0` |
| Full Phase 3/4/5 pytest | pass | isolated vendored run: 21 passed in 4750.31s; `pytest_vendor_full.log`, exit 0 |
| Frozen splits | pass | `results/phase3/split_integrity.json` |
| Trivial floors | pass | `results/trivial_floor/metrics.json` |
| External baselines | partial | `results/external_baselines/`; mat-agent rerun is live for phase-ID MLP; Hf DPA E0 helper hit a PyTorch/e3nn compatibility load failure and was terminated |
| SiO₂ QH | running | remote PID 1996888, log under the isolated workspace |

Conclusion: **not yet Phase 6 ready**. The remaining blockers are concrete environment/runtime gates and the still-running QH calculation; no benchmark-model training was started.

## Deviations from design

The first full pytest attempt was blocked by a missing package; a vendored pytest install was then made only inside the isolated remote workspace. The corrected full run passed all 21 tests. No frozen split, metric definition, or head policy was changed.
