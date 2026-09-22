# Phase 6 readiness gate

| gate | result | evidence |
|---|---|---|
| Environment sync | pass (thu-GenSi locked) | `results/phase6_readiness/env_lock.md`; isolated thu-GenSi `/share/jzr/codex_fes_overnight_20260922_1790052610710037915`, conda export and 183-line pip freeze recorded. |
| Python syntax | pass | `compile.log`, `compile_exit:0` |
| Full Phase 3/4/5 pytest | pass | isolated vendored run: 21 passed in 4750.31s; `pytest_vendor_full.log`, exit 0 |
| Frozen splits | pass | `results/phase3/split_integrity.json` |
| Trivial floors | pass | `results/trivial_floor/metrics.json` |
| External baselines | pass (structural N/A folds recorded) | `results/external_baselines/`, including `phase_id_mlp_matagent.md`; phase-ID MLP temp-extrap skill −0.858, LOPO/LOSO N/A by split construction; Hf DPA E0 helper remains explicitly unavailable |
| SiO₂ QH | **pending** (technical clarification recorded) | `results/phase2/sio2_qh_gate_audit.md` and `results/phase2/sio2_qh_report.md`; source-level cutoff behavior and actual `[1,1,1]` cell edges are recorded, while all phases remain `qh_reliable=false`. |
| Backup resource (V100) | backup — balance pending | Direct SSH probe timed out; authenticated Bohr query reports 0 CNY; not a readiness gate. |

Conclusion: **Phase 6 is not opened.** Readiness prerequisites are locked to thu-GenSi; the SiO₂ QH gate remains pending, and V100 is backup-only pending recharge. No benchmark-model training was started.

## Deviations from design

The first full pytest attempt was blocked by a missing package; a vendored pytest install was then made only inside the isolated remote workspace. The corrected full run passed all 21 tests. No frozen split, metric definition, or head policy was changed.
