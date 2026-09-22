# Phase 6 readiness gate

| gate | result | evidence |
|---|---|---|
| Environment sync | partial | isolated thu-GenSi `/share/jzr/codex_fes_overnight_20260922_1790052610710037915`; V100 probe timed out again (`dcwq1547908.bohrium.tech:22`, 2026-09-22 15:35 UTC). Authenticated Bohr resource query lists V100 SKUs, but account balance remains 0 CNY, so no paid replacement job was submitted. |
| Python syntax | pass | `compile.log`, `compile_exit:0` |
| Full Phase 3/4/5 pytest | pass | isolated vendored run: 21 passed in 4750.31s; `pytest_vendor_full.log`, exit 0 |
| Frozen splits | pass | `results/phase3/split_integrity.json` |
| Trivial floors | pass | `results/trivial_floor/metrics.json` |
| External baselines | pass (structural N/A folds recorded) | `results/external_baselines/`, including `phase_id_mlp_matagent.md`; phase-ID MLP temp-extrap skill −0.858, LOPO/LOSO N/A by split construction; Hf DPA E0 helper remains explicitly unavailable |
| SiO₂ QH | pass (with expected unreliability) | isolated `results/phase2_sio2_sse_pbe/sio2/` raw QH CSVs, `qh_comparison.json`, and three `*_imaginary_diagnosis.json` files; all phases `qh_reliable=false` |

Conclusion: **not yet Phase 6 ready**. All computational/data gates are green; the sole remaining gate is environment sync to V100, blocked by SSH timeout and zero Bohrium balance. No benchmark-model training was started.

## Deviations from design

The first full pytest attempt was blocked by a missing package; a vendored pytest install was then made only inside the isolated remote workspace. The corrected full run passed all 21 tests. No frozen split, metric definition, or head policy was changed.
