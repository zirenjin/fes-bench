# E1–E14 status — 2026-09-22 overnight queue

| item | status | evidence | blocker / note |
|---|---|---|---|
| E1 | complete | `results/phase5/e1_tc_table.md`, `sio2_reanalysis.json` | frozen-checkpoint re-evaluation persisted |
| E2 | complete | `results/phase5/e2_leakage_audit.md` | train/heldout calibration audit persisted |
| E3 | complete | `results/phase5/e3_recovery.json` | synthetic recovery accepted |
| E4 | complete | `results/trivial_floor/metrics.json`, `skill_scores.md` | zero and pairwise constant floors; skills all negative for E1 |
| E5 | partial | `results/external_baselines/metrics.json`, `skill_scores.md` | Bartel/interp/global floors run; phase-ID MLP unavailable because remote torch cannot load `libtorch_global_deps.so` |
| E6 | running | remote `/share/jzr/codex_fes_overnight_20260922_1790052610710037915/results/phase2_sio2_sse_pbe.log` | three-phase SSE-PBE QH process is live; final curves pending |
| E7 | pending | `configs/qh/sio2_qh_domains_sse_pbe.yaml` | report/plot follows E6 completion |
| E8 | partial | `results/phase6_readiness/compile.log`, `pytest.log` | compile passed; pytest package absent in all checked runtimes |
| E9 | complete (preflight) | `results/phase5/sio2_second_truth_source_survey.md` | source/format surveyed, no download per instruction |
| E10 | complete | `results/phase3/split_integrity.json` | frozen split hashes verified |
| E11 | complete for active domain | `data/{sio2,hf}/`, `results/phase1_delivery.md` | CaSiO3 remains source hold |
| E12 | partial | `results/phase6_readiness/readiness.md` | Phase 6 readiness awaits QH and pytest evidence |
| E13 | complete | `results/table2_delta_g_amplitude.md` | pair ΔG amplitude table generated |
| E14 | partial | this report and `results/overnight_20260922.md` | overnight queue still has live QH work |

## Deviations from design

No new benchmark model was trained, no frozen split was changed, and no head policy was changed. Missing energies, broken torch runtime, and absent pytest are recorded as explicit holds rather than substituted.
