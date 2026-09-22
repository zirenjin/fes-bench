# E1–E14 status — 2026-09-22 overnight queue

| item | status | evidence | blocker / note |
|---|---|---|---|
| E1 | complete | `results/phase5/e1_tc_table.md`, `sio2_reanalysis.json` | frozen-checkpoint re-evaluation persisted |
| E2 | complete | `results/phase5/e2_leakage_audit.md` | train/heldout calibration audit persisted |
| E3 | complete | `results/phase5/e3_recovery.json` | synthetic recovery accepted |
| E4 | complete | `results/trivial_floor/metrics.json`, `skill_scores.md` | zero and pairwise constant floors; skills all negative for E1 |
| E5 | complete (with structural N/A folds) | `results/external_baselines/metrics.json`, `skill_scores.md`, `phase_id_mlp_matagent.md` | Bartel/interp/global floors run; mat-agent phase-ID MLP gives temp-extrap skill −0.858; LOPO/LOSO N/A because held-out targets have no training rows |
| E6 | complete (raw QH) | `results/phase2/sio2_qh_report.md`, remote `results/phase2_sio2_sse_pbe/sio2/qh_summary.json` | three-phase SSE-PBE QH CSVs and comparison are complete; all minimum frequencies are negative |
| E7 | complete | `results/phase2/sio2_qh_report.md`, isolated `results/phase2_sio2_sse_pbe/*_imaginary_diagnosis.json` | all three phases have non-Γ-local negative modes and `qh_reliable=false` |
| E8 | complete | `results/phase6_readiness/compile.log`, `pytest_vendor_full.log` | compile passed; vendored full suite 21/21 passed in 4750.31 s |
| E9 | complete (preflight) | `results/phase5/sio2_second_truth_source_survey.md` | source/format surveyed, no download per instruction |
| E10 | complete | `results/phase3/split_integrity.json` | frozen split hashes verified |
| E11 | complete for active domain | `data/{sio2,hf}/`, `results/phase1_delivery.md` | CaSiO3 remains source hold |
| E12 | partial | `results/phase6_readiness/readiness.md` | tests are green; Phase 6 readiness still awaits QH and V100 environment sync |
| E13 | complete | `results/table2_delta_g_amplitude.md` | pair ΔG amplitude table generated |
| E14 | partial | this report and `results/overnight_20260922.md` | QH and reliability diagnosis finished; V100 environment sync remains blocked |

## Deviations from design

No new benchmark model was trained, no frozen split was changed, and no head policy was changed. Missing energies, broken torch runtime, and absent pytest are recorded as explicit holds rather than substituted.
