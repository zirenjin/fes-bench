# Benchmark status — 2026-10-02

This file records the E1–E14 inventory and the out-of-list overnight work. No frozen split or head policy was changed.

| item | status | evidence | blocker or note |
|---|---|---|---|
| E1 | complete | `results/phase5/e1_tc_table.md`, `sio2_reanalysis.json` | Frozen-checkpoint re-evaluation persisted. |
| E2 | complete | `results/phase5/e2_leakage_audit.md` | Train/held-out calibration audit persisted. |
| E3 | complete | `results/phase5/e3_recovery.json` | Synthetic recovery accepted. |
| E4 | complete | `results/trivial_floor/metrics.json`, `skill_scores.md` | Zero and pairwise constant floors are recorded. |
| E5 | complete with structural N/A folds | `results/external_baselines/metrics.json`, `skill_scores.md` | LOPO/LOSO have no eligible held-out training phase. |
| E6 | complete, raw QH | `results/phase2/sio2_qh_report.md`, remote QH summary | Three-phase SSE-PBE QH outputs are persisted. |
| E7 | pending technical gate | `results/phase2/sio2_qh_gate_audit.md`, imaginary-mode diagnostics | Non-Γ negative modes imply `qh_reliable=false`. |
| E8 | complete | `results/phase6_readiness/compile.log`, `pytest_vendor_full.log` | The recorded vendored suite passed. |
| E9 | complete | `results/phase5/sio2_second_truth_source_survey.md` | Source and format survey completed without an unrequested download. |
| E10 | complete | `results/phase3/split_integrity.json` | Frozen split hashes verified. |
| E11 | complete for the active domain | `data/{sio2,hf}/`, `results/phase1_delivery.md` | CaSiO3 remains on source hold. |
| E12 | partial | `results/phase6_readiness/readiness.md`, `env_lock.md` | thu-GenSi is locked; the SiO2 QH gate remains pending. |
| E13 | complete | `results/table2_delta_g_amplitude.md` | Pair ΔG amplitude table generated. |
| E14 | partial | `results/overnight_20260922.md`, QH gate audit | Raw QH and reliability diagnosis are complete; the technical gate remains pending. |

## Out-of-list overnight work

| line | status | evidence or blocker |
|---|---|---|
| A residual oracle and diagnostics | complete | `result/experiments/residual_oracle/`, thermo-form floor, README checks, and canonical pair checks are pushed. |
| B T3 temperature extrapolation | complete with gate failure recorded | All planned raw runs are present; only runs passing the train-region pair gate enter the main table. Polynomial is excluded; tlog is retained with its held-out metrics. |
| C T2 QH residual ablation | complete | Three-seed train/test evidence and controlled configuration are pushed. |
| D1 DPA-vs-reference-DFT snapshots | partial / completed where executable | Hf and SiO2 inference completed in an independent thu-GenSi directory; Ti/Zr remain `inference_unavailable` because the local runtime lacks deepmd and the independent remote directory lacked complete archives. |
| D2/D3 VASP static energies | partial, committed | Existing Expanse outputs were read only and assembled locally; all structure symmetries are audited, but final-adopted-converged output is absent. Pushed in `2afe4ae`. |
| E DeltaAI | not run | No accessible endpoint or credentials were available. |
| F integrity and reporting | complete | `src/notebooks/run_checks.py` exits with zero warnings; this status file and the dated overnight report are included in the final reporting commit. |

## Scope and safety

The authorized T3, T2, and VASP-static workspaces remain local; the T3 and DFT result commits are pushed, while the authorized temporary helper remains unstaged. No remote checkout received Git commands or workspace edits.
