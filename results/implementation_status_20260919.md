# G(T) benchmark implementation-status audit — 2026-09-21

This is a status audit, not a completion claim.  Each row is tied to a
persisted artifact or an explicit missing-input hold.

| Requirement area | Current evidence | Status |
|---|---|---|
| Independent repository and canonical loader | `fes_bench.data.load`, schema tests, and `README.md` | Implemented |
| SiO2/Hf source tables and crossing audit | `data/{sio2,hf}/reference_crossings.json`; `results/phase1_delivery.md` | Implemented for active domain |
| Ti/Zr representatives | `results/phase1_delivery.md`, `results/phase1_thu_audit_20260920.log` | Held: structures absent; fresh thu audit preserves the hold |
| CaSiO3 absolute `G(T,P)` tables | `results/phase1/casio3_source_hold.md` | Held: only relative/source-incomplete data |
| QH production baseline | `results/phase2_delivery.md`; `data/hf/*/{fqh.csv,phonon_report.json}`; `results/phase2_materialization_20260921_v100_r17.log` | Held: Hf diagnostic curves are canonicalized for provenance but both reports set `qh_reliable:false`; V100 r17 full regression passed; no invalid FQH promoted |
| Frozen temp/LOPO/LOSO splits | `splits/*.json`, `results/phase3/split_integrity.json` | Implemented; embedded hashes pass |
| Reference and noise evaluator acceptance | `tests/test_eval_roundtrip.py`; V100 r12–r14 and r18 logs; `results/phase3_thu_hf_revalidate_20260920.log` | Implemented; exact reference passthrough, 5 meV noise run, and crossing-local MAE metric are verified |
| ΔG curves and metric tables | `results/reference/*/{metrics.json,summary.md,delta_g_curves.*}` | Implemented for all three splits |
| F_QH physical-baseline head and fixed ablations | authorized `feat/fes-head` diff; `configs/models/`; `results/phase4_delivery.md`; `results/phase4/v100_upstream_fes_20260921.log` | Implemented; scoped upstream FES file passes 33 tests on V100; whole-repository pytest remains out of scope |
| `c_system` calibration primitive | `fes_bench/models/calibrate.py`, `tests/test_calibrate.py` | Implemented; not used in a production prediction absent valid FQH data |
| Bartel, interpolation, phase-ID MLP fixtures | `fes_bench/baselines/`, Phase 5 tests; V100 r15 | Implemented; no production benchmark fit run |
| E3 recovery | `results/phase5/e3_recovery.json`; V100 r15 recheck | Implemented; qh-residual continuous-tlog vs repr-only MLP, MAE ratio=0.8203 |
| E1/E2 existing-checkpoint re-evaluation | `results/phase5/e1_e2_reanalysis_hold.md` | Held: Hf checkpoint absent; SiO2 legacy four-phase domain differs from canonical three-phase split |

## V100 evidence

`results/phase3/v100_r12_r14.log` records the original decisive isolated runs;
`results/phase3/v100_r18_noise_20260921.log` records the fresh noise and
crossing-local-MAE revalidation.  The
latest Phase 3 verifier run (r14) returned `14 passed in 5.00s`.  The
dedicated thu Hf-only revalidation regenerated all three split kinds and
returned exact-zero reference metrics.  Phase 5 r11 returned `11 passed in
6.81s` for baseline fixtures and E3 recovery.

## Remaining conditions before a complete benchmark claim

1. Supply/source-audit Ti/Zr representative structures and CaSiO3 absolute
   phase-resolved `G(T,P)` data, or retain those systems outside the minimal
   benchmark.
2. Produce QH curves that pass the stated reliability gate, or label the
   affected ablations as unavailable rather than physical baselines.
3. Provide matching frozen checkpoints and frame mappings for the canonical
   SiO2/Hf E1/E2 re-evaluation.
4. Run any desired deepmd-kit tests outside the scoped
   `source/tests/pt/test_fes.py` file; that complete FES file now passes on
   V100, while focused eager/JIT/DeepEval tests also pass on thu-GenSi.
5. Establish the independent repository's initial Git commit if commit-based
   provenance is required in addition to the already verified split hashes.

## Deviations from design

The active benchmark deliberately remains narrower than the paper plan until
the missing source inputs and physical-QH gate are resolved.  No missing
structure, FQH curve, absolute free energy, checkpoint result, or commit hash
is substituted with a surrogate.
