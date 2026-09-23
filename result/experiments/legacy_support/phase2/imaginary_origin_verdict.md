# SiO2 beta-quartz imaginary-mode origin verdict

## Result

The beta-quartz sanity comparison identifies a substantial model/PES contribution to the earlier imaginary-mode result.  The calculation used the canonical fixed-cell `structure.extxyz`, unit-cell finite-displacement force constants (`supercell=[1,1,1]`), a `12x12x12` gamma-centered mesh, and `0.01 A` displacements.  No FES or TI data entered the calculation.

| model | checkpoint/head | force max (eV/A) | minimum (THz) | negative modes | negative-mode fraction | negative q-point fraction |
|---|---|---:|---:|---:|---:|---:|
| frozen pretrained | `DPA-3.1-3M.pt`, `Domains_Alloy` | 0.126303 | -7.72070 | 75,851 / 632,772 | 11.9871% | 100.0000% |
| SiO2 PES-finetuned | `model.ckpt-5000.pt`, head `None` | 0.000871 | -1.95616 | 1,718 / 632,772 | 0.271504% | 99.8848% |

Fine-tuning reduces the negative-mode count by 97.7 percentage points relative to the frozen model (about 97.7% fewer negative modes) and reduces the force residual by two orders of magnitude.  The remaining negative modes are sparse in mode count but occur at almost every q-point, so this is not a claim of a fully stable relaxed phase.  It is nevertheless a decisive model-sensitivity result: the frozen `Domains_Alloy` head's broad instability cannot be assigned to physical beta-quartz instability from this test.

## Verdict

**Model/PES issue; the earlier pretrained-head phonon conclusion must be redone with the SiO2-adapted PES (or an equivalently validated PES).**  This fulfills the P0-3(a) stop condition.  No Phase 6 work is started.

## Provenance and deviations

- Fine-tuned checkpoint: `/GenSIvePFS/users/zirenj/fes_experiment_pbe_d3bj/pes_dpa31_fes_full_singlehead_a100/model.ckpt-5000.pt`.
- The checkpoint was trained for 5,000 steps on the audited SiO2 PES dataset; its training config has `pref_v=0`, so no virial label was silently substituted.
- This is a fixed-cell harmonic diagnostic, not a production QH free-energy calculation.  Negative-mode threshold was `-0.05 THz`.
- Raw machine output is retained in `sio2_quartz_pes_comparison.json`.
