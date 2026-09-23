# SiO₂ β-quartz head check — stop decision

Run: isolated thu-GenSi directory
`/GenSIvePFS/users/zirenj/fes-bench/runs/phase2_sio2_quartz_head_20260922`.
Checkpoint: `DPA-3.1-3M.pt`; supercell `2×2×2`; Gamma-centered mesh
`[12,12,12]`; displacement `0.01 Å`; ASR and force-constant symmetry both
disabled.  Negative threshold is `−0.05 THz`, with the Hf diagnostic's 5%
conservative-BZ-radius Γ neighborhood.

| head | single-point E (eV/atom) | max force (eV/Å) | min frequency (THz) | negative modes | negative q points | Γ-neighborhood negative modes |
|---|---:|---:|---:|---:|---:|---:|
| `Domains_Alloy` | −2.66506732 | 0.229281 | −1.926950 | 13,577 / 632,772 (2.1456%) | 868 / 868 (100%) | 20 / 13,577 (0.1473%) |
| `Domains_SSE_PBE` | −4.08071262 | 0.099727 | −0.975320 | 235 / 632,772 (0.03714%) | 183 / 868 (21.08%) | 0 / 235 |

## Decision

`Domains_Alloy` fails the Hf-derived head criterion: negative modes occur at
every sampled q point and are overwhelmingly outside the Γ neighborhood.
This is an all-zone out-of-domain result, not a Γ-local acoustic/ASR artifact.

Per the objective, stop here.  Do **not** switch the SiO₂ policy to
`Domains_Alloy`, do **not** rerun the three SiO₂ relaxations or QH, and do not
roll back the existing policy.  The existing `Domains_SSE_PBE` row is retained
only as the requested same-structure control; it also has nonzero dispersed
negative q points, so it is not being promoted as a new physical baseline.

Raw evidence: `head_diagnosis.json` and `sio2_alloy_head_check.md` under the
remote run directory above; the standalone driver is
`scripts/diagnose_sio2_quartz_remote.py`.

## Deviations from design

- None in the decision rule: the authorized `−0.05 THz` threshold, 5% Γ
  neighborhood, 2×2×2 displacement setup, and no-ASR/no-fc-symmetry protocol
  were used.  P0-2 and all later items were intentionally not entered because
  P0-1 failed.
