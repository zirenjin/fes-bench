# Deviations and corrections to PROTOCOL.md

## D1: test entry belongs to a group (recorded 2026-10-31, before any group-level MLIP result)

PROTOCOL.md states that neither pipeline-test entry belongs to a polymorph group. That is wrong for mp-10009 (GaTe), which forms the two-member group GaTe {mp-10009, mp-542812}. mp-1000 (BaTe) is in no group.

- **What was seen.** Single-phase numbers for mp-10009 only:
  - δZPE = −3.05 meV/atom;
  - δF_vib(600 K) ≈ −20 meV/atom;
  - absolute δE₀ = −63.6 meV/atom.
- **What was not seen.** Nothing for mp-542812, and no pair quantity.
- **Consequence.** No rule was changed. The GaTe group stays in the analysis, flagged, and every result is also reported with GaTe excluded.

## D2: compute platform and data transfer (operational, not analytic)

- The upload from the local machine to thu-GenSi was too slow (about 16 KB/s), so the files were downloaded on thu-GenSi directly from Alexandria.
- The 1,850 files are byte-identical on both machines (SHA-256 list in `data/RAW_SHA256SUMS.txt`).

## D3: added reference-consistency rule (recorded before any selection or MLIP result was seen)

`select_entries.py` adds one DFT-side rule not in PROTOCOL.md: the stored free energy must be reproduced from the stored force constants to ≤ 0.1 meV/atom. Entries that fail are listed, not repaired. This checks only the internal consistency of the reference and uses no MLIP information.

## D4: stability-threshold sensitivity (recorded after DFT-only selection statistics, before any MLIP result)

The frozen rule (minimum DFT mesh frequency ≥ −0.05 THz) removes 1,035 of the 1,850 group entries, although all of them are flagged stable in the Alexandria `pbe.csv`.

DFT-only distribution of the minimum mesh frequency (THz):

| Range | Entries |
|---|---:|
| ≥ −0.05 | 815 |
| [−0.1, −0.05) | 40 |
| [−0.2, −0.1) | 78 |
| [−0.5, −0.2) | 244 |
| [−1, −0.5) | 292 |
| < −1 | 381 |

Most exclusions are real imaginary branches on the dense mesh. For these, phonopy omits imaginary modes, so the harmonic F_vib is not defined. The dataset's `stable` flag is therefore more lenient than a full-mesh check.

- **Primary analysis.** Unchanged frozen rule: 239 groups, 521 entries, 351 pairs; 80 duplicates removed.
- **Added secondary sensitivity.** Threshold −0.2 THz, reported separately.

No MLIP quantity had been computed for any group entry when this was recorded. MLIP is run for all 1,850 entries, so that both selections use the same computations.
