# Deviations from PROTOCOL.md (SHA-256 0b382dba…)

## D1. Spectrum binning (recorded 2026-11-02, before any full-data route result)

**What the protocol said.** Bins cover 0.05–60 THz, and model modes ≤ 0.05 THz are excluded, "as phonopy does".

**What was found.**
- The frozen fidelity gate failed on the first 20 complete formulas, a smoke run used only to test the code. Re-integration deviated from the stored F(1000 K) by up to 1.7 meV/atom in about 11% of entries.
- In phonopy 3.5 the thermal-property cutoff defaults to **0 THz**: every mode with ν > 0 is included, with no renormalisation (`phonopy/phonon/thermal_properties.py`). The stored F therefore includes the 0–0.05 THz modes, and these are large in entries with imaginary branches.

**Change.**
- Re-integration bins: 480 log bins over 1e-4–60 THz, containing every model mode with ν > 0, which matches phonopy. Modes below 1e-4 THz go into the first bin.
- Route C's stage-1 ln-ratio target uses only band-matched modes where both the model and DFT frequencies exceed 0.05 THz, kept as a separate weight array. This stops near-zero acoustic modes from dominating the target.
- The fidelity gate itself (within 0.05 meV/atom at 1000 K for at least 99% of entries) is unchanged.

**Affected runs.** The spectrum run was stopped and restarted from scratch. The smoke-run outputs (`smoke/`) are discarded and not reported.

## D2. Gate reference and how route C is applied (recorded 2026-11-02, before any full-data route result)

**What was found.** After D1, binning reproduced the stored F within 0.002–0.03 meV/atom for F0, but deviated by up to 0.08 meV/atom for MatterSim, including on a stable entry (mp-22865). The tail audit had already shown the same thing: recomputing MatterSim F_th differs from the stored value by about 0.08 meV/atom, while DPA matches exactly. MatterSim GPU forces are not bitwise reproducible.

**Change.**
- The fidelity gate now compares re-integrated bins with the **same-run** phonopy F. This is a test of the binning, which is what the gate is for.
- Recompute reproducibility against the stored F is reported separately.
- Route C is applied as a difference, F_C = F_stored + [F_bins(corrected) − F_bins(uncorrected)]. Route C then differs from "none" only by its spectral correction, not by recompute noise or binning error.
