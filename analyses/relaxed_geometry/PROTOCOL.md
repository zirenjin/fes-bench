# Relaxation channel: how much phase-relative error does the MLIP's own equilibrium geometry add?

Frozen before any MLIP relaxation is run. Zero new DFT; inference only (thu-GenSi GPU).

## Motivation

The F1 diagnosis found large MLIP forces on the undisplaced DFT supercells (0.009–0.049 eV/Å on four entries). The MLIP's own equilibrium geometry therefore differs from DFT's.

All earlier Alexandria tests evaluated the MLIP at the DFT geometry, so this channel was invisible to them. It is the static analogue of the relaxation term in δΔG(T,P) = δΔE₀ + δΔF_vib(T, V_ref) + δΔG_relax.

## Entries and model

- **Entries.** All 521 entries of the primary selection (`fes-alexandria-polymorph-audit-20261030/metadata/selection/FINAL_GROUPS.csv`), 239 formulas.
- **Model.** DPA-3.1-3M, head `MP_traj_v024_alldata_mixu`. The head is the same as F0. The `Omat24` head is a later sensitivity check, not part of this protocol.

## Procedure, per entry

1. **Start** from the stored DFT unit cell.
2. **Relax** positions and cell at P = 0 with the MLIP:
   - ASE `FrechetCellFilter` plus `FixSymmetry` (symprec 1e-3), so the relaxation keeps the space group and the phase identity;
   - FIRE optimiser, fmax 1e-3 eV/Å, at most 3,000 steps.
3. **Record:**
   - static energy before and after relaxation;
   - relative volume change;
   - maximum atomic displacement after removing the cell change;
   - space group before and after (spglib, symprec 1e-3);
   - whether the relaxation converged;
   - the RMS of the MLIP force on the undisplaced DFT supercell (R₀).
4. **Phonons at the MLIP geometry,** by the same path as `mlip_vs_dft_phonons.py`:
   - phonopy with the entry's stored supercell and primitive matrices, displacement 0.01 Å (phonopy default ± rule), symprec 1e-4;
   - mesh density 100;
   - free energy per atom on the 0–1500 K grid;
   - dense-mesh minimum frequency.

## Pair quantities (meV/atom; MLIP minus DFT unless stated)

- **e_fixed(T).** The existing fixed-geometry error: ΔF_MLIP@DFT-geometry − ΔF_DFT.
- **δ_relax(T).** ΔF_MLIP@MLIP-geometry − ΔF_MLIP@DFT-geometry. It splits into:
  - a static part, δ_relax,static = Δ(E_relaxed − E_at-DFT-geometry);
  - a vibrational part, the remainder.
- **e_pipeline(T) = e_fixed(T) + δ_relax(T).** This is the error of the self-consistent harmonic MLIP pipeline (each method at its own static minimum). No thermal expansion is included on either side.

## Analyses (pre-specified; descriptive, no success criterion)

All use formula-weighted means with a cluster bootstrap over formulas (2,000 draws), at 300, 1000 and 1500 K.

- **R1. Size.** RMS of δ_relax (total, static and vibrational) against RMS of e_fixed.
- **R2. Common mode.** Cancellation fraction c of the relaxation term, single-phase against pair.
- **R3. Interaction.** The uncentred correlation of δ_relax with e_fixed, and the fraction of pairs where |e_pipeline| < |e_fixed|. This asks whether relaxation tends to cancel or to add to the fixed-geometry error.
- **R4. Decisions.**
  - Non-degenerate wrong-ordering rate (|ΔF_DFT| ≥ 1 meV/atom) for e_pipeline against e_fixed.
  - Coverage at α = 0.05 under the held-out-formula abstention rule.
- **R5. Label-free warning.** Spearman correlation of |e_pipeline| with two features that need no phonon labels:
  - the pair difference in relaxation energy, |Δ(E_relaxed − E_at-DFT)|;
  - the pair maximum of R₀.

## Exclusions

These are reported, never silently dropped:
- entries whose relaxation does not converge;
- entries whose space group changes;
- entries with an MLIP-geometry minimum frequency below −0.05 THz.

Pairs that contain such an entry are excluded from R1–R5 and counted separately. The rate of MLIP-geometry dynamical instability is itself reported as a failure mode.
