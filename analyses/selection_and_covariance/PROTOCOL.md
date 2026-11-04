# Thermodynamic error-covariance audit: Alexandria PBE polymorphs, same-geometry DFT vs DPA (protocol)

Written and hashed before any MLIP result was inspected. MLIP numbers seen so far are limited to the two
pipeline-test entries (mp-1000 BaTe, mp-10009 GaTe). Neither belongs to a polymorph group used here, and
neither informed any rule below.

## Scope

- **Task.** Fixed-volume harmonic Helmholtz free energy at the DFT geometry:
  F(T) = E₀ + F_vib(T), split into E₀, ZPE = F_vib(0 K), and F_th(T) = F_vib(T) − F_vib(0 K).
- **Not covered.** QHA, thermal expansion, Gibbs free energies, anharmonicity, MLIP geometry relaxation (a separate later experiment). Results must not be described as Gibbs/QHA results.
- **Reference.** Alexandria PBE phonon benchmark (Loew et al. 2025; CC BY 4.0): VASP PBE static energies and finite-displacement forces as stored.
- **MLIP.** DPA-3.1-3M (SHA-256 86dd3a80…).
  - Primary head: `MP_traj_v024_alldata_mixu`.
  - Sensitivity head: `Omat24`.
  - Platform: thu-GenSi A100, environment `mat-agent` (deepmd dev build, phonopy 2.29). All entries are computed on this single platform.
  - Measured run-to-run noise floor: about 0.06 meV/atom in F_vib(600 K) for a single phase.
- **Shared pipeline.** Both sides go through identical phonopy processing: stored displacements → produce_force_constants → mesh density 100 → thermal properties on a 0–1500 K grid in 50 K steps.
- **Known confound.** MPtrj likely contains relaxations of these equilibrium structures, so static-energy errors may be optimistically small. Displaced configurations are not in MPtrj.

## Entry and group selection (DFT-side information only)

1. **Groups.** `metadata/POLYMORPH_GROUPS.csv`: reduced formula, at least two entries flagged stable in `pbe.csv`, formula resolved through MP OPTIMADE.
2. **Dynamical stability.** The minimum DFT mesh frequency must be ≥ −0.05 THz. Small negative acoustic values are numerical.
3. **Duplicates.** Within a formula, entries with the same space-group number, |ΔE₀| < 1 meV/atom and |ΔV/V| < 0.5% are treated as one structure; the lowest mp-id number is kept.
4. **Final groups.** A formula is kept only if at least two distinct entries survive steps 2–3.
5. **Failures.** Entries whose MLIP or DFT pipeline fails are reported and excluded. No substitution is made.

## Quantities

For each entry (meV/atom; MLIP minus DFT):
- δE₀;
- δZPE;
- δF_th(T) at T ∈ {300, 600, 1000, 1500} K;
- δF(T) = δE₀ + δZPE + δF_th(T).

For each unordered pair (i, j) within a formula, the same components of the difference:
- δΔE₀;
- δΔZPE;
- δΔF_th(T);
- δΔF(T).

## Questions and metrics

**Q1. Cross-phase common mode.** For each component and temperature, compare the RMS pair error with √2 × the RMS single-phase error. The covariance-cancellation fraction is c = 1 − Var(δΔ) / [Var(δ_i) + Var(δ_j)].

**Q2. Cross-component cancellation (the HfO₂ question).** Per pair and temperature, with a = δΔE₀ + δΔZPE (the 0 K anchor) and b = δΔF_th(T):
- the fraction of pairs where a and b have opposite signs;
- the fraction where |a + b| < max(|a|, |b|);
- Spearman correlation of a and b across pairs.

**Q3. Partial-correction harm.** Per pair and temperature, compare |δΔF| under three conditions:
- no correction: a + b;
- anchor corrected only: b;
- thermal corrected only: a.

Report the fraction of pairs where each partial correction increases |δΔF|, and where it changes the sign of the predicted ΔF ordering relative to DFT. The ordering criterion is sign(ΔF) at fixed volume; pairs with |ΔF_DFT| < 1 meV/atom are reported separately as near-degenerate.

**Q4. Temperature stability.** Q1–Q3 at each temperature, and whether the sign of the a–b relation is the same at 300, 600, 1000 and 1500 K.

## Statistics

- Inference is clustered by formula. Temperatures are not independent samples.
- Uncertainty: cluster bootstrap over formulas (2,000 resamples).
- Two weightings are reported: each formula weighted equally (primary), and each pair weighted equally (secondary).
- Large polytype families (CdI₂, ZnS, SiO₂, …) are additionally reported with those formulas excluded.

## Interpretation limits

- This is a single-volume harmonic subtask on mostly non-magnetic compounds.
- A covariance pattern shows exploitable error structure. It does not show that the structure can be predicted or corrected with limited information; that needs a separate held-out test.
