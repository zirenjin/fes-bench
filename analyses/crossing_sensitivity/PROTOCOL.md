# Transition-temperature sensitivity, error-type triage and tail audit

Frozen before any of the quantities below are computed.
- **Status.** Descriptive, on previously analysed data: the Alexandria fixed-geometry harmonic pairs, with F0 (DPA-3.1-3M, MP_traj head) and MatterSim-v1.0.0-5M. No new DFT. The only new model evaluations are re-runs of the existing phonon path for the tail-audit entries listed in part 3.
- **Purpose.** Turn "transition prediction is hard" into a per-pair accuracy budget, and decide which pairs need an offset correction, which need a response correction, and which are outside the harmonic approximation.

## Part 1. Root-finding sensitivity

**Scope.** Every primary pair (351 pairs, 239 formulas) whose DFT reference ΔG(T) = G_A − G_B, with G = E₀ + F_harm(T), has a zero crossing in 50–1500 K. The **first** such crossing is T_c^ref, the same definition as the Tc arm of the anchored-response protocol. ΔS_ref = −∂ΔG_ref/∂T at T_c^ref is taken from the 50 K grid gradient, linearly interpolated.

**Accuracy budget.** For a tolerance τ ∈ {50, 100} K, the allowed local free-energy error is ε_τ = |ΔS_ref| · τ. The model's local error is e_Δ = ΔG_model(T_c^ref), using the unanchored ("none") model curve. Reported:
- the distribution of ε_50;
- the fraction of pairs with |e_Δ| ≤ ε_τ;
- the same under the static-ΔE₀ arm, as secondary.

**Checking the linear relation.** The model transition T_c^mod is the model crossing in 50–1500 K nearest to T_c^ref; it is missing if the model has no crossing. Compared:
- the actual δT_c = T_c^mod − T_c^ref;
- the first-order δT_lin = e_Δ / ΔS_ref;
- the response-aware δT_resp = e_Δ / ΔS_mod(T_c^ref).

Agreement is reported as a Spearman correlation and as the fraction of pairs with |δT_c − δT_lin| ≤ 50 K. A missing model crossing is reported as **missed**, not imputed.

## Part 2. Error-type triage (per model; classes applied in this order, mutually exclusive)

| Class | Rule |
|---|---|
| **A. Applicability** | Either entry has a dense-mesh minimum frequency < −0.05 THz at the DFT geometry, in DFT or in the model being triaged |
| **B. Correct** | A model crossing exists with \|δT_c\| ≤ 50 K |
| **C. Response error** | After the constant Tc anchor (the model curve shifted so ΔG_model(T_c^ref) = 0), the sign is wrong at some grid point in 50–1500 K with \|T − T_c^ref\| ≥ 100 K. This includes a wrong-sign ΔS_mod. A constant offset cannot fix it. |
| **D. Location error** | Neither B nor C: the constant Tc anchor makes every grid point at \|T − T_c^ref\| ≥ 100 K correct, so an offset suffices |

- **Within D,** also report whether the deployable static-ΔE₀ arm alone brings |δT_c| ≤ 50 K.
- **Sub-flag within C:** whether ΔS_mod(T_c^ref) has the wrong sign.
- **Uncertainty.** All fractions are given with a 95% formula-cluster bootstrap (2000 draws).

## Part 3. Tail audit (fixed list, no follow-up selection)

**Entries.** All entries of CsCl, CsBr, LiI and AgI: the four formulas with the largest MatterSim excess at 1000 K in `fes-mattersim-replication-20261102/analysis/PAIR_EXCESS_INFLUENCE.txt`.

**Checks, for DFT, F0 and MatterSim on the identical stored displacements:**
1. **Structure identity.** Space group (spglib, symprec 1e-3), atoms in the unit cell, volume per atom.
2. **Supercell and long range.**
   - Supercell lattice lengths and the shortest periodic-image distance, compared with the model cutoffs (MatterSim 5 Å, DPA 6 Å).
   - NAC is not applied on either side. This is a shared reference approximation, recorded but not counted as model error.
3. **Force constants.**
   - Relative FC error against DFT.
   - The acoustic sum-rule violation (max |Σ_j Φ_ij|) before symmetrisation.
   - The change in F_th at 1000 K when `symmetrize_force_constants()` is applied to both sides.
4. **Spectrum processing.**
   - The fraction of mesh modes that are imaginary, and of modes below 0.5 THz.
   - ⟨ln ν⟩ over real modes.
   - Whether δF_th(1000 K) ≈ 3 k_B T δ⟨ln ν⟩, i.e. whether the error is a frequency-scale (softening) error.
   - The share of F_th from modes below 0.5 THz.
5. **Numerical convergence.** F_th(1000 K) at mesh 100 against mesh 200.

**Reading.** A case is called anharmonic or soft-mode-driven only if checks 1–5 show no processing or identity problem **and** the error is concentrated in low-frequency modes. Otherwise the case is labelled with the problem found.
