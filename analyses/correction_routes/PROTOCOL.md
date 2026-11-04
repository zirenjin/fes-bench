# Four correction routes for phase-relative free energies (development validation)

Frozen before any route is fitted or any new spectrum is computed.
- **Status.** All data have been seen before, so every result is **development validation**.
- **Weights and budget.** MLIP weights are frozen (F0 = DPA-3.1-3M MP_traj head; MS = MatterSim-v1.0.0-5M). No new DFT, and no potential is retrained.
- **New model evaluations.** The only new evaluations re-run the frozen models on the stored displacements, to record mode spectra.

**Goal.** Lower the phase-relative free-energy error near transitions and improve phase ordering, without creating tail failures. A lower generic or single-phase MAE does not count as success.

## Data

- **Pairs.** The Alexandria primary set: 521 entries, 239 formulas, 351 same-composition pairs. Fixed DFT geometry, harmonic. G_p(T) = E₀,p + F_harm,p(T), with ZPE included, on T = 0–1500 K in 50 K steps.
- **Spectra** (new; `scripts/spectra.py`).
  - For each entry and model, the mesh-100 frequencies are taken for DFT and the model on the same q-points, sorted per q (band-matched).
  - Modes are binned by model ln ν into 240 log bins over 0.05–60 THz. Each bin stores its weight per atom, the mean ln ν_model, and the mean (ln ν_DFT − ln ν_model).
  - Model modes ≤ 0.05 THz are excluded, as phonopy does. Their weight is recorded.
  - **Gate.** Re-integrating the bins must reproduce the stored model F(T) within 0.05 meV/atom at 1000 K for at least 99% of entries. Otherwise route C is reported as not runnable.

## Information settings (identical for all routes)

- **S0, no reference information.** E₀ and F come from the model.
- **S1, common static-ΔE₀ control.** The DFT static energy replaces the model's for every phase, including test phases. F still comes from the model.

No target-formula reference label is ever an input: no T_c, no reference F or E in S0, no reference spectrum.

## Features (per phase; the same set for routes A–C of a model)

- **Model outputs.** E₀, F at 0/300/600/1000/1500 K, S at 300/1000 K, and spectral summaries: ⟨ln ν⟩, sd(ln ν), the fraction of modes below 0.5 and 1 THz, the 99th-percentile ln ν, f_min and the excluded-mode weight.
- **Structure** (the given DFT geometry). Volume per atom, density, atoms per primitive cell, space-group number and crystal-system index.
- **Composition** (formula-level). Mean and standard deviation of Z, mass and Pauling electronegativity, and the number of elements.
- **Within-formula centred copies** (x_p minus the formula mean) of every model-output and structure feature.

## Routes (MLIP weights frozen in all)

**A. Location only.** c_p is constant in T.
- Target: the formula-centred mean over 50–1500 K of the residual r_p(T) = G_ref,p − G_model,p. This minimises the uniform pair-curve loss exactly.

**B. Location + uniform frequency scale.** c_p(T) = α_p + δ_p U_p(T), where U_p = F_p − T ∂F_p/∂T is the model's vibrational internal energy. This is the first-order correction for scaling every frequency by e^δ.
- Stage-1 target δ̃_p: a per-phase least-squares fit to r_p(T) and to 100 K·(the entropy residual), with a free constant.
- Stage 2: α as in A, on the residual left after the cross-fitted δ.

**C. Spectral correction, then re-integration.**
- Correction: ln ν′ = ln ν + Σ_k θ_p,k B_k(ln ν), where the B_k are four piecewise-linear hats at 0.5, 2, 6 and 20 THz (constant beyond the end knots).
- Stage-1 target θ̃_p: a weighted least-squares fit of the band-matched ln-ratio over bins, plus a 1e-3 ridge. This is a spectral target, not a free-energy target.
- The corrected F(T) is re-integrated exactly from the bins.
- Stage 2: α as in A, on the residual left after the cross-fitted spectral correction.

**D. Model combination.** ΔG_D = w_f ΔG_F0 + (1 − w_f) ΔG_MS, with **one weight per formula** shared by all its phases. Variants:
- **D1:** w = 0.5.
- **D2:** a global w on the grid 0, 0.1, …, 1.
- **D3:** w_f = sigmoid(v·z_f). The formula-level features z_f are the two models' pair-ΔG disagreement at 300/1000/1500 K, their ⟨ln ν⟩ disagreement, and composition. The fit minimises the pair-curve loss with an L2 penalty.

**Learners for the A–C stage targets.** Ridge regression on standardised features, or gradient boosting with depth 2, learning rate 0.05 and subsample 0.8.
- The learner choice and its hyperparameters are tuned in the inner loop: ridge α ∈ {0.1, 1, 10, 100}; boosting {100, 300} trees; D2's w; D3's penalty ∈ {0.1, 1, 10}.
- **Predictions for test formulas are centred within each formula**, because only relative corrections are applied.

## Validation

- **Outer loop.** 5 folds by formula. Every phase and temperature of a formula stays in one fold; the seed is fixed. Three repeats (seeds 1, 2, 3).
- **What is reported.** Primary numbers are from repeat 1. Stability is the range over the three repeats.
- **Inner loop.** 4 grouped folds within the outer-training formulas. All tuning happens here, and the inner out-of-fold predictions supply the cross-fitted stage-1 outputs and the conformal thresholds.
- **Comparison.** All routes are evaluated on the same outer folds, against "none" (the uncorrected model) in the same information setting.

## Metrics (outer out-of-fold predictions)

1. **Near transition (primary 1).** On the 101 pairs with a reference crossing in 50–1500 K: the formula-weighted RMS of e_Δ(T_c^ref), and the formula-weighted fraction within the ±50 K budget ε₅₀ = \|ΔS_ref\|·50 K.
2. **Ordering (primary 2).** The formula-weighted wrong-ordering rate at 600, 1000 and 1500 K on non-degenerate pairs (\|ΔG_ref\| ≥ 1 meV/atom), and its mean over the three temperatures.
3. **Pair error.** Pair ΔG RMS at 300, 600, 1000 and 1500 K.
4. **Response.** Pair δΔS RMS at 600, 1000 and 1500 K.
5. **Events.** Reference crossings hit within ±100 K, missed crossings, false crossings, and the median \|δT_c\|.
6. **Tail.**
   - The 95th percentile of \|e(1000 K)\|.
   - The number of pairs with max_T \|e\| > 20 meV/atom.
   - The number of pairs made worse than "none" by more than 10 meV/atom at any of 300/600/1000/1500 K.
7. **Coverage.** Conformal abstention at α = 0.05, with thresholds from the inner out-of-fold residuals of the training formulas, at 600, 1000 and 1500 K.
8. **Fit and transfer.** Training-set against outer out-of-fold R² of the stage targets, and the pair-curve loss on training against test.

**Uncertainty.** Each route minus "none" is compared on a formula-cluster bootstrap (2000 draws) with the same draws for both.

**Development success** requires, against "none" with the same model and setting, all three of:
- metric 1's RMS is lower, with a 95% interval excluding 0;
- the mean wrong-ordering rate is not higher (point estimate);
- at most 2% of pairs (7 or fewer) are worsened by more than 10 meV/atom, and the 95th percentile of \|e(1000 K)\| is not higher.

Each route is labelled success, partial (some criteria met) or failure. Each failure is attributed to one of three causes:
- the correction does not transfer (good training fit, no out-of-fold gain);
- the model cannot fit even the training formulas;
- the correction creates tails.
