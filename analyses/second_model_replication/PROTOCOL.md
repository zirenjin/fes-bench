# MatterSim model replication on the existing Alexandria pairs

Frozen before any MatterSim force or energy on these entries exists. This is **model replication on previously analysed data, not structure-zero-shot confirmation.** MatterSim's pre-training overlap with these Alexandria structures is unknown.

**Question.** Do the capability differences seen with DPA (F0) also appear in an independent model? The differences are:
- single-phase versus pair error;
- common-mode cancellation;
- temperature-dependent dominance of anchor versus thermal error;
- ordering errors;
- pipeline failure rate;
- the anchored-response pattern.

## Model and settings

- **Checkpoint.** MatterSim-v1.0.0-5M, SHA-256 e3df9fa708725e3d453140646c7d1838324b347a3d1214cf1440522146f872b5. Package `mattersim` 1.2.0, torch 2.10.0+cu128, `MatterSimCalculator(device="cuda")`. No fine-tuning, no relaxation.
- **Evaluation.** `scripts/mattersim_vs_dft_phonons.py`, which is `mlip_vs_dft_phonons.py` with only the evaluator swapped. It uses the identical stored displaced supercells at the DFT geometry, the identical phonopy path (mesh density 100) and T = 0–1500 K in 50 K steps. The output keys are identical.
- **Sample.** The same 521 primary-selection entries (239 formulas, 351 pairs) as the F0 audit. They are not re-selected on MatterSim results.

## Scoring (reused unchanged)

1. **`analyze_covariance.py` (Q1–Q4).**
   - Single-phase and pair RMS.
   - Cancellation fraction c, for δE₀, δZPE and δF_th.
   - Anchor–thermal relation.
   - Ordering under none, anchor-corrected and thermal-corrected.
2. **`abstention_check.py`.** Residual structure, and coverage at α = 0.05.
3. **`anchored_alexandria.py` and `supp_same_subset.py`.** Arms none, static, Tc and T₀ ref. The Tc arm is applied on the pairs where the **reference** has a crossing in 50–1500 K, so the subset is the same as for F0.
4. **Failure rate.** Entries that fail, and entries whose dense-mesh MatterSim minimum frequency is below −0.05 THz at the DFT geometry.

**Primary metric.** As frozen: pair δΔF_th RMS at 300, 600, 1000 and 1500 K (zero-point energy excluded), formula-weighted with a cluster bootstrap. It is reported next to F0.

## Pre-specified replication checks

Each check is marked replicated or not for MatterSim, alongside the F0 value:

| Check | Criterion |
|---|---|
| C1. Common-mode cancellation | δF_th cancellation c > 0.5 at 600, 1000 and 1500 K |
| C2. Decisions not implied by MAE | Uncorrected non-degenerate wrong-ordering rate ≥ 5% at every T |
| C3. Temperature-dependent value of the static anchor | Anchor-only correction reduces wrong ordering at 300 K by more than at 1500 K (absolute points) |
| C4. Anchored response near the decision region | On the reference-crossing subset, the Tc arm's wrong-ordering rate is below the uncorrected rate at 600 and 1000 K |

## Sensitivity subsets (compatibility, not performance)

- **All 239 formulas** (primary; this matches the F0 scoring).
- **Excluding 8 formulas** that contain a Hubbard-U-scheme transition metal without O or F (FeP₄, FeS₂, AsCoS, Li₇N₄V, Li₃S₄V, MoN₄Sr₃, N₃WY, S₂W). Their spin treatment may differ between the reference and the model's training data.
- **Excluding all 71 formulas** with a potentially magnetic element.

No formula is in the MP +U scheme (U-element with O or F).

## Not in scope here

- HfO₂ QHA runs only afterwards, and only after a phonon-validity check. If the strict tier does not hold, lenient-tier curves are not used for model ranking.
- MACE models are a later addition and are not treated as a controlled stepwise intervention.
