# Transition-temperature sensitivity, error-type triage and tail audit: result

- **Protocol.** `PROTOCOL.md`, SHA-256 e1d865a3…, frozen 2026-11-02T17:25Z, before any quantity below was computed.
- **Status.** Descriptive, on previously analysed Alexandria fixed-geometry harmonic pairs. There is no new DFT; the only new model evaluations are re-runs of 11 tail entries.
- **Outputs.** `results/TC_SENSITIVITY.json` (with a per-pair `.pairs.csv`) and `results/TAIL_AUDIT.json`.
- **Scripts.** `scripts/tc_sensitivity.py`, `scripts/tail_audit.py`, `scripts/dump_curves.py`.
- **Deviation.** The first tail-audit run had a units error in two diagnostic columns (h/k_B written as 0.048 instead of 47.99 K/THz) and divided Alexandria's per-primitive-cell free energy by the conventional-cell atom count. Both were fixed and all 11 entries re-run. No protocol-defined quantity changed.

## Part 1. Root-finding sensitivity (101 reference-crossing pairs, 81 formulas)

| Quantity | F0 | MatterSim |
|---|---|---|
| \|ΔS_ref\| at T_c^ref, median (IQR), µeV/(atom K) | 7.5 (2.7–13.9) | same pairs |
| Allowed local error for ±50 K, ε₅₀ = \|ΔS_ref\|·50 K, median (IQR), meV/atom | **0.38 (0.13–0.69)** | same |
| Actual \|e_Δ(T_c^ref)\|, unanchored, median (IQR), meV/atom | **2.7 (0.8–6.4)** | **2.0 (0.9–4.5)** |
| Within the ±50 K budget | 12% [5, 19] | 7% [2, 13] |
| Within the ±100 K budget | 18% [11, 27] | 15% [8, 23] |
| Within the ±50 K budget, static-ΔE₀ arm | 16% [9, 24] | 18% [11, 27] |
| Model crossing found in 50–1500 K | 57% | 61% |
| Median \|δT_c\| where found | 175 K | 183 K |

**How well the linear relation explains the observed shifts** (pairs with a model crossing; F0 n = 60, MatterSim n = 63):
- With the model's own slope, δT_resp = e_Δ/ΔS_mod is within 50 K of the actual shift in **92% and 90%** of pairs.
- With the reference slope, δT_lin = e_Δ/ΔS_ref, the figures are 60% and 52%; the Spearman correlation is 0.92 and 0.87.

The observed transition-temperature errors are therefore almost entirely the root-finding relation: the local free-energy error divided by the slope. The slope error itself matters as well. The typical local error is 5–7 times the ±50 K budget.

## Part 2. Error-type triage (same 101 pairs)

| Class | F0 (pairs; formula-weighted share) | MatterSim |
|---|---|---|
| A. Applicability (either entry < −0.05 THz in DFT or the model) | 19; 16% [9, 25] | 27; 24% [16, 34] |
| B. Correct (\|δT_c\| ≤ 50 K) | 8; 8% [3, 14] | 6; 6% [1, 11] |
| C. Response error (a constant Tc anchor still leaves wrong signs ≥ 100 K away) | 7; 9% [2, 15] | 11; 10% [4, 17] |
| … of which the ΔS_mod sign is wrong | 6 of 7 | 11 of 11 |
| D. Location error (an offset suffices) | **67; 67% [57, 77]** | **57; 60% [49, 70]** |
| … D fixed to ±50 K by the deployable static-ΔE₀ arm | 11 of 67 | 11 of 57 |

The class labels are partly model-specific. Both models put 44 pairs in D. 11 of F0's D pairs are A for MatterSim, and 7 are C for MatterSim.

**Reading.**
- **The dominant failure is location, not response.** In 60–67% of transition pairs the relative-entropy sign is right, and a correct constant offset at the transition would fix the ordering away from it. But that offset must be accurate to about 0.4 meV/atom (median).
- **The deployable static-ΔE₀ correction fixes only 16–19% of the location cases.** The offset at T_c includes integrated thermal error, not only the static term.
- **"A constant suffices" is a statement about the error type, not about the information available.** The class-D test uses the reference T_c.
- **Response errors (about 10%) are almost always a wrong-sign ΔS.** No offset can repair these.
- **16–24% of transition pairs are outside the harmonic validity gate.**

## Part 3. Tail audit (CsCl, CsBr, LiI, AgI; 11 entries)

**Checks that pass for every entry:**

| Check | Result |
|---|---|
| Structure identity | The spglib space group matches the Alexandria record for every entry. CsCl and CsBr are B2 (Pm-3m) against B1 (Fm-3m). |
| Processing | Our DFT F(600 K) matches Alexandria's own processed value within 0.11 meV/atom. The recomputed F_th equals the stored values. |
| Supercell | The shortest periodic image is 12.0–17.5 Å, longer than the per-layer cutoffs (5 and 6 Å). |
| NAC | Not applied on either side, and not in the yaml. This is a shared approximation of the reference, not a model–reference difference. |
| Force constants | The acoustic sum-rule violation is ≤ 4.4e-5 eV/Å². Symmetrising changes F_th(1000 K) by ≤ 0.25 meV/atom. |
| Mesh convergence | Moving from mesh 100 to mesh 200 changes F_th(1000 K) by ≤ 0.2 meV/atom. |

**Where the error is.**
- **It is a frequency-scale error.** δF_th(1000 K) ≈ 3 k_B T δ⟨ln ν⟩ within 0.1–3.2 meV/atom for every entry. The largest gap is F0 CsCl-B2, whose imaginary modes are excluded from the average.
- **For MatterSim, it is concentrated in one polymorph per pair:**

| Formula | Phase | ⟨ln ν⟩, DFT / F0 / MatterSim | Modes < 0.5 THz, DFT / F0 / MatterSim | δF_th at 1000 K, F0 / MatterSim (meV/atom) |
|---|---|---|---|---|
| CsCl | B2 | 0.65 / 0.51 / 0.72 | 0.5 / 3.9 / 0.4% | −34.1 / +16.5 |
| CsCl | B1 | 0.65 / 0.60 / **0.32** | 0.6 / 0.9 / **12.4%** | −13.2 / **−83.9** |
| CsBr | B2 | 0.42 / 0.39 / 0.41 | 0.7 / 0.9 / 1.0% | −6.6 / −2.9 |
| CsBr | B1 | 0.41 / 0.20 / **0.09** | 1.0 / 1.6 / **14.7%** | −51.2 / **−83.4** |
| LiI | P6₃mc | 1.20 / 1.16 / **1.02** | 0.6 / 0.8 / 2.1% | −8.8 / **−43.1** |
| AgI | P6₃mc / F-43m | 0.60 / 0.55 / 0.51 (both) | 1–2 / 1–2 / 2–4% | about −10 / about −21 (both) |
| AgI | P4/nmm | 0.43 / 0.40 / 0.45 | 8.3 / 11.1 / 6.5% | −6.5 / +6.3 |

**Reading.**
- **No processing, identity, supercell or convergence problem was found.**
- **CsCl-B1 and CsBr-B1 under MatterSim.** 12–15% of the modes fall below 0.5 THz against about 1% in DFT, with no imaginary modes. That is a model-predicted near-instability that the DFT harmonic spectrum does not have.
- **These are model softening, not physical anharmonicity.** By the protocol's rule the two cases are low-frequency-concentrated model softening. They are not evidence of anharmonicity in the real systems, because the reference is harmonic and stable here.
- **LiI-P6₃mc is a broader softening of one phase.**
- **AgI.** MatterSim softens the two tetrahedral phases equally, so that error cancels between them. The pair errors involve P4/nmm, which it does not soften.
- **All four cases share one pattern: phase-specific softening.** The model treats one polymorph of a composition differently from the others. This is the same lack of common-mode error that lowered MatterSim's cancellation fraction.
- **Attribution by a 0.5 THz bin is not additive here** (modes move between bins), so no share-of-error number is given.

## What this implies for "how to optimise"

1. **The accuracy requirement is set by ΔS at the transition: about 0.4 meV/atom for ±50 K at the median pair.** Current models are 5–7 times above it, and a lower single-phase MAE does not change that ratio by itself.
2. **Most transition failures need offset information near the decision region, not a better response.** The static ΔE₀ offset is not enough. A deployable source would be a reference query near the model's own predicted T_c. That needs a new frozen protocol and reference computation.
3. **About 10% of transition pairs need the response fixed (a wrong-sign ΔS); about 20% are outside the harmonic gate.** Neither group should be counted as a calibration problem.
4. **The large tail errors are phase-specific model softening, with no processing artefact found.** Polymorph-specific spectral disagreement could be checked as a reference-free flag, for example via ⟨ln ν⟩ or the low-mode fraction compared across polymorphs or across models. It is untested here and would need its own protocol.
