# MatterSim model replication on the Alexandria pairs: result

- **Protocol.** `PROTOCOL.md`, SHA-256 d2804a39…, frozen 2026-11-02T05:41Z.
- **Model.** MatterSim-v1.0.0-5M, SHA-256 e3df9fa7… (verified on GenSi before running).
- **Status.** Model replication on previously analysed data. Pre-training overlap is unknown. No re-selection on MatterSim results.
- **Run.** 521 of 521 entries; no failures. The scorers are reused unchanged.
- **Outputs.** `analysis/covariance`, `analysis/ABSTENTION.json`, `analysis/ANCHORED.json`, `analysis/SUPP_SAME_SUBSET.json`, and the sensitivity runs `analysis/cov_excl8`, `analysis/cov_excl71`, `analysis/F0_cov_*`.

## Primary metric and error structure (239 formulas, 351 pairs)

**Weighting.** The RMS values are pooled over entries or pairs, as `analyze_covariance.py` computes them. PROTOCOL.md calls the primary metric formula-weighted; the formula-weighted RMS is reported in the post-hoc section below. The cancellation fraction c and the ordering rates are formula-weighted, with a 95% formula-cluster bootstrap.

| Quantity | F0 (DPA) | MatterSim |
|---|---:|---:|
| Single-phase δE₀ RMS | 46.1 | 45.9 |
| Pair δE₀ RMS | 7.19 | **5.52** |
| Single-phase δZPE RMS | 2.93 | **1.08** |
| Single-phase δF_th RMS, 300 / 600 / 1000 / 1500 K | 3.22 / 7.86 / 14.32 / 22.53 | **2.88 / 6.07 / 10.37 / 15.77** |
| **Pair δF_th RMS** (primary), 300 / 600 / 1000 / 1500 K | **1.77 / 4.00 / 7.15 / 11.25** | 2.63 / 5.59 / 9.67 / 14.91 |
| δF_th cancellation c at 1000 K | 0.71 (0.65–0.77) | 0.56 (0.49–0.63) |
| Relative-entropy error RMS at 600 / 1000 / 1500 K (µeV/(atom K)) | 7.8 / 8.2 / 8.5 | 10.2 / 10.5 / 10.8 |
| Uncorrected wrong ordering, 300 / 600 / 1000 / 1500 K | 11.6 / 10.6 / 11.4 / 12.2% | 11.8 / 8.7 / 11.5 / 11.6% |
| Imaginary modes at the DFT geometry (< −0.05 THz), of 521 entries | 55 | 67 (DFT itself: 13) |

Units are meV/atom unless marked. MatterSim's δF_th cancellation fraction is 0.56 at every temperature.

**The reading (revised after the post-hoc checks below).**
- MatterSim has 26–30% lower single-phase thermal free-energy error than F0, and less of its error is common to competing polymorphs: c is 0.56 against 0.71.
- Its pooled pair RMS is 35–49% higher. Paired bootstrap intervals for that ratio include 1, and the excess is concentrated in a few formulas: CsCl and CsBr (one pair each) carry 84–96% of it.
- In 63–66% of formulas MatterSim has the smaller pair error; the median per-formula ratio is 0.68–0.70.
- Uncorrected ordering-error rates do not differ detectably (paired intervals about ±5 points). The pattern holds on both compatibility subsets:

| Subset | c at 1000 K, F0 / MatterSim | Pair δF_th RMS at 1000 K, F0 / MatterSim |
|---|---|---|
| excluding 8 transition-metal formulas | 0.70 / 0.55 | 7.2 / 9.8 |
| excluding 71 potentially magnetic formulas | 0.69 / 0.55 | 7.8 / 11.0 |

The c difference holds on both subsets. The pooled pair-RMS ranking on these subsets has not been checked for tail influence.

This is a cross-model observation, not an intervention: the two models differ in training data and architecture. It does **not** show a general ranking reversal. What it shows is that the pair-error ranking depends on the aggregation: a tail-sensitive pooled RMS ranks F0 better, while the typical formula ranks MatterSim better. Single-phase accuracy predicts neither the tail nor the cancellation fraction.

## Pre-specified replication checks

| Check | F0 | MatterSim | Replicated |
|---|---|---|---|
| C1. Cancellation c > 0.5 at 600 / 1000 / 1500 K | 0.70 / 0.71 / 0.71 | 0.56 / 0.56 / 0.56 (95% CI lower bound about 0.48–0.49) | point-estimate criterion met; the interval crosses 0.5, so c > 0.5 is not statistically established for MatterSim |
| C2. Uncorrected wrong ordering ≥ 5% at every T | 10.6–12.2% | 8.7–11.8% | yes. This shows that ordering errors persist; it is not by itself evidence that single-phase and pair metrics are not interchangeable. |
| C3. The anchor-only correction helps more at 300 K than at 1500 K | 9.0 against 1.9 points | 10.8 (11.8 → 1.0) against 4.2 (11.6 → 7.4) points | yes |
| C4. Tc arm below uncorrected on the reference-crossing subset at 600 and 1000 K | 1.8 < 16.4; 6.1 < 16.7 | 12.7 < 14.5; 6.1 < 22.7 | point-estimate criterion met; at 600 K the margin is 1.8 points with n = 55, so this is not a strong replication |

Further checks:
- **Anchor–thermal relation.** It replicates: opposite-sign fraction 0.57–0.58 and Spearman −0.18 against F0's 0.55–0.57 and −0.19.
- **Dominance switch.** It replicates. For MatterSim, thermal-only correction already beats anchor-only correction at 1500 K (6.6 against 7.4% wrong).

## Anchored-response arms (same pairs and points as F0)

| Arm | F0 wrong, 600 / 1000 / 1500 K | MatterSim wrong |
|---|---|---|
| none (all 351 pairs) | 11.2 / 10.6 / 10.8% | 9.0 / 11.3 / 11.4% |
| static ΔE₀ | 5.4 / 7.0 / 9.8% | 4.3 / 5.3 / 9.4% |
| T₀ ref, 300 K (diagnostic) | 2.9 / 4.9 / 7.1% | 1.8 / 3.2 / 5.1% |
| Tc subset (101 pairs): none | 16.4 / 16.7 / 7.1% | 14.5 / 22.7 / 15.5% |
| Tc subset: Tc | 1.8 / 6.1 / 4.8% | 12.7 / 6.1 / 9.5% |
| Tc subset: T₀ ref | 7.3 / 9.1 / 6.0% | 5.5 / 7.6 / 10.7% |

- **Tc arm side correctness.** F0 0.92, with 90% of pairs fully correct; MatterSim 0.88, with 86% fully correct.
- **The near-decision-region advantage seen for F0 does not replicate cleanly.** For MatterSim, the 300 K reference anchor beats the Tc anchor at 600 K (5.5 against 12.7%). It loses at 1000 K and 1500 K only narrowly (7.6 against 6.1%; 10.7 against 9.5%).
- **Static ΔE₀.** It again helps mainly at low temperature.

## Residual structure

The MatterSim residuals are weakly anticorrelated with the true thermal gap: Pearson −0.11 to −0.14, against about 0 for F0. MatterSim tends to compress phase-to-phase thermal differences. Heteroscedasticity is similar (Spearman |e|–|D| 0.38–0.40).

## What replicated, and what did not

| Replicated in an independent model | Model-specific or not replicated |
|---|---|
| A substantial residual pair error after cancellation | **How much** cancels: 0.71 for F0, 0.56 for MatterSim |
| Ordering errors of about 10% that are not tracked by single-phase accuracy | The size of the Tc-anchor advantage at 600 K |
| The low-temperature value of a static anchor | |
| The temperature-dependent dominance switch | |
| Near-chance anchor–thermal coupling | |
| Imaginary-mode failures at the DFT geometry | |

**Single-phase phonon accuracy did not predict the pair-error tail or the cancellation fraction in this pair of models.** The pooled pair-RMS "inversion" is real as computed, but it is carried by a few formulas (see below), so it is not a broad ranking reversal.

## Post-hoc checks (requested 2026-11-02; descriptive, not in PROTOCOL.md)

**Scripts and outputs.**
- `scripts/paired_model_contrast.py` produces `analysis/PAIRED_CONTRAST.json`.
- `scripts/pair_excess_influence.py` produces `analysis/PAIR_EXCESS_INFLUENCE.txt`.
- Per-entry minimum frequencies come from the GenSi npz files and are in `analysis/FMIN_ENTRIES.csv`.

**Method.** Both models are scored on the same 2000 formula-cluster draws. The metric definitions are those of `analyze_covariance.py`, and the point estimates reproduce the tables above.

**1. Paired bootstrap** (all 239 formulas). Ratios and differences are MatterSim relative to F0, with 95% intervals.

| T | Single-phase RMS ratio | Pair RMS ratio (pooled) | Pair RMS ratio (formula-weighted) | Δc | Δ wrong ordering (points) | Draws with pooled inversion |
|---|---|---|---|---|---|---:|
| 300 K | 0.89 [0.64, 1.14] | 1.49 [0.98, 2.02] | 1.47 [0.86, 2.03] | −0.12 [−0.21, −0.03] | +0.2 [−5.0, 5.4] | 79% |
| 600 K | 0.77 [0.57, 0.98] | 1.40 [0.92, 1.88] | 1.37 [0.81, 1.91] | −0.14 [−0.23, −0.05] | −1.9 [−6.7, 3.0] | 94% |
| 1000 K | 0.72 [0.55, 0.90] | 1.35 [0.91, 1.81] | 1.32 [0.81, 1.83] | −0.15 [−0.24, −0.06] | +0.2 [−5.8, 6.0] | 93% |
| 1500 K | 0.70 [0.53, 0.87] | 1.33 [0.90, 1.77] | 1.29 [0.82, 1.78] | −0.16 [−0.25, −0.07] | −0.5 [−5.7, 4.5] | 92% |

What the intervals establish:
- **Established:** MatterSim's lower single-phase error at 600–1500 K, and its lower cancellation fraction at every temperature.
- **Not established:** its higher pair RMS. The interval includes 1 at every temperature, under both weightings.

**2. Influence.**
- **Concentration.** CsCl and CsBr (one pair each) carry 84–96% of MatterSim's excess pooled sum of squares; adding LiI brings it to 99–117%. A share above 100% means the remaining formulas favour MatterSim.
- **Dropping the largest-error formulas** (largest in either model, so not selected on the contrast) brings the pooled ratio to 0.92–1.09 for three formulas and 0.71–0.76 for ten.
- **Per formula,** MatterSim's pair error is larger in only 34–37% of the 239 formulas. The median per-formula ratio is 0.68–0.70.

**3. Common harmonically stable subset.**
- **Gate.** The dense-mesh minimum frequency at the DFT geometry must be ≥ −0.05 THz for DFT, F0 and MatterSim alike. This is the threshold already used above, fixed before the check was run.
- **Failures.** DFT fails 13 entries, F0 55 and MatterSim 67; 114 fail at least one.
- **Coverage.** 364 of 521 entries (70%), 219 of 351 pairs (62%) and 171 of 239 formulas (72%) remain.

| T | Single-phase RMS F0 / MS (ratio) | Pair RMS F0 / MS, pooled (ratio) | Pair RMS ratio (formula-weighted) | c F0 / MS (Δc) | Wrong ordering F0 / MS |
|---|---|---|---|---|---|
| 600 K | 6.94 / 5.19 (0.75 [0.48, 0.98]) | 3.75 / 5.22 (1.39 [0.90, 1.76]) | 1.25 [0.75, 1.55] | 0.69 / 0.54 (−0.16 [−0.26, −0.05]) | 11.9 / 9.7% |
| 1000 K | 12.83 / 9.02 (0.70 [0.47, 0.92]) | 6.61 / 8.94 (1.35 [0.87, 1.72]) | 1.21 [0.73, 1.51] | 0.70 / 0.53 (−0.18 [−0.28, −0.06]) | 10.8 / 13.4% |
| 1500 K | 20.30 / 13.83 (0.68 [0.46, 0.89]) | 10.21 / 13.59 (1.33 [0.85, 1.69]) | 1.19 [0.72, 1.49] | 0.71 / 0.53 (−0.18 [−0.28, −0.06]) | 10.4 / 12.8% |

The pattern on the common stable subset is the same as on all formulas: lower single-phase error and lower c for MatterSim, and a pooled pair-RMS ratio above 1 whose interval includes 1. CsCl (min frequency −0.33 THz) leaves the subset. CsBr, LiI and AgI pass the gate and stay in it, so the subset does not remove the tail.

**Revised claim for the paper.**
- **Supported:** MatterSim, though more accurate per phase, has a smaller fraction of error common to competing polymorphs.
- **Supported:** its pair error has a heavier tail, concentrated in a few alkali- and silver-halide pairs. The pooled pair-RMS ranking of the two models therefore depends on those formulas.
- **Not supported:** a robust, broad reversal of the model ranking between the single-phase and pair tasks.

HfO₂ QHA for MatterSim is deferred, per protocol, until a phonon-validity check is done.
