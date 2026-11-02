# Thermodynamic error-covariance audit: result (primary head `MP_traj_v024_alldata_mixu`)

## Status and scope

- **Protocol.** `PROTOCOL.md` was frozen before any group-level MLIP result (SHA-256 85f16e6b…, 2026-10-31T06:40Z). Deviations D1–D4 are in `DEVIATIONS.md`; none changed an analysis rule after MLIP results were seen.
- **Task.** Fixed-volume harmonic Helmholtz F = E₀ + F_vib at the DFT geometry. This is not a Gibbs, QHA or anharmonic result.
- **Reference.** Alexandria PBE phonon benchmark (Loew et al. 2025, CC BY 4.0).
- **Model.** DPA-3.1-3M on thu-GenSi (A100).
- **Run.** 1,849 of 1,850 entries computed; mp-730190 failed out-of-memory and is not in any final group.
- **Selection.** DFT-only rules: dynamically stable (minimum mesh frequency ≥ −0.05 THz), duplicates removed, at least two phases per formula.
  - Primary: **239 formulas, 521 entries, 351 pairs.**
  - Sensitivity, excluding formulas with five or more phases: 235 formulas, 295 pairs.
  - Sensitivity, threshold −0.2 THz: 286 formulas, 446 pairs.
- **Statistics.** Formula-weighted means with cluster-bootstrap 95% intervals over formulas.
- **Agreement.** All three selections agree to within the intervals.

The DFT-side selection already exposed one data-contract finding: 1,035 of the 1,850 entries flagged stable in the dataset have imaginary branches on a dense mesh. 673 of them are below −0.5 THz.

## Q1: common-mode error across phases (cancellation in same-composition differences)

| Component | Single-phase RMS (meV/atom) | Pair RMS (meV/atom) | Cancellation fraction c (95% CI) |
|---|---:|---:|---|
| δE₀ | 46.1 | 7.2 | 0.84 (0.78–0.89) |
| δZPE | 2.9 | 0.9 | 0.77 (0.71–0.82) |
| δF_th, 300 K | 3.2 | 1.8 | 0.68 (0.61–0.74) |
| δF_th, 600 K | 7.9 | 4.0 | 0.70 (0.64–0.76) |
| δF_th, 1000 K | 14.3 | 7.2 | 0.71 (0.65–0.77) |
| δF_th, 1500 K | 22.5 | 11.3 | 0.71 (0.64–0.77) |

- Single-phase δE₀ includes per-element energy-gauge offsets. These cancel exactly within a composition, so part of its high c is trivial.
- The vibrational cancellation is not trivial. About 70% of the single-phase MLIP thermal error is shared by competing polymorphs of the same formula, at every temperature. This is shared per-formula error, not a global softening: the single-phase MLIP-vs-DFT F_th slope is only 1.02–1.04 (see the shrinkage check below).
- The remainder is still large: the pair δF_th RMS is 4 meV/atom at 600 K and 11 meV/atom at 1500 K.

## Q2: cross-component cancellation (the HfO₂ question)

Here a = δΔE₀ + δΔZPE is the 0 K anchor error and b = δΔF_th(T) the thermal error. The values are the same to two decimals at 300, 600, 1000 and 1500 K.

- **Opposite signs.** 0.55–0.57 of pairs (CI about 0.49–0.63), barely above the 0.5 expected by chance.
- **Correlation.** Spearman on the symmetrized set is −0.19 at every temperature. The uncentred correlation is −0.06 to −0.08.

**The anchor and thermal errors are almost unrelated.** When they cancel, as in HfO₂ m → t, it is mostly coincidence and not a systematic error structure of this MLIP.

## Q3: partial correction and decision errors

Pair RMS (meV/atom): a = 7.1 at all temperatures. b grows with temperature: 1.8, 4.0, 7.2 and 11.3 at 300, 600, 1000 and 1500 K.

| T | Anchor-only correction increases \|error\| | Thermal-only correction increases \|error\| | Wrong ordering, no correction | Wrong ordering, anchor corrected | Wrong ordering, thermal corrected |
|---|---|---|---|---|---|
| 300 K | 18% | 49% | 11.6% | **2.6%** | 11.5% |
| 600 K | 29% | 44% | 10.6% | **4.7%** | 9.4% |
| 1000 K | 36% | 35% | 11.4% | 7.5% | 10.2% |
| 1500 K | 41% | 30% | 12.2% | 10.3% | **7.8%** |

Ordering uses non-degenerate pairs (|ΔF_DFT| ≥ 1 meV/atom; 278–297 pairs). The 95% intervals are in `analysis/primary/RESULTS.json`.

- Partial corrections make an individual pair worse often: 18–49% of pairs, depending on the component and the temperature.
- In aggregate, correcting the dominant component reduces ordering errors.
- **The dominant component switches with temperature.**
  - Below about 1000 K the static anchor dominates. Correcting only ΔE₀ + ΔZPE cuts wrong orderings from 11.6% to 2.6% at 300 K.
  - Above about 1000 K the thermal term dominates, and the static correction loses most of its value.

## Q4: temperature stability

The weak negative anchor–thermal relation has the same sign and size at every temperature. The cancellation fractions in Q1 are flat from 600 to 1500 K.

## Implications (descriptive; one MLIP, one dataset, one fixed-volume harmonic task)

1. A phase-relative treatment already gets most of the available cancellation: about 70% of the vibrational error is common mode. A method built around "cross-phase correlation" adds little beyond the plain MLIP ΔF.
2. No exploitable anchor–thermal coupling was found. A joint-correction rule based on predicting cancellation is not supported on this data. That cancellation is mostly fortuitous is itself the finding.
3. The actionable structure is the **temperature-dependent dominance of error components**. This points to a minimal-information rule, still to be tested on held-out systems: a static DFT ΔE₀ query has high decision value at low temperature, and vibrational information becomes necessary above about 1000 K.

**Correction to implication 3.** The "anchor corrected" arm in Q3 corrects ΔE₀ + ΔZPE, and ZPE needs DFT phonons. With static DFT ΔE₀ only, wrong orderings are 3.6 / 4.9 / 7.5 / 11.4% at 300 / 600 / 1000 / 1500 K, against 2.6 / 4.7 / 7.5 / 10.3% for ΔE₀ + ΔZPE. The low-temperature conclusion still holds.

## Omat24 head sensitivity (protocol analysis, `analysis/omat24/`)

Same 239 formulas and 351 pairs. Nothing is missing (mp-730190 is out-of-memory again and is in no group).

- **Q2 replicates.**
  - Opposite-sign fraction 0.53–0.55.
  - Spearman −0.15; uncentred correlation −0.06 to −0.07.
- **Q3 replicates the dominance switch.** Wrong orderings (none / anchor corrected / thermal corrected):

  | T | None | Anchor corrected | Thermal corrected |
  |---|---:|---:|---:|
  | 300 K | 14.3% | 1.0% | 13.8% |
  | 1000 K | 11.2% | 6.3% | 9.5% |
  | 1500 K | 14.5% | 10.6% | 7.9% |

- **Q1 differs in size.** The vibrational cancellation fraction is lower, about 0.47 (δE₀ 0.66), so common-mode cancellation depends on the head. Pair RMS: a = 5.6; b = 1.4 / 3.2 / 5.9 / 9.4 meV/atom.

## Exploratory: coupled anchor/thermal correction (post-hoc, not in PROTOCOL.md)

Script: `scripts/coupled_correction_check.py`. Outputs: `analysis/exploratory/COUPLED_CHECK_{primary,omat24}.json`.

Couplings are fitted leave-one-formula-out, through the origin. The metric is the same non-degenerate wrong-ordering rate.

- **Cross term.** 2⟨ab⟩/⟨(a+b)²⟩ is −0.04 to −0.07 (primary) and −0.03 to −0.08 (Omat24). The anchor–thermal covariance term carries under 10% of the pair error variance.
- **M3s: static ΔE₀ plus a learned coupling for the remainder.** Fitted slopes are −0.06 to −0.23. Against static-only correction the gain is none or marginal:
  - primary, 1500 K: 11.4 → 9.2%, CI on the difference −4.4 to −0.5 points;
  - all other temperatures and both heads: CIs include 0.
  - RMS does not improve.
- **M3z: zero target labels, MLIP-only pair features.** No improvement over no correction. At 1000 K (primary) it is worse, by 0.0 to +3.5 points.
- **M3sz: static ΔE₀ plus an MLIP-feature rescaling of the remainder.**
  - RMS improves at every temperature. The fitted thermal error is about +0.3 × the MLIP's own ΔF_th. The shrinkage check below shows this is regression shrinkage, not a scale bias.
  - Ordering improves only at 300 K: 3.6 → 1.7% (primary) and 2.2 → 0.5% (Omat24). Both CIs exclude 0.
  - At higher temperatures the ordering CIs include 0.

**Reading.** On this data a coupled correction does not beat correcting the dominant component, because there is almost no coupling to exploit. The M3sz gain is a shrinkage estimator: it pulls a noisy MLIP ΔF_th toward zero. It is not a physical scale error and not an anchor–slope coupling.

**Pending.** A held-out check is still needed of whether the dominance crossover can be predicted without reference labels.

## Exploratory: scale bias or shrinkage? (post-hoc; `scripts/shrinkage_check.py`, `analysis/exploratory/SHRINKAGE_CHECK.txt`)

All fits are through the origin. Pair ΔF_th is computed at the DFT geometry.

| Head | Single-phase slope, MLIP on DFT | Pair slope, MLIP on DFT | Pair slope, DFT on MLIP | Pair r |
|---|---|---|---|---|
| Primary | 1.02–1.04 | 0.97–1.00 | 0.69–0.70 | 0.82–0.83 |
| Omat24 | 1.01 | 1.00–1.03 | 0.76–0.78 | 0.87–0.89 |

**MLIP pair ΔF_th has no scale bias.** Its slope on the DFT value is about 1, so it behaves like the DFT value plus scatter that is roughly independent of the signal.

The reverse slope of about 0.7 is regression dilution. It equals the signal fraction of the variance, Var(DFT)/Var(MLIP), at the population level. α ≈ 0.3 is therefore a property of this population's spread of ΔF_th relative to the MLIP scatter. It should not be read as phonon physics, and it transfers only to target populations with a similar spread.

## Exploratory: residual structure and held-out-formula abstention (post-hoc; `scripts/abstention_check.py`, `analysis/exploratory/ABSTENTION_*.json`)

### Residual structure (ε = MLIP − DFT, pair level)

- **Uncorrelated with the true gap.** The symmetrized Pearson correlation of ε with ΔF_DFT is between −0.04 and +0.07 for both the thermal and the total gap, on both heads.
- **Not independent of it.** ε is heteroscedastic: the Spearman correlation of |ε| with |ΔF_DFT| is +0.24 to +0.34. For the total gap, the RMS of ε rises from about 4 to 9–17 meV/atom from the smallest to the largest |gap| tertile.
- **One thermal error per pair.** The correlation of ε_th across temperatures is 0.91–0.997. The thermal error is close to a single per-pair quantity growing with T, so per-temperature rates are not independent evidence.
- **Not interpretable: the between-formula variance share.** Most formulas contribute a single pair, and the pairs within a formula are linearly dependent (e_ik = e_ij + e_jk). The cluster bootstrap over formulas remains the valid unit.

### Abstention

Rule: for a held-out formula, the threshold is the (1 − α) quantile of |ε| over the other formulas. A pair is decided only if |ΔF_MLIP| exceeds the threshold. The metric uses non-degenerate pairs.

| Head, arm | T | No abstention, wrong | α = 0.05 (T-aware): coverage | Selective risk | P(decided ∧ wrong) |
|---|---|---:|---:|---:|---:|
| Primary, raw MLIP | 300 K | 13.0% | 44% | 0.8% | 0.4% |
| Primary, raw MLIP | 1500 K | 10.8% | 30% | 1.1% | 0.3% |
| Primary, + static DFT ΔE₀ | 300 K | 3.5% | 70% | 0.0% | 0.0% |
| Primary, + static DFT ΔE₀ | 1500 K | 9.8% | 31% | 3.3% | 1.0% |
| Omat24, raw MLIP | 300 K | 15.1% | 47% | 1.5% | 0.7% |
| Omat24, + static DFT ΔE₀ | 300 K | 2.1% | 75% | 0.5% | 0.4% |
| Omat24, + static DFT ΔE₀ | 1500 K | 9.8% | 35% | 2.9% | 1.0% |

(The no-abstention rates differ slightly from Q3, which is formula-weighted. These are pair-weighted.)

- **The bound held everywhere.** P(decided ∧ wrong) ≤ α in every cell (16 arm-temperature cells × 3 α × 2 threshold modes × 2 heads). The largest value is 0.034 at α = 0.10. The bound is loose, because a wrong decision needs |ε| > |ΔF_MLIP|.
- **The cost is coverage.** At α = 0.05, raw MLIP decides only 30–47% of non-degenerate pairs. A static DFT ΔE₀ raises this to 70–75% at 300 K, but not at 1500 K (31–35%).
- **T-aware thresholds matter when the anchor is supplied.** With static ΔE₀ at 300 K, the T-aware threshold covers 0.70 against 0.41 for the pooled threshold. At 1500 K the pooled threshold under-protects: joint wrong is 0.027–0.034 against 0.014–0.020 at α = 0.10.
- **Ranking by |ΔF_MLIP| is close to ranking by the true |ΔF_DFT|.** AURC is within about 1 point, sometimes better. Note that sorting by the true gap is not the optimal ranking.

## Exploratory: frequency-band decomposition of the F0 thermal error (post-hoc; `scripts/freqband_decomp.py`, `scripts/freqband_analyze.py`, `analysis/exploratory/FREQBAND_MP_traj.json`)

Method:
- Modes are matched by sorted index at each q-point.
- Each mode's thermal free-energy error is assigned to the band of its DFT frequency.
- Closure against the fixed-geometry run at 1000 K: median 0.012, p95 0.052 meV/atom.
- Sample: 521 entries, 351 pairs.

| Band (THz) | Mode fraction | Mean δν/ν (MLIP − DFT) | Pair-error share, 300 K | 1000 K | 1500 K | Band cancellation c, 1000 K |
|---|---:|---:|---:|---:|---:|---:|
| < 1 | 0.02 | −24.2% | 0.19 | 0.26 | 0.30 | 0.14 |
| 1–2 | 0.09 | −6.8% | 0.37 | 0.28 | 0.26 | 0.48 |
| 2–4 | 0.20 | −4.0% | 0.26 | 0.22 | 0.20 | 0.80 |
| 4–8 | 0.28 | −3.6% | 0.16 | 0.18 | 0.18 | 0.85 |
| 8–16 | 0.28 | −3.2% | 0.02 | 0.05 | 0.05 | 0.96 |
| > 16 | 0.13 | −2.4% | 0.00 | 0.01 | 0.01 | 0.97 |

- **The MLIP is softer in every band, and the relative softening grows sharply at low frequency.**
- **Single-phase error is spread roughly in proportion to mode count.** At 1000 K the shares are 0.07 / 0.15 / 0.23 / 0.28 / 0.24 / 0.03.
- **Pair (phase-relative) error is dominated by the softest modes.** Below 2 THz, 11% of modes carry 54–56% of the pair error. The high-frequency error is almost entirely common mode (c 0.96–0.97). The low-frequency error is phase-specific (c 0.14–0.48).
- **A 19-entry preview had suggested the opposite; it was a small-sample artefact.**
- **Implication for fine-tuning objectives.** An unweighted force-response loss spends most of its weight on modes whose errors cancel between phases. Weighting by contribution to phase-relative free energy (about 1/ω², concentrated below about 2 THz) targets the part that does not cancel.

### Robustness and corrections (post-hoc; `scripts/freqband_robust.py`, `analysis/exploratory/FREQBAND_ROBUST.json`)

All values are formula-weighted, with 95% cluster-bootstrap intervals over formulas.

- **Correction: what the shares measure.** The band shares above are cross-term-inclusive squared-error attributions, sum(c_b x)/sum(x²). They are not independent shares, and they do not equal the error that fixing a band would remove.
- **The direct counterfactual.** This is the mean-square pair error removed if the < 2 THz mode errors alone were zero, 1 − ⟨(x − c_<2)²⟩/⟨x²⟩:
  - 0.67 [0.46, 0.81] at 300 K;
  - 0.62 [0.37, 0.78] at 1000 K;
  - 0.63 [0.38, 0.79] at 1500 K.

  That is roughly a 40% RMS reduction, with wide intervals.
- **Low-frequency cutoff.** Dropping modes with ν_DFT < 0.1 or < 0.3 THz on both sides leaves all numbers unchanged to 0.01. The result is not driven by near-zero acoustic modes or by the cutoff.
- **Acoustic branches.** Dropping the three lowest branches at every q-point leaves 5.85 of the 7.14 meV/atom pair RMS. The < 2 THz removal is then 0.45 [0.24, 0.61] at 1000 K. Acoustic branches carry part of the effect, but low-lying optical modes also contribute.
- **Matching-free binning** (each side's modes binned by their own frequency) is **not interpretable** as an attribution. Systematic softening moves MLIP modes across bin edges, which creates large opposite-sign bin terms (shares > 1, negative removal). The sorted-index matching is therefore not validated by this check. It remains a spectral attribution, not a statement about the same physical mode.
- **Retracted wording.** "An unweighted force-response loss spends most of its weight on modes whose errors cancel" is a hypothesis; no loss or gradient decomposition was done. The data support a low-frequency-sensitive objective as a direction. They do not establish 1/ω² as the right weight.
- **Status of these analyses.** They used all 239 primary formulas. Any objective designed from them (for example F1c) is mechanism development on the current test split, not an untouched confirmation.
