# Anchored temperature-response test (development; F0 only)

Frozen before any anchored quantity is computed on any of the three datasets. All three datasets have been seen in earlier analyses, so every result here is **development evidence, not confirmation**.

**Question.** Given a realistically obtainable anchor, over what range does F0's own temperature response support correct phase-competition decisions?

**Model.** F0 only (DPA-3.1-3M, head `MP_traj_v024_alldata_mixu`, or the generic-descriptor readout for SiO₂). F1 and F1b are not candidates, so there is no model selection.

## Arms and information

Each arm adds a per-phase constant, G_p^corr(T) = G_p^F0(T) + b_p. The table states which pair differences b_p − b_q each anchor identifies.

| Arm | Information supplied | Role | Identified biases |
|---|---|---|---|
| none | no target reference label | baseline | none |
| static | ΔE₀ of each phase under the reference Hamiltonian (two reference energies per pair) | low-cost correction | every phase with a reference E₀: b_p = E₀,p^ref − E₀,p^F0 |
| Tc | the reference's first coexistence temperature T_c for one specified pair (A, B), at the dataset pressure | conditional prediction, **reference-sourced coexistence calibration** | b_A − b_B only, such that ΔG_corr(T_c) = 0. Other phases stay unanchored, except where they are reference-labelled sources (SiO₂, below). |
| T₀ ref | reference ΔG at T₀ for every phase against a fixed base phase | **diagnostic control** with high information. In the harmonic task it costs as much as the whole reference curve. | all phases |

The anchored crossing in the Tc arm is guaranteed by construction. It is marked "provided" and never counted as a predicted hit or as T_c accuracy.

## Metrics, per arm, at the evaluation temperatures

**Pair residual and slope.**
- e(T) = ΔG_corr(T) − ΔG_ref(T), as RMS.
- Relative-entropy error, −de/dT, as RMS. Constant anchors leave it unchanged, so it is reported once per dataset as the response quality.

**Ordering.** Non-degenerate pairs only (|ΔG_ref| ≥ 1 meV/atom).
- Wrong-ordering rate.
- Coverage at α = 0.05 under the held-out-formula abstention rule (Alexandria only, where there are enough formulas).

**Tc arm only.**
- **Side correctness.** The fraction of evaluation temperatures with |T − T_c| ≥ 100 K at which the corrected sign matches the reference.
- **Extra crossings.** Predicted crossings of the anchored pair not within ±100 K of a reference crossing.
- **Unprovided events.** Reference crossings other than the provided one: hit within ±100 K, missed, or false.

**Events (all arms).** Reference against predicted crossings in the evaluation window, matched within ±100 K, with the provided event excluded from hits.

**Applicability coverage.** The fraction of pairs on which each arm can be applied. For example, the Tc arm needs a reference crossing in the window. Pairs without one are reported, not dropped.

## Datasets

### 1. Alexandria (fixed-geometry harmonic, PBE)

- **Pairs.** All primary-selection pairs: 239 formulas, 351 pairs. Grid 0–1500 K in 50 K steps.
- **Evaluation temperatures.** 600, 1000 and 1500 K.
- **Arms.** None; static (E₀ from the Alexandria DFT energy); Tc (first crossing of ΔF_DFT(T) in 50–1500 K, linearly interpolated, per pair); T₀ ref at 300 K.

### 2. HfO₂ (direct QHA, F0, against DFT-QHA)

**Validity.** The minimum frequencies were seen before this protocol was written; no G(T) was. Two pre-specified tiers are reported side by side:
- **strict:** dense-mesh minimum frequency ≥ −0.1 THz at V/V₀ = 1;
- **lenient:** ≥ −2.0 THz at V/V₀ = 1, the level of the DFT's own m-HfO₂ artefact (Alexandria mp-352: −1.74 THz).

Only phases passing a tier enter that tier's comparisons, and failure coverage is reported.

**Window and evaluation temperatures.** 300–2500 K; 600, 1000, 1500 and 2000 K.

**Arms.**
- None.
- Static: **not available**, because the reference package gives no separate static energies. The 0 K G includes ZPE and is not substituted.
- Tc: pair t−m at its reference first crossing. This identifies only b_t − b_m; other pairs are reported as not identified.
- T₀ ref at 300 K: every phase against m.

### 3. SiO₂ (ds62 tier; frozen generic-DPA readout, β-cristobalite held out)

**Configurations.** Arm B (DPA) under both scalings, for all four source configurations. Primary: coordinate scaling with `low_quartz_beta+tridymite_p63mmc` and `high_all_three`, as in VERDICT.

Source phases are reference-labelled inputs (b = 0). Only the target cristobalite is predicted, so a single anchor identifies b_crist, and every cristobalite pair is then identified.

**Window and evaluation temperatures.** Native window 851–2499 K; evaluation at 1000, 1500 and 2000 K.

**Arms.**
- None.
- Static: **not available** (ds62 tier).
- Tc: cristobalite − quartz at its reference first crossing in the window.
- T₀ ref at 851 K, the lowest grid temperature, since 300 K is outside the window.

Envelope events are scored on the VERDICT pool `quartz_beta;tridymite_p63mmc;cristobalite_beta`.

## Reporting

- Each dataset is reported separately; nothing is pooled across datasets.
- Failures of every arm are attributed separately to:
  - **response** (slope error);
  - **identity or applicability** (not applicable, or an invalid QHA phase);
  - **under-constrained calibration** (pairs not identified).
- Ordering alone is never used to attribute a failure.
