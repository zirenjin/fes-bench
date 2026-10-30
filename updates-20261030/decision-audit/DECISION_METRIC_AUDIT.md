# Post-hoc decision-metric audit: which error metrics describe phase-transition hits and misses?

2026-10-30. **Post-hoc descriptive evidence only. Not a confirmation, not pre-registered.** No new DFT, MD, TI, training or remote job; no model was refit, reselected or rescored.

## Data and integrity

Every row is one (frozen model, candidate pool, reference event). Matching uses the frozen identity-aware rule: same from/to phases, within 100 K.

| Source | Reference tier | Predictions | Verification |
|---|---|---|---|
| SiO₂ A/B/C, Si/O task checkpoint (96 models, 4 held-out targets) | DFT PBE-D3(BJ), DaRUS-4999 | regenerated with the frozen `run_abc.py` | gap MAE ≤ 4.5×10⁻¹¹ relative; 0 event-list mismatches |
| SiO₂ A/B/C, generic DPA-3.1-3M (B, C arms) | same | regenerated with the frozen `run_abc_generic.py` | ≤ 1.2×10⁻¹⁴; 0 mismatches |
| SiO₂ cristobalite LOPO and Ti LOSO, R export (C-DPA, C-simple, QH-like) | DFT PBE-D3(BJ); DFT PBE, DaRUS-3582 | frozen `point_predictions.csv` (outer test) | read directly |
| H1b, 3 GPa and 10 GPa (72 models) | **ds62 assessed** | regenerated with the frozen `score_h1b.py` | selection file byte-identical; 0 mismatches |
| Hf temperature tail, validation-region event | DFT PBE, DaRUS-3582 | frozen `MULTISEED_POINT_PREDICTIONS.csv` | read directly |

**Row counts.** 498 rows, of which 146 duplicate an earlier prediction (for example, scalings that cannot change a single-block readout). 352 unique rows remain.

**They rest on only 6 physically distinct crossings:**

| Tier | System | Crossing | T_ref |
|---|---|---|---|
| DFT | SiO₂ | quartz→cristobalite | 1,543 K |
| DFT | SiO₂ | quartz→tridymite | 1,582 K |
| DFT | Hf | hcp→bcc | 1,920 K |
| DFT | Ti | hcp→bcc | 1,014 K |
| ds62 | SiO₂ | coesite→quartz | 1,193 K |
| ds62 | SiO₂ | stishovite→coesite | 1,745 K |

**Excluded, with reasons**

- **Ti/Zr/SiO₂ temperature tails.** Every reference crossing lies in the training region and no test tail contains one, so crossing metrics are N/A.
- **Hf tail.** The 1,919 K event lies in validation. It is kept as a separate stratum and flagged: for learned T2/T3 the validation region drove checkpoint selection.
- **R export `source_only_family_mean` (Ti).** The frozen file carries the documented phase-order bug, so it is excluded. `zero_delta_G` is excluded because it is degenerate by design.
- **Carbon.** See the last section; no temperature-dependent G curve exists.

## Metrics compared

| Metric | Definition |
|---|---|
| Global MAE | ΔG MAE of the event pair over the scored window |
| e₀ | Signed and absolute gap error at T_ref |
| ΔS_ref | Reference relative entropy at T_ref |
| Anchor-only linearised shift | δT = e₀ / ΔS_ref |
| Anchor + response shift (**post-hoc**) | e₀ / ΔS_pred, using the model's own slope at T_ref |
| Competitor margin | At T_ref, and (**post-hoc**) at the predicted transition T* |

**Hit/miss criteria.** Each criterion uses the same 100 K tolerance, converted to energy through each event's own ΔS_ref. No threshold was fitted to outcomes.

| Criterion | Predicts a hit when |
|---|---|
| MAE | MAE ≤ 100 K · \|ΔS_ref\| |
| Anchor-only | \|δT\| ≤ 100 K |
| Anchor-only + competitor at T_ref | the above, and no competitor intrudes at T_ref |
| AR (post-hoc) | anchor + response shift ≤ 100 K, with the correct slope sign |
| ARC (post-hoc) | AR, and no competitor intrudes at T* |

The threshold-free comparison is the Mann–Whitney AUC: the probability that a hit row has a smaller metric value than a miss row.

## Main table (unique rows; tiers kept separate)

| Tier / stratum | Rows | Distinct crossings | Hits | Acc. MAE | Acc. anchor-only | Acc. + competitor at T_ref | Acc. AR (post-hoc) | Acc. ARC (post-hoc) |
|---|---|---|---|---|---|---|---|---|
| **DFT PBE-D3(BJ)**: SiO₂ LOPO development folds | 242 | 2 | 45 | 0.719 | 0.648 | 0.702 | 0.756 | **0.880** |
| DFT PBE-D3(BJ): SiO₂ LOPO, R export | 3 | 1 | 2 | 0.667 | 1.000 | 1.000 | 1.000 | 1.000 |
| **DFT PBE**: Ti LOSO, R export | 3 | 1 | 0 | 1.000 | 0.667 | 0.667 | 1.000 | 1.000 |
| DFT PBE: Hf tail, validation event (6 of 8 rows selection-exposed) | 8 | 1 | 2 | 0.625 | 0.625 | 0.625 | 1.000 | 1.000 |
| **ds62 assessed**: H1b 3 GPa | 32 | 1 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| ds62 assessed: H1b 10 GPa | 64 | 1 | 0 | 1.000 | 0.875 | 0.938 | 1.000 | 1.000 |

**Threshold-free AUC** (SiO₂ development folds, the only stratum containing both hits and enough misses):

| Metric | AUC |
|---|---|
| Global MAE | 0.771 |
| \|e₀\| | 0.762 |
| \|anchor-only shift\| | 0.774 |

These are indistinguishable.

**Root location** (same stratum, the 151 rows with a predicted root of the event pair):
- The anchor-only linearised shift tracks the actual root shift with Spearman 0.964 and median \|difference\| 0.5 K.
- MAE ranks the absolute root shift equally well (Spearman 0.971), because these readout errors are dominated by a constant offset.

## What explains the failures (unique rows)

| Failure class | Where | Rows | Captured by |
|---|---|---|---|
| **Response (slope) error**: anchor-only shift within tolerance, but the predicted slope moves the root more than 100 K or has the wrong sign | Hf T3 (3 seeds); Ti C-simple (1); H1b 10 GPa (4); SiO₂ development (13) | 21 | AR |
| **Competitor intrudes at T_ref** | H1b 10 GPa (4); SiO₂ development (13) | 17 | competitor terms |
| **Competitor crosses between T_ref and the predicted root**: pair root within tolerance, but the adjacent reference crossing quartz→tridymite lies only about 40 K above quartz→cristobalite, and a delay beyond that spacing hands the envelope to tridymite | SiO₂ development | 59 | ARC (margin at T*) |
| **Affine three-phase degeneracy**: two-source readout predicts the held-out tridymite as an exact affine combination of the source pair, so at the source crossing it ties both (margin 0 to rounding) and the envelope switches to the wrong phase at the right temperature | SiO₂ development (tridymite targets, sources {quartz, cristobalite}) | 29 | **not captured by any error metric**; set by readout geometry |
| **Gross offset** | H1b 3 GPa (all), H1b 10 GPa, Ti, Hf T2 | many | every criterion; MAE suffices |

Excluding the 29 affine-degenerate rows, ARC classifies all 213 remaining development rows correctly; MAE gets 81.7%.

## Answers to the three questions

**1. Do decision-relevant metrics explain the observed event errors better than MAE?**

Only in the post-hoc, three-part form, and only partly.
- **Anchor-only metrics are no better than MAE.** The anchor-only gap error or shift does no better than global MAE: its AUC is equal, and as a classifier it is worse (0.648 vs 0.719 in the development folds). Errors in these readouts are mostly constant offsets, so MAE and the anchor error carry the same information.
- **What MAE cannot see is captured by two additions:** the response (slope) error at the crossing, and the competitor margin evaluated where the predicted transition lands. With both, classification rises to 0.880 overall, and to 1.000 outside the affine-degenerate class, against MAE's 0.719.
- **Caveat on the improvement.** Both additions were introduced after inspecting the failures (Hf T3 and the SiO₂ tridymite-spacing rows). The improvement is therefore descriptive, not validated.
- **Root location is not evidence for decision metrics over MAE.** The anchor-only formula locates actual roots to a median of 0.5 K, but MAE ranks root displacement equally well.

**2. Which cases do not fit?**
- **Affine-degenerate failures (29 rows).** The wrong phase wins at the correct temperature because of a three-way tie forced by the two-source affine readout. Pair-level error metrics cannot predict this; it is a readout-geometry property.
- **Slope-error failures.** They break the anchor-only linearisation (Hf T3: shift ≈ 40 K predicted, no root).
- **Uninformative strata.** H1b 3 GPa and Ti contain only misses with large errors, so every metric "explains" them trivially and they cannot discriminate between metrics. Hf's two hits are deterministic analytic controls; the six learned rows are selection-exposed.

**3. Is the case count enough for a general conclusion?**

**No.**
- The 352 unique rows rest on 6 distinct physical crossings. The 242 development rows come from just two SiO₂ crossings, re-scored across overlapping folds, pools and feature sets.
- Rows are highly correlated and are not independent samples.
- Only one stratum contains enough hits and misses to compare metrics at all.
- The post-hoc refinements were tuned to the same failures they now explain.

What the audit supports is a mechanistic description of how transitions are lost in these cases: offset, response, competitor spacing, and readout degeneracy. It does not support a general law, or the claim that decision-relevant metrics outperform MAE in general.

## Carbon (separate; no transition analysis)

- **Why no transition analysis.** The archived carbon experiment (`fes-carbon-entropy-transfer-20261029`) uses only fixed-volume harmonic S_vib(300 K) labels for graphite and diamond: 14 points, PBE-D3(BJ). Its own report states that the four-temperature QHA G table (16 points) was not used and that no Tc or phase-event score is appropriate. Carbon therefore does not enter the hit/miss table.
- **Verified entropy-interpolation comparison.** Within phase, leave-one-volume-out, MAE in µeV atom⁻¹ K⁻¹:

  | Phase | Quadratic V/V₀ baseline | DPA | Linear volume |
  |---|---|---|---|
  | Graphite | 0.00467 | 0.00743 | 0.0893 |
  | Diamond | 0.00372 | 0.00920 | 0.1129 |

- **Cross-phase transfer.** It fails with large bias (DPA 15.27 for graphite→diamond, 22.53 for diamond→graphite) and with wrong or too-shallow volume response.
- **Caveats.** These are development results on one system and seven correlated volumes per phase. Label precision below phonon-convergence uncertainty is not established.

## Files

- `DECISION_EVENT_TABLE.csv`: one row per model × pool × reference event; non-computable fields carry `N/A: <reason>`.
- `DECISION_SUMMARY.csv`: per tier and stratum.
- `curves/*_CHECK.json`: regeneration checks.
- Scripts in `scripts/`: `regen_abc_curves.py`, `regen_h1b_curves.py`, `decision_table.py` (post-hoc additions are marked in the code).
