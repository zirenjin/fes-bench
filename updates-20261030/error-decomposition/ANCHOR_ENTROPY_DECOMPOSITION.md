# DPA → DPA+volume: anchor–integrated-entropy error decomposition

## Question and protocol

This is post-fit analysis of the frozen three-source cristobalite development case. No model, feature scaling, ridge value, or target fit was changed. It compares B (DPA) with the predeclared C primary arm (DPA+volume, equal total group weight) using saved target phase curves.

The common anchor temperature is the already specified reference quartz→cristobalite boundary, `T0 = 1542.7087048 K`. For each directed phase gap, `eG(T) = predicted ΔG(T) − reference ΔG(T)` and `δΔS(T) = predicted ΔS(T) − reference-derived ΔS(T)`. We evaluate

`eG(T) = eG(T0) − ∫[T0,T] δΔS(t) dt`.

The entropy term is integrated in µeV/(atom K) and converted to meV/atom. Local entropy errors use the same nearest-21-temperature Thermo4 derivative convention as the primary analysis. Numerical decomposition closure is below `8e-8 meV/atom` for all pairs, so the saved curves are reproduced to numerical precision. The reference entropy remains a derivative of supplied G, not an independent phonon result.

## Pairwise summary

Gap directions are explicitly shown. Anchor and ΔG errors are in meV/atom; ΔS errors are in µeV/(atom K). MAEs and biases cover 851–2499 K.

| Arm | Directed gap | Anchor error at T0 | ΔG MAE | ΔG bias | ΔS MAE | ΔS bias | ΔG error at 2499 K | Integrated entropy contribution at 2499 K |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| B DPA | cristobalite − quartz | +0.41394 | 0.31488 | +0.30321 | 0.36339 | +0.36339 | −0.12603 | −0.53997 |
| C equal groups | cristobalite − quartz | +0.44336 | 0.40817 | +0.40817 | 0.25127 | −0.15151 | +0.41810 | −0.02526 |
| B DPA | tridymite (P63/mmc) − quartz | +0.00147 | 0.00115 | +0.00115 | 0.01130 | +0.01130 | +0.00141 | +0.00006 |
| C equal groups | tridymite (P63/mmc) − quartz | +0.00149 | 0.00116 | +0.00116 | 0.01130 | +0.01120 | +0.00147 | +0.00010 |
| B DPA | cristobalite − tridymite (P63/mmc) | +0.41247 | 0.31511 | +0.30321 | 0.36326 | +0.36210 | −0.12922 | −0.54169 |
| C equal groups | cristobalite − tridymite (P63/mmc) | +0.44186 | 0.40817 | +0.40817 | 0.25731 | −0.15270 | +0.41485 | −0.02702 |

The integrated contribution column is `−∫δΔS`; it adds to the anchor error to reconstruct the total ΔG error. For the cristobalite−quartz gap, C’s anchor is 0.02941 meV/atom worse than B. More consequentially, by 2499 K B’s integrated entropy term is −0.540 meV/atom, partially cancelling its positive anchor error; C’s term is only −0.025 meV/atom, so its gap error remains +0.418 meV/atom. At 851 K the paired effect goes the other way: C’s integrated contribution is 0.334 meV/atom more negative than B’s, leaving C with a smaller total gap error there. Thus the ΔG MAE tradeoff reflects both the 0.029 meV/atom anchor change and the temperature-dependent sign/accumulation of entropy error, not a uniform loss across the interval.

## Why the competing phase remains in the envelope

At `T0`, the reference quartz and cristobalite are tied. The reference tridymite−quartz gap is `+0.29402 meV/atom`. Both readouts recover that source–source gap closely: B predicts `+0.29549` (error `+0.00147`) and C predicts `+0.29551` (error `+0.00149`).

The target–competitor gap is the failure point. The reference cristobalite−tridymite gap at `T0` is `−0.29402 meV/atom`, meaning cristobalite is below tridymite. B predicts `+0.11845` and C predicts `+0.14784 meV/atom`, reversing the ordering by `+0.41247` and `+0.44186 meV/atom`, respectively. C therefore worsens the target’s anchor relative to the competing phase by `0.02939 meV/atom`; the much more accurate quartz–tridymite gap does not fix this target-specific offset. The model can consequently retain the source-supported quartz→tridymite crossing, while cristobalite enters too late or not as the correct stable-envelope transition.

At `T0` quartz is still the lowest phase in both predicted envelopes; the pairwise ordering error is latent there and affects the sequence as temperature rises. This is why the decomposition considers all three pair gaps instead of attributing the envelope only to the target–quartz curve.

## B→C paired changes

For cristobalite−quartz, C−B changes the anchor error by `+0.02941 meV/atom`, the full-support ΔG MAE by `+0.09329 meV/atom`, and ΔS MAE by `−0.11212 µeV/(atom K)`. The integrated entropy contribution changes by `−0.33386 meV/atom` at 851 K and `+0.51471 meV/atom` at 2499 K. The sign change across the interval explains why the lower ΔS MAE does not imply a better integrated gap over the full support.

For cristobalite−tridymite, corresponding changes are nearly the same: anchor `+0.02939 meV/atom`, ΔG MAE `+0.09306 meV/atom`, ΔS MAE `−0.10595 µeV/(atom K)`, and high-temperature integrated contribution `+0.51467 meV/atom`. For tridymite−quartz, all changes are negligible. This localizes the B→C difference to the target phase’s relative placement and temperature response, not to the source–source comparison.

## Interpretation limits

This decomposition clarifies how the selected readout changes the developed cristobalite curves; it does not show that volume supplies independent thermodynamic information or that DPA lacks it. The effect is conditional on the frozen source phases, equal-group scaling and linear ridge readout. Every SiO₂ fold previously inspected remains development evidence; this analysis is not confirmation.

## Files

- `ANCHOR_ENTROPY_PAIR_SUMMARY.csv`: per-arm and per-pair anchor, error, slope and endpoint summaries.
- `ANCHOR_ENTROPY_PAIRED_B_TO_C.csv`: direct C−B contrasts by pair.
- `ANCHOR_ENTROPY_ERROR_CURVES.csv`: full-temperature reconstruction and closure residuals.
- `decompose_anchor_entropy.py`: reproducible post-fit analysis.
