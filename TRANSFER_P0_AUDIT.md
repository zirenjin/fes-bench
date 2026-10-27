# Transfer audit: target-phase participation, pair order and topology

## Findings

1. The held-out SiO₂ phase contributes substantially to the predicted stable envelope. Adding cristobalite changes 917/1649 decisions for each C readout; source-only accuracy is 41.9648%, compared with 97.5743% for the full envelope. The hypothesis that the high score is supplied entirely by source phases is therefore rejected by these archived predictions.
2. Both readouts nevertheless couple all three pair crossings near 1582.52 K. The reference has two distinct pair crossings and no cristobalite–tridymite root. A high envelope score does not certify the complete pairwise ordering.
3. The Ti family-control outer pair scorer used hcp−bcc coefficients with bcc−hcp labels. This was a real scoring bug. Its envelope construction used the intended phase order. Aligning the pair predictions gives MAE 15.8708 meV/atom, replacing the incorrect 12.0170 value.
4. A common numerical tie rule makes pair-sign and two-phase envelope accuracy agree in all 37 non-degenerate two-phase case/domain/method checks. The 693-point Ti target support contains three reference ties; the decision denominator is therefore 690.

## SiO₂ participation

| Method | Reference target-stable states | Predicted target-stable states | Changed decisions | Source-only accuracy | Full accuracy | Missed target-stable states |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| C-DPA | 957 | 917 | 917 | 41.9648% | 97.5743% | 40 |
| C-simple | 957 | 917 | 917 | 41.9648% | 97.5743% | 40 |

The reference makes cristobalite uniquely stable at 957 states. The models recover 917 of those, with no false target-stable states. Conditional target-stable recall is 95.8203%; the 40 missed states lie before the delayed predicted boundary. The target phase is doing real work in the envelope. Its independent relative placement and crossing topology remain inaccurate.

## Pair crossings

| Readout | Pair | Reference roots (K) | Predicted roots (K) |
| --- | --- | --- | --- |
| C-DPA | cristobalite_beta − quartz_beta | 1542.7087047758675 | 1582.5215448575386 |
| C-simple | cristobalite_beta − quartz_beta | 1542.7087047758675 | 1582.5234471589101 |
| C-DPA | cristobalite_beta − tridymite_p63mmc | None | 1582.5168893830298 |
| C-simple | cristobalite_beta − tridymite_p63mmc | None | 1582.5234471589104 |
| C-DPA | quartz_beta − tridymite_p63mmc | 1582.3215510581538 | 1582.5233631987269 |
| C-simple | quartz_beta − tridymite_p63mmc | 1582.3215510581538 | 1582.5234471589104 |

## Why the crossings are coupled in these readouts

The SiO₂ fit uses only two source phase structures, hence one source feature contrast x. The readout solves a regularized linear coefficient problem of the form

`min_W ||xᵀ W − y||²_H + α ||W||²_F`, with `α > 0`.

Decompose W into a component parallel to x and an orthogonal component. The orthogonal component cannot change the data term and increases the ridge penalty, so it vanishes in exact arithmetic. W therefore has rank at most one. For any phase contrast δz, the predicted curve is a scalar times one common temperature function. Its derivative is scaled by the same scalar. Extra frozen descriptor dimensions or a different intermediate layer do not remove this constraint when the source design and linear readout are unchanged.

The saved C-simple readout is numerically rank one. C-DPA is nearly rank one; small floating-point components remain. We report the actual residual and singular values rather than calling its roots exactly identical.

| Readout | Target/source gap factor | Max deviation from scaled source gap (meV/atom) | Second/first singular value | Root spread (K) |
| --- | ---: | ---: | ---: | ---: |
| C-DPA | 1.390582452 | 2.46061857e-05 | 1.51744708e-05 | 0.0064738157 |
| C-simple | 1.142850005 | 1.73472348e-15 | 1.97650601e-23 | 2.27373675e-13 |

The scalar projection is an after-the-fact diagnostic. It is never used to change predictions or scores. The reconstructed phase energies use one common quartz gauge at each temperature, with no phase-specific alignment. All original pair curves are checked for cycle consistency.

This isolates a constraint of the fitted mapping and the source design. It does not establish absent information in the underlying DPA representation. The current model can transfer the shared source temperature shape and a structural scale factor; it cannot independently place all four target temperature coefficients in the rank-one limit.

## Ti correction and decision convention

The saved family readout defines `delta_coeff = mean(source hcp coefficients − source bcc coefficients)`. Source pair scoring uses hcp−bcc and was consistent. Outer scoring sorts the phases as bcc,hcp but previously reused delta_coeff without negating it. The independent scorer now aligns the coefficient difference to the pair labels. All 693 affected outer pair predictions are corrected. The original files and hashes are preserved.

For decision accuracy, a reference gap with magnitude ≤1e−12 eV/atom is a numerical tie and is excluded. A predicted tie on an identifiable reference state is an abstention and is not correct. In a multiphase envelope, the two lowest energies define the corresponding tie rule. A predictor that is fully degenerate, such as ΔG=0, has no phase-decision score. This numerical tolerance is not a physical uncertainty interval. Event matching keeps the archived 100 K tolerance, ordered one-to-one matching and hit-conditioned Tc error.

| Model | Corrected MAE (meV/atom) | V2 pair sign = envelope | Identifiable states |
| --- | ---: | ---: | ---: |
| C-DPA | 9.139169 | 51.1594% | 690 |
| C-simple | 7.026457 | 25.5072% | 690 |
| QH_like_numerical | 169.708807 | 25.5072% | 690 |
| source_only_family_mean | 15.870798 | 26.5217% | 690 |
| zero_delta_G | 6.696825 | NA: fully degenerate | NA |

## Zero-gap relative skill

`skill = 1 − MAE(model) / MAE(zero gap)` on identical target-related pair states. No cross-task or cross-material total score is introduced. In Ti, C-simple and C-DPA have negative MAE skill; neither beats the zero-gap numerical control. The zero-gap control remains unable to choose a stable phase. In SiO₂, the corresponding skills are +0.7461 and +0.3175.

## Priority after the audit

- Completed: independent target participation, source-only/full envelopes, pair-order repair, consistent tie handling and zero-gap skill.
- Next: add independent admitted source contrasts, then compare the C readouts and formally trained T2/T3 on the same outer manifests. A layer-only probe on the current one-contrast SiO₂ fold cannot unlock an independent target temperature shape under the same linear ridge mapping.
- Next: use M₀/M₁/M₂ only on a compatible reference layer with an accepted physical QHA input and matched selection.
- Retained: research-review website layout, expandable compute and source evidence, and explicit NR/RF states. No new training or DFT job was run for this audit.

## Reproduction and versioning

The frozen point source is S00623 with SHA-256 `abbee6a01bab1ba3ae0ad511805bcfe0f18b89cf3ba1a9255ab4bf6e4294feea`. The family readout hash is `36fe9aac42924d2da71b75b82a475426f1a81b4c4740672edc095ab8c5492571`. `postprocess/audit_transfer_p0.py` independently builds the revised scores. The original `data/review_data.json` stays unchanged; the active score overlay is `data/transfer_score_audit.json`. Current transfer figures, the failure landscape and the fit-audit family-control curve use this V2 layer. Original values remain in comparison columns.

- [All two-phase consistency checks](data/TI_PAIR_ENVELOPE_CONSISTENCY.csv)
- [Target participation](data/SIO2_TARGET_PARTICIPATION.csv)
- [Relative skill](data/TRANSFER_RELATIVE_SKILL.csv)
- [Corrected pair metrics](data/TRANSFER_PAIR_METRICS_V2.csv)
- [Corrected envelope metrics](data/TRANSFER_ENVELOPE_METRICS_V2.csv)
- [Envelope figure](figures/sio2_envelope_audit.pdf)
