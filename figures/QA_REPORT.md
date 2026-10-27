# Scientific-figure revision: QA record

The original 71 scientific figure exports passed the PDF collision and glyph checks. 55 multi-panel figures passed the measured 1.5 pt alignment gate; 16 single-panel figures were not applicable. Minimum effective PDF text size: 8 pt. All 71 include PDF, editable-text SVG, 300 dpi PNG, and a figure-specific source CSV.

## Visual and scientific review

The complete figure series and all panels were inspected using Python-generated contact sheets, with individual full-size checks of the evidence matrix, matched-tail overview, paired-transfer plot, carbon P–T grid, Hf reconstruction, derivative legend, and manuscript counts. Page screenshots confirm the new figures are connected to index.html.

- Evidence links are binary relationships, not evidence-strength scores.
- Reference matrix values preserve separate units and admission layers. Source-pair candidates are not automatically accepted truth.
- Tail dots show every available row; unavailable methods remain blank. Circle: seed 11 or deterministic control; diamond: seed 23; triangle: seed 37.
- Outer-fold headline plot uses only target-related pairs. Pair sign accuracy and stable-envelope accuracy stay separate. Event counts and undefined Tc errors remain in a table.
- Carbon bounds are deterministic numerical sensitivity, not a confidence interval. Unresolved sign is outlined explicitly.
- Ti remains a single-phase comparison layer. All 155 T/P states are retained.
- Historical Hf retains all 1252 temperature rows and its original protocol. It is not pooled with later training recipes.
- Derivative figures retain all 405 archived rows across 24 exact selector partitions; ranges are archived seed minima/maxima. Distinct 3/4-term forms also use dashed/solid lines and square/circle marks.
- Compute counts are ledger rows, not unique jobs or admitted labels. Manuscript counts are placeholder occurrences, not unfinished experiment counts.

## Resolved collision findings

Two initial bar-plot grid lines crossed numeric labels. The unnecessary grid lines were removed, the PDFs were re-exported, and the final audits passed. No final collision warning remains.

## Reviewed source-preflight advisories

The source heuristic reports four non-blocking advisories: PNG rather than TIFF, 300 rather than 600 dpi, width expressed as 183/25.4 rather than a literal inch value, and vertical evidence-ID ticks. These are intentional for an offline evidence page: PDF/SVG are the scientific vector deliverables; PNG is a screen preview; every recorded final width is 183 mm; rotated IDs were inspected and collision-audited. This package does not claim journal production submission clearance.

## Per-figure checks

| Figure | Rows | Alignment | Collision | Minimum PDF text (pt) |
|---|---:|---|---|---:|
| evidence_map | 19 | NOT APPLICABLE | PASS | 8.0 |
| reference_inventory | 8 | NOT APPLICABLE | PASS | 8.0 |
| coverage_targets | 4 | NOT APPLICABLE | PASS | 8.0 |
| tail_cross_system_60_20_20_batch8_3000_ti_test | 11 | PASS | PASS | 8.0 |
| tail_cross_system_60_20_20_batch8_3000_ti_train | 11 | PASS | PASS | 8.0 |
| tail_cross_system_60_20_20_batch8_3000_ti_validation | 11 | PASS | PASS | 8.0 |
| tail_cross_system_60_20_20_batch8_3000_zr_test | 11 | PASS | PASS | 8.0 |
| tail_cross_system_60_20_20_batch8_3000_zr_train | 11 | PASS | PASS | 8.0 |
| tail_cross_system_60_20_20_batch8_3000_zr_validation | 11 | PASS | PASS | 8.0 |
| tail_cross_system_60_20_20_batch8_3000_sio2_test | 6 | PASS | PASS | 8.0 |
| tail_cross_system_60_20_20_batch8_3000_sio2_train | 6 | PASS | PASS | 8.0 |
| tail_cross_system_60_20_20_batch8_3000_sio2_validation | 6 | PASS | PASS | 8.0 |
| tail_hf_60_20_20_batch2_3000_hf_train | 8 | PASS | PASS | 8.0 |
| tail_hf_60_20_20_batch2_3000_hf_validation | 8 | PASS | PASS | 8.0 |
| tail_hf_60_20_20_batch2_3000_hf_test | 8 | PASS | PASS | 8.0 |
| tail_hf_60_20_20_batch8_10000_selected8000_hf_train | 1 | PASS | PASS | 8.0 |
| tail_hf_60_20_20_batch8_10000_selected8000_hf_validation | 1 | PASS | PASS | 8.0 |
| tail_hf_60_20_20_batch8_10000_selected8000_hf_test | 1 | PASS | PASS | 8.0 |
| tail_overview | 28 | PASS | PASS | 8.0 |
| transfer_overview | 6 | PASS | PASS | 8.0 |
| transfer_sio2_lopo_cristobalite_beta_fit | 4 | PASS | PASS | 8.0 |
| transfer_sio2_lopo_cristobalite_beta_source_validation | 4 | PASS | PASS | 8.0 |
| transfer_sio2_lopo_cristobalite_beta_source_id_test | 4 | PASS | PASS | 8.0 |
| transfer_sio2_lopo_cristobalite_beta_outer_test | 12 | PASS | PASS | 8.0 |
| transfer_metals_loso_ti_fit | 10 | PASS | PASS | 8.0 |
| transfer_metals_loso_ti_source_validation | 10 | PASS | PASS | 8.0 |
| transfer_metals_loso_ti_source_id_test | 10 | PASS | PASS | 8.0 |
| transfer_metals_loso_ti_outer_test | 5 | PASS | PASS | 8.0 |
| carbon_pt_grid | 16 | NOT APPLICABLE | PASS | 8.0 |
| carbon_0_0_1 | 4 | NOT APPLICABLE | PASS | 8.0 |
| carbon_0_0_0 | 4 | NOT APPLICABLE | PASS | 8.0 |
| carbon_5_0_1 | 4 | NOT APPLICABLE | PASS | 8.0 |
| carbon_5_0_0 | 4 | NOT APPLICABLE | PASS | 8.0 |
| carbon_10_0_1 | 4 | NOT APPLICABLE | PASS | 8.0 |
| carbon_10_0_0 | 4 | NOT APPLICABLE | PASS | 8.0 |
| carbon_15_0_1 | 4 | NOT APPLICABLE | PASS | 8.0 |
| carbon_15_0_0 | 4 | NOT APPLICABLE | PASS | 8.0 |
| ti_all | 155 | PASS | PASS | 8.0 |
| ti_0_0 | 31 | PASS | PASS | 8.0 |
| ti_1_0 | 31 | PASS | PASS | 8.0 |
| ti_2_0 | 31 | PASS | PASS | 8.0 |
| ti_3_0 | 31 | PASS | PASS | 8.0 |
| ti_4_0 | 31 | PASS | PASS | 8.0 |
| hf_prediction | 1252 | NOT APPLICABLE | PASS | 8.0 |
| hf_components | 1252 | NOT APPLICABLE | PASS | 8.0 |
| shape_01 | 24 | PASS | PASS | 8.0 |
| shape_02 | 24 | PASS | PASS | 8.0 |
| shape_03 | 24 | PASS | PASS | 8.0 |
| shape_04 | 24 | PASS | PASS | 8.0 |
| shape_05 | 24 | PASS | PASS | 8.0 |
| shape_06 | 24 | PASS | PASS | 8.0 |
| shape_07 | 24 | PASS | PASS | 8.0 |
| shape_08 | 24 | PASS | PASS | 8.0 |
| shape_09 | 24 | PASS | PASS | 8.0 |
| shape_10 | 21 | PASS | PASS | 8.0 |
| shape_11 | 21 | PASS | PASS | 8.0 |
| shape_12 | 21 | PASS | PASS | 8.0 |
| shape_13 | 21 | PASS | PASS | 8.0 |
| shape_14 | 21 | PASS | PASS | 8.0 |
| shape_15 | 21 | PASS | PASS | 8.0 |
| shape_16 | 6 | PASS | PASS | 8.0 |
| shape_17 | 6 | PASS | PASS | 8.0 |
| shape_18 | 6 | PASS | PASS | 8.0 |
| shape_19 | 12 | PASS | PASS | 8.0 |
| shape_20 | 12 | PASS | PASS | 8.0 |
| shape_21 | 12 | PASS | PASS | 8.0 |
| shape_22 | 3 | PASS | PASS | 8.0 |
| shape_23 | 3 | PASS | PASS | 8.0 |
| shape_24 | 3 | PASS | PASS | 8.0 |
| compute_outcomes | 253 | NOT APPLICABLE | PASS | 8.0 |
| manuscript_gaps | 11 | NOT APPLICABLE | PASS | 8.0 |






## English first-principles additions

Two additional Python figures use archived observations. `explainer_fit_controls` keeps rounded report values and marks Hf as a separate protocol. `explainer_decision_example` retains all four selected SiO₂ method/seed observations. Neither figure invents data or a confidence interval. Both passed the final PDF collision and 5 pt glyph-floor checks; minimum observed size is 8 pt. The two-panel decision figure passed the measured alignment gate.

The English page retains original metric values. Three Hf derivative plots now combine complementary learned and analytic model rows after the pair-alias audit. The visible scientific captions are English. Original documents remain in their original language. The browser checks found no JavaScript error, missing request or unexpected Chinese narrative text.

## Failure-landscape additions

Three new variants map 10 explicit cases × 9 methods. Numeric cells keep reported seed ranges; NR/NA/RF/NE are distinct data states. Colour uses the worst reported seed. Ranges in the MAE panel use two lines to keep labels inside cells. Hatch fills were removed after the initial PDF geometry check; final exports passed collision, alignment and glyph checks. Minimum landscape text size is 7 pt. No new model calculation or causal ablation was run.

## Derivative provenance correction

The Hf batch-2 pair aliases now share three domain figures. All six complementary models and all windows remain; no metric values are recomputed. Final PDFs passed collision and glyph checks (minimum 8 pt) and measured panel alignment. The active derivative selector contains seven protocol cases with three domain views each. See DERIVATIVE_PROVENANCE_AUDIT.md.

## Hf phase-order audit

The four-panel Hf audit retains all 1,252 temperature pairs and independently compares phase exports, pair differences, sign-counterfactual residuals and static/vibrational terms. It passed measured panel alignment, final PDF collision and glyph checks (minimum 7 pt). This verifies arithmetic and provenance, not physical QHA validity.

## Optimization Sanity Check additions

Seven new source-fit/test distribution figures retain all matched paired points and separate seeds. The horizontal limits match between panels. No new fitting, target alignment, convergence threshold or layer probe is introduced. All final PDFs passed collision, glyph-floor and measured panel-alignment checks. See OPTIMIZATION_AND_PROBING.md for source roles and matrix ranks.

## Transfer score revision V2

The target-phase envelope audit and all nine transfer figure variants passed final PDF collision, glyph and measured alignment gates. Current scores align the Ti family coefficient order and exclude reference ties from decision denominators. All 37 two-phase case/domain/model checks have identical pair-sign and envelope accuracy. Original frozen evidence remains unchanged; derived V2 tables retain original-value columns.
