# Relative-derivative provenance audit

## Finding

The two Hf pair labels denote the same oriented difference, bcc minus hcp. They are complementary model rows within `hf_60_20_20_batch2_3000`, not two independent experiments. The derivative builder creates `bcc−hcp` by subtracting phase-resolved curves and retains `bcc_minus_hcp` for direct pair exports. The old plotter grouped raw strings and split these models into two views.

- `bcc−hcp`: T2 and T3, each summarized over three seeds.
- `bcc_minus_hcp`: four analytic controls, each with one archived seed identifier. These deterministic controls are not extra stochastic replicates.
- All three domains retain all three smoothing windows. No numerical metric has been changed or averaged across models.
- The long-run batch-8 Hf protocol remains separate from the batch-2 protocol.

The directly exported Thermo3 and Thermo4 references match the phase-subtracted T2 seed-11 reference at every one of 1,252 (domain, temperature) points per model, to less than 1e-12 eV/atom. This confirms orientation, evaluation support and reference agreement. The builder assigns the same split protocol; it does not establish identical training budgets across the learned and analytic methods.

## Counting hierarchy

The 405 summarized rows contain four systems, six system–directional-pair combinations and seven protocol–system–pair cases. Each case has Train, Validation and Test views: 21 views in total. The seven cases are not seven independent materials. Pairs within SiO₂ share phases. Hf protocols share a material. Domains and smoothing windows are repeated analyses, not independent material replicates. Seed counts remain per model and must not be summed over domains or windows.

No duplicate keys were found after canonicalizing `(protocol, system, domain, directional pair, model, window)`. Thus there are no metric rows to discard. The error was in display partitioning. The upstream summary averages only within model/domain/pair/window over seed rows; it does not pool material or domain scores. The current failure-landscape builder does not consume this derivative table. No overall derivative score or independent-material count was computed by this portal.

## Reproducibility

`data/DERIVATIVE_CANONICAL_ROWS.csv` preserves each raw pair label, source ID, source hash and original CSV line. `data/DERIVATIVE_PROVENANCE_AUDIT.json` records counts and pointwise checks. The original source table and evidence snapshot are unchanged. The pair aliases are merged only within the original protocol and domain. Other protocols remain distinct.

The upstream `derive_relative_shape_diagnostics.py` and `summarize_shape_sensitivity.py` explain the two construction paths. The portal plot script now normalizes the verified alias before grouping and rejects duplicate model/window keys.
