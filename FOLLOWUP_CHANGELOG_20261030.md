# First-screen redesign · 30 October 2026

The page now opens with the research question and one featured finding instead of navigation and project counts.

- **Hero.** Title, the question "When lower prediction error fails to produce better phase decisions", and a short statement of what FES-Bench measures.
- **Featured finding (development evidence, not independently confirmed).** SiO₂, held-out β-cristobalite, three source phases: DPA + volume vs DPA under equal-block feature scaling.
  - ΔS MAE 0.363 → 0.254 µeV atom⁻¹ K⁻¹ (−30%).
  - ΔG MAE 0.315 → 0.408 meV/atom (+30%).
  - Stable-phase accuracy 83.0% → 42.0%.
  - The scope note states that the trade-off does not appear under the primary coordinate scaling or with the generic DPA-3.1-3M checkpoint.
- **Featured figure.** `figures/hero/HERO_SIO2_DPA_vs_DPAV.*`, drawn by `postprocess/plot_hero_sio2.py` from archived predictions. The curves come from the post-hoc regeneration of the frozen A/B/C run, which matches its archived metrics and event lists. The script asserts the plotted accuracy, gap MAE and event identities against the frozen `PRIMARY_COMPARISON.csv` and `EVENTS.csv` before drawing. Source data: `figures/hero/HERO_SIO2_SOURCE_DATA.csv`.
- **Latest experiments.** Three entries (SiO₂ development folds, Carbon nonlinear baseline, H1b external stress test) and a link to the post-hoc decision-metric audit. Their result files are copied under `updates-20261030/` so links work offline.
- **Scope and reference status.** The previous count cards and section navigation now sit below the latest experiments. Card "2 completed outer folds" is replaced by "3 held-out transfer settings scored".
- **Dates.** The page shows "updated 30 October 2026". The snapshot line and footer distinguish the page date from the core evidence snapshot (2026-10-26) and name the later, separately versioned packages.
- **Unchanged.** Data snapshot files.
- **Checks.** Rendered with headless Edge at 1440 px, and in a 390 px frame; no horizontal page overflow. This is a local static-site source bundle and is not deployed.

## Stable-phase band comparison (added later on 30 October 2026)

- `figures/stable-phase/FIG_A_STABLE_PHASE_MAIN.*`: one protocol only. The rows are held-out β-cristobalite, three source phases and the pool {β-quartz, P6₃/mmc, β-cristobalite} at 851–2499 K: DFT reference, Volume-only, DPA, DPA + volume (primary coordinate scaling), and DPA + volume under the equal-block sensitivity scaling. Volume-only and DPA are identical under both scalings. Right-hand columns show ΔG MAE, stable-phase accuracy and transitions (hit/miss/false). Red strips mark the wrong-phase intervals.
- `figures/stable-phase/FIG_B_SOURCE_SUPERVISION.*`: two sources {β-quartz, P6₃/mmc} vs three (+ C222₁) on the same pool, A/B/C, coordinate scaling. The other two-source subsets were archived on different candidate pools and are not shown.
- **Verification.** Every accuracy, gap MAE and predicted event is asserted against the frozen `PRIMARY_COMPARISON.csv`, `ENVELOPE_METRICS.csv` and `EVENTS.csv`. Wrong-phase intervals and identities are in `STABLE_PHASE_ERROR_INTERVALS.csv`; per-temperature labels are in `STABLE_PHASE_SOURCE_DATA.csv`. Script: `postprocess/plot_stable_phase_comparison.py`.
- **Excluded by design.** Thermo4 (uses target thermal labels) and T2/T3/Phase-ID (temperature-tail protocols): different information budgets.

## Featured figures on the page (later on 30 October 2026)

- **Layout.** The featured section shows three figures directly: the cristobalite-margin figure beside the headline numbers, then the two stable-phase band figures, each with a short caption to its right (below it on narrow screens).
- **Captions.** No figure numbers or panel letters at this stage. The longer explanatory paragraph and scope note were folded into the captions.
- **Band figures.** Regenerated without in-image titles; data and assertions are unchanged.

## Homepage as one argument (evening of 30 October 2026)

- **Removed.** The "Latest experiments" three-card block, its nav link and the matching sentence in the snapshot note. Their content now appears as the figures below or in Methods & Reproducibility.
- **Order.** What happened (SiO₂ margin and stable-phase bands, unchanged) → *Why does it happen?* → *Does it generalize?* → *Is DPA actually necessary?* → collapsible **Methods & Reproducibility**. Captions are cut to 2–3 sentences; conditions, checkpoints, data tiers and scripts moved into the Methods table.
- **Why does it happen?** New `figures/error-decomposition/ERROR_DECOMPOSITION_SIO2.*` (`postprocess/plot_error_decomposition.py`). It plots the archived post-fit decomposition `e_ΔG(T) = e_ΔG(T0) − ∫ e_ΔS dt` for the cristobalite − quartz gap, DPA vs DPA + volume (equal-block). A green band shows the reference margin of cristobalite over its nearest competitor (≤ 0.317 meV/atom).
  - The script checks the archived curves against the regenerated predictions (max deviation 1.0×10⁻⁶ meV/atom).
  - It checks anchor, 2499 K, MAE and integral values against `ANCHOR_ENTROPY_PAIR_SUMMARY.csv` to 1e-9.
  - The "error inside the band" rule matches the predicted cristobalite-stable labels at all but one temperature point (DPA, at the boundary).
  - The caption states that this is an error diagnostic based on a thermodynamic identity, not an independently established causal mechanism.
  - The archived decomposition files are copied to `updates-20261030/error-decomposition/`.
- **Does it generalize?** New page version of the H1b curves, `figures/h1b/H1B_DG_CURVES_PAGE.*` (`postprocess/plot_h1b_curves.py`), drawn from the archived regenerated fits (the regeneration check must pass). It has no in-image title, and the 10 GPa axis is held on the crossing region (two-source fit marked off scale).
  - It replaces the archived figure's grey-curve note "72 total", which counted both panels, with per-panel counts (22 and 46 others; 24 and 48 fits).
  - The script asserts: no sign change at 3 GPa (24/24); 10 GPa crossings at 1354 K (DPA, both scalings) and 1368 K (DPA + volume, coordinate), i.e. 377–391 K early (difference of rounded temperatures, as in RESULT.md; unrounded 377.1–391.8 K).
  - Qualifier on the page: "External cross-tier stress test · PATTERN_ABSENT · Not an independent DFT-level confirmation."
- **Is DPA actually necessary?** New bar chart `figures/carbon-control/CARBON_NEGATIVE_CONTROL.*` (`postprocess/plot_carbon_control.py`). It shows same-phase S_vib(300 K) MAE for the quadratic V/V₀ fit vs the DPA readout, asserted against the FOLLOWUP_RESULTS.md table. It is labelled a negative control that neither tests nor supports the SiO₂ mechanism.
- **Checks.** All local links and anchors in the new block resolve. Rendered with headless Edge at 1440 px (Methods collapsed and expanded) and in a 390 px frame; no horizontal page overflow. No new computation was run.

## Figure-led homepage, second pass (late 30 October 2026)

- **Removed.**
  - The overview section "What the completed experiments establish", with its three result cards and its three navigation cards. Their unique claims moved to *Additional controls*. The "Evaluation scope and statistical units" text and the snapshot note moved into Methods & Reproducibility.
  - The old S2 section (`#sio2-intervention`, intervention ladder plus curvature follow-up). Its result is now the target-information figure; the follow-up figures and data are linked from *Additional controls* and Methods.
  - The cristobalite-margin hero figure and the two separate band figures are no longer on the page; they repeated the same stable-phase results. Their files remain and are linked from Methods.
  - `explainer.js` no longer assumes `#question-path` exists.
- **Lower Gibbs error, worse phase decisions.** New `figures/stable-phase/FIG_PHASE_DECISIONS.*`, drawn by `postprocess/plot_stable_phase_comparison.py`.
  - Layout: one DFT reference row, then two-source and three-source groups of Volume only, DPA, DPA + volume (per-dimension) and DPA + volume (equal-block). Shared temperature axis and phase colours; red wrong-phase bars; event ticks; ΔG MAE, accuracy and hit/miss/false columns. The two DPA rows are highlighted.
  - "Per-dimension" is the frozen run's "coordinate" scaling. The script asserts every row against the frozen tables. It also asserts:
    - Volume only and DPA are identical under both scalings;
    - DPA 0.853 → 0.315 meV/atom and 97.57 → 82.96%;
    - the two-source DPA event is quartz → cristobalite, while three-source DPA gives quartz → tridymite → cristobalite.
- **How error cancellation changes phase decisions.** Panels renamed to match: DPA, and DPA + volume (equal-block). The figure labels T₀, the anchor error, the signed entropy-error integral and the final gap error. The caption sits below the figure and defines the green band mathematically. It also states that ΔS is a derivative of the same Gibbs curves and that this is a decomposition of a development case, not an established causal mechanism.
- **How much target information is needed to repair a phase prediction?** New `figures/target-information/TARGET_INFORMATION_S2.*`, drawn by `postprocess/plot_target_information.py` from the archived S2 summary.
  - It is a grouped point plot, not a line, of gap MAE and transition-temperature error for C-simple and C-DPA. Practical few-shot budgets use filled markers; the oracle diagnostics use open markers on a grey background.
  - The script asserts: one label raises MAE for both readouts; two labels bring MAE below 0.1 meV/atom; the 11-point oracle reaches a transition error below 1.3 K; every condition has one identity-matched hit and no false events.
  - The page labels it "Development diagnostic — not externally validated calibration performance".
- **Does the development pattern survive external transfer?** H1b page figure redrawn at full width:
  - legends moved below the panels;
  - background fits lighter;
  - T₀ ± 2σ bands (±56 K, ±43 K) added;
  - both three-source crossings labelled (1354 K, 1368 K).
  
  The caption gives the corrected result, the ds62 tier, the lack of independent confirmation, and the three factors that cannot be attributed separately. The status cards were removed.
- **Additional controls.** Carbon bar chart, Γ-curvature query, known-phase temperature tails, and other source-pair and sensitivity experiments.
- **Methods & Reproducibility.** Moved to the top of the page-end `#methods` section, with rows for every figure plus the evaluation-scope text and the snapshot note.
- **Checks.**
  - No duplicate IDs. All local files and anchors resolve; the `#E…` register links are filled in by JS, as before.
  - Rendered at 1440 px (Methods collapsed and expanded) and in a 390 px frame; no horizontal page overflow.
  - No new computation.

## Margin panel restored in the first figure (late 30 October 2026)

- `FIG_PHASE_DECISIONS.*` now has two panels on a shared temperature axis.
  - Top: the original cristobalite-margin curves for the three-source reference, DPA and equal-block DPA + volume (G(cristobalite) − min G(quartz, tridymite)).
  - Bottom: the merged two- and three-source stable-phase bands, unchanged.
- The script additionally asserts that the DPA margin crosses zero at its 1824 K tridymite → cristobalite event, and that the equal-block margin stays positive throughout.
- The caption now describes both panels.
