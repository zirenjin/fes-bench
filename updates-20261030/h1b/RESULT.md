# H1b result: cross-tier OOD phase-transfer stress test on ds62 SiO₂

2026-10-30. Single scoring run under `H1B_PROTOCOL.md` (frozen 13:56:46) with deviations D1–D4 (`DEVIATIONS.md`, each recorded before scoring). Implementation hashes are in `IMPLEMENTATION_FREEZE.txt`; the selection hash `4cc64cd0…` was written before target scoring.

**Evidence tier.** This is the ds62 assessed-tier external stress test. It is **not** a DFT-level result and **not** an independent confirmation of H1, which remains not done.

## Verdict

| Test | Verdict | R1b gap gain | R2b decision reversal | R3b anchor-dominated |
|---|---|---|---|---|
| **Primary: 3 GPa, coesite target** | **PATTERN_ABSENT** | false | false | true |
| Secondary: 10 GPa, stishovite target | PATTERN_ABSENT | true | false | true |
| Pressure-sensitivity structures | BLOCKED (gate E2, D4); not scored | — | — | — |

The DFT-level development pattern (extra source supervision improves the gap but reverses the phase decision) **does not appear**. At 3 GPa, no readout produces any phase transition. At 10 GPa, no readout recovers the transition within tolerance (see the correction below).

> **Correction (2026-10-30, post-hoc figure review).** An earlier summary said that *no* readout predicts *any* transition and that every model predicts one phase throughout. That is true only for the 3 GPa test (24/24 fits). At 10 GPa, the three-source {q, trd, coe} readouts predict the **correct** transition (stv→coe) too early: arm B at 1,354 K (both scalings, identical model) and arm C at 1,368 K (coordinate scaling), against the reference 1,745 K. That is 377–391 K early, outside both the 100 K tolerance and the ±43 K reference uncertainty, so the frozen scoring of miss plus false event stands. The first-order anchor estimate predicts a shift of about −329 K (anchor error +20.3 meV/atom divided by ΔS_ref −61.9 µeV atom⁻¹ K⁻¹), consistent with the observed −391 K.

## What happened (primary test)

Source sets compared: two-source {quartz, tridymite} and three-source {quartz, tridymite, cristobalite}. Arm B, coordinate scaling.

| | Two-source | Three-source |
|---|---|---|
| coe−q gap MAE (meV/atom) | 6.28 | 8.71 |
| Zero-gap skill | −1.68 | −2.72 |
| Relative entropy MAE (µeV atom⁻¹ K⁻¹) | 4.76 | 6.52 |
| Stable accuracy (σ-filtered) | 46.2% (45.7%) | 53.8% (54.3%) |
| Events hit/miss/false | 0/1/0 | 0/1/0 |

- **No readout beats predicting ΔG = 0.** Across all 24 configuration × scaling × arm combinations, zero-gap skill is negative everywhere (best −0.07, equal-block C on {tridymite, cristobalite}).
- **Constant-phase predictions.** Stable accuracy is always 46.2% or 53.8%. These are the window fractions on each side of the reference crossing, so every model predicts one phase over the whole 500–2000 K window.
- **Volume, DPA and DPA+volume all fail.**
- **Anchor decomposition (three-source).** The reference crossing T₀ = 1193.3 K is shallow (ΔS_ref = −4.44 µeV atom⁻¹ K⁻¹; 2σ reference Tc uncertainty ±56 K). The anchor error at T₀ is +9.78 meV/atom, which by first-order arithmetic displaces the crossing by about 2,200 K. R3b is "true" only in this trivial sense: the anchor error is about 20 times what a 100 K tolerance allows. That is gross misplacement, not the subtle reversal seen in development.

## Secondary test (10 GPa, stishovite)

- **Three-source {q, trd, coe}:** gap MAE 22.9 meV/atom (zero-gap skill +0.29) against 716.5 for two-source {q, coe}. Stable accuracy rises from 17.0% to 73.9%, with no hit and one false event.
- **Best configuration:** some reach positive skill (best +0.60, equal-block C on {q, crst}). Still, no configuration in any arm recovers the stv→coe transition within tolerance at T₀ = 1745.4 K (±43 K).
- **Early transitions:** only the three-source B and C readouts change phase at all, in the correct direction but 377–391 K early. This is the expected first-order consequence of their anchor error, not a separate failure mode.

## Post-hoc figures (descriptive, not pre-registered)

The scoring run saved metrics but not curves. `posthoc/posthoc_figures.py` regenerates the deterministic fits by importing the frozen `score_h1b.py` unchanged. It proceeds only if the regenerated `FROZEN_SELECTION.csv` is byte-identical to the frozen one (SHA-256 `4cc64cd0…`, confirmed) and the gap MAEs match (maximum difference 0.0); see `posthoc/REGENERATION_CHECK.json`.

- `posthoc/FIG_H1B_DG_CURVES.pdf`: reference ±2σ against predicted ΔG at 3 and 10 GPa. The 3 GPa reference shows a kink near 1,550 K, from ds62's Landau α–β quartz term, which the smooth four-term basis cannot represent. It lies above T₀ and does not explain the missed sign change.
- `posthoc/FIG_H1B_OOD.pdf`: structure-level and atom-level OOD diagnostics.

## Out-of-distribution status (label-free, recorded before scoring)

| Target | Nearest source / median source–source distance | Atomic NN distance q50 / q95 (source leave-one-out q95 = 0.25) | Affine-hull weights (q, trd, crst) |
|---|---|---|---|
| Coesite | 1.86 | 0.15 / 0.25 | 1.14, 0.72, −0.86 |
| Stishovite | 10.98 | 0.80 / 1.07 | 2.51, −4.00, 2.50 |

Coesite sits at the edge of source coverage. Stishovite's six-coordinated environments lie far outside it.

## Interpretation, bounded

- **What it shows.** Under assessed-tier labels and for phases outside the task checkpoint's training distribution, the frozen linear readout does not transfer. It does no better than ΔG = 0 for coesite and recovers neither transition within tolerance (at 10 GPa the best readouts cross in the right direction but 377–391 K early).
- **What it cannot show.** Which of three confounded factors drives the failure:
  - representation extrapolation (out of distribution);
  - the cross-tier label change (ds62 versus PBE-D3);
  - the declared state mismatch (near-0 GPa, high-T, compressed source representatives used for 3 and 10 GPa labels).

  No further runs were made to separate them; the protocol allows one run.
- **For H1.** This result neither supports nor refutes the DFT-level mechanism. It shows the mechanism is not observable in this stress setting, because the readout fails at a coarser level first. Independent DFT-level confirmation of H1 remains open.

## Gates and integrity

| Gate | Result |
|---|---|
| E1, source re-embedding vs archive (patched checkpoint) | ≤ 2.4×10⁻⁷, pass |
| E2, primary targets, supercell invariance (D1 rule) | coesite 6.6×10⁻⁸; stishovite 2.1×10⁻⁷; pass |
| E2, sensitivity set | 9.26 GPa stishovite supercell exceeded memory; set blocked (D4) |
| R1, ds62 reproduction (stored precision, D2) | ≤ 0.000499 meV/atom, pass |
| R2, unique T₀ | 1193.3 K (3 GPa), 1745.4 K (10 GPa), pass |
| Budget | 0 SU; no DFT, TI or MD; local CPU only |

## Files

- `scoring/primary/`: `VERDICT.json`, `FROZEN_SELECTION.csv`, `SELECTION_FREEZE.json`, `PAIR_METRICS.csv`, `ENVELOPE_METRICS.csv`, `EVENTS.csv`, `ANCHOR_DECOMPOSITION.csv`
- `scoring/sensitivity/BLOCKED.json`
- `stage0/STAGE0_EMBEDDING_AND_OOD.json`
- `reference/GATES_R.json` and `reference/REFERENCE_DS62.npz`
- `audit/` (structure and topology audits)
- `H1B_PROTOCOL.md`, `DEVIATIONS.md`, `IMPLEMENTATION_FREEZE.txt`
