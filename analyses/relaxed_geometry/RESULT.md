# Relaxation channel: result (head `MP_traj_v024_alldata_mixu`)

Protocol: `PROTOCOL.md`, SHA-256 463ff3b4…, frozen 2026-11-01T23:29Z. All analyses are pre-specified and descriptive. Results are in `results/RELAX_RESULTS.json`; per-entry files are in `results/relax_npz.tgz`.

## Exclusions

These are failure modes of the self-consistent harmonic MLIP pipeline, reported and not dropped silently.

| Reason | Entries |
|---|---:|
| Dynamically unstable at the MLIP geometry (minimum frequency < −0.05 THz) | **44 of 521 (8.4%)** |
| Space group changed despite FixSymmetry | 6 |
| Failed (GPU out of memory, mp-730190, same entry as before) | 1 |
| Not converged | 0 |

All 44 unstable entries are dynamically stable at the DFT geometry by selection. On these, the MLIP pipeline cannot produce a harmonic free energy at all.

**Used:** 470 entries, 295 pairs, 203 formulas.

## Geometry

- **Volume change.** Median +0.29%, with the 10th–90th percentile range −0.54% to +1.48%.
- **Static relaxation energy.** Median −0.29 meV/atom.

The geometry change is small, but its effect on the vibrational free energy is not (below).

## R1–R4 (formula-weighted, 95% cluster-bootstrap CI)

| T | RMS e_fixed | RMS δ_relax | static part | vibrational part | RMS e_pipeline | Correlation of δ_relax with e_fixed | Wrong ordering, fixed → pipeline |
|---|---|---|---|---|---|---|---|
| 300 K | 6.54 | 1.80 | 0.82 | 1.69 | 6.58 | −0.11 | 11% → 11% |
| 1000 K | 8.48 | 5.52 | 0.82 | 5.54 | 7.74 | −0.41 | 12% → 15% |
| 1500 K | 11.30 | 8.34 | 0.82 | 8.38 | 9.61 | −0.51 | 11% → 11% |

RMS values are in meV/atom.

- **Size.** At 1000–1500 K the RMS of the relaxation term reaches 65–75% of the RMS of the fixed-geometry error. This is a magnitude comparison, not a share of the total: the two terms partly cancel, so they do not add as an error budget. The term is almost entirely vibrational: a small geometry change moves the phonons.
- **Common mode.** The relaxation-term cancellation fraction is 0.56–0.68, weaker than the fixed-geometry vibrational cancellation (about 0.70).
- **Interaction.** The term is anticorrelated with e_fixed (−0.41 / −0.51), so the pipeline RMS is somewhat below the fixed-geometry RMS at high T. But |e_pipeline| < |e_fixed| holds in only 44–53% of pairs, and decisions do not improve. Coverage at α = 0.05 is 0.44 → 0.41 at 1000 K.

## R5: label-free warning

- The pair difference in static relaxation energy correlates with |e_pipeline|: Spearman 0.26 / 0.34 / 0.33. This is computable from the DFT structure plus the MLIP alone, with no phonon labels. It is a candidate risk feature only; the held-out test below does not show a useful gain.
- The R₀ magnitude does not correlate (−0.07 to +0.02).

## Reading

- Fixed-geometry tests do not report the error of the actual MLIP workflow.
  - About 8% of DFT-stable phases fail outright at the MLIP geometry.
  - The relaxation term is as large as the fixed-geometry thermal error at high T. It partly offsets that error on average, but not reliably for a given pair.
- Earlier fixed-geometry conclusions (temperature-dependent dominance, abstention coverage) should be re-stated for the pipeline before they are used as benchmark claims.

## Coverage over all pairs and conditional-subset bias (post-hoc; `scripts/relax_coverage.py`, `results/RELAX_COVERAGE.json`)

All numbers here are pair-weighted.

- **Evaluable pairs.** 295 of 351 pairs (84%). 36 of 239 formulas lose all their pairs.
- **The successful subset is easier.**

  | T | Fixed-geometry RMS, all pairs | Evaluable pairs | Failed pairs |
  |---|---:|---:|---:|
  | 300 K | 7.17 | 6.07 | 11.32 |
  | 1000 K | 9.71 | 7.70 | 16.70 |
  | 1500 K | 12.91 | 10.12 | 22.48 |

  Values are meV/atom. Median |e| differs only modestly (3.3 against 4.0 meV/atom at 1000 K), so the difference is in the tail. Pipeline RMS values in the table above are therefore conditional on the easier subset. They should only be compared with fixed-geometry values on the same subset, as R1 does.
- **Decisions over all non-degenerate pairs** (284–297), with failed pairs as non-decisions:
  - the pipeline decides 84–85% of pairs;
  - wrong decisions are 10.9 / 12.7 / 10.4% of all pairs at 300 / 1000 / 1500 K;
  - the worst-case bound, counting failures as wrong, is 27 / 28 / 26%.
- **Wrong ordering is not enriched among the failed pairs** at fixed geometry: 13 / 7 / 13% against 13 / 11 / 10% for evaluable pairs.
- **R5 held-out check** (leave-one-formula-out conformal, α = 0.05, evaluable pairs). Normalising the abstention threshold by the static-relaxation feature lowers coverage at every T:

  | T | Coverage (plain → normalised) | Selective risk (plain → normalised) |
  |---|---|---|
  | 300 K | 0.44 → 0.39 | 1.0% → 1.1% |
  | 1000 K | 0.40 → 0.37 | 2.1% → 1.1% |
  | 1500 K | 0.34 → 0.31 | 0 → 0 |

  It does not improve the risk–coverage trade-off. It is not yet a validated warning signal, and the R₀ magnitude is not either.
- **Scope.** This is the MLIP-relaxed harmonic pipeline. It is not a finite-temperature phase diagram with thermal expansion or anharmonicity.
