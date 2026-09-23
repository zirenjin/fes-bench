# ΔG MAE convention audit

## Resolution

The reported **2.56 meV/atom is not an in-domain aggregate**. It is the single `quartz_beta__tridymite_c2221` cell (`0.002562650 eV/atom`) in the historical four-phase file `fes_static_only_4phase_extensive_gauge_c_a100/all_pair_metrics.json`, evaluated on 1649 full-grid points. It therefore cannot be compared to the current three-phase E1 aggregate as an ‘ID MAE’. The apparent contradiction is a mislabeled/detached table cell plus changed phase set/checkpoint, not a train-vs-test, pair-vs-phase, or grid convention difference.

## Matched re-evaluation

Current E1 uses the same ten historic four-phase checkpoints but restricts evaluation to quartz, cristobalite, and P6₃/mmc tridymite, on the original 851–2499 K / 1649-point grid. Its reported aggregate is the equal-weight mean of the three pairwise full-grid MAEs.

| basis | full grid (meV/atom) | frozen train window 851–1442 K (meV/atom) | test window 1443–2499 K (meV/atom) |
|---|---:|---:|---:|
| polynomial | 12.8260095 ± 0.0300687 | 18.1287186 ± 0.0417663 | 9.8560911 ± 0.0235818 |
| tlog_polynomial | 6.5145488 ± 0.0067310 | 9.8282543 ± 0.0108535 | 4.6586229 ± 0.0044701 |

These recomputed full-grid values reproduce E1 (polynomial ≈12.826 meV/atom; T-log ≈6.515 meV/atom). Neither matched train-window value (≈18.13 and ≈9.83 meV/atom) is 2.56 meV/atom. ΔG cancels the system gauge, so the formal-train versus heldout calibration convention also cannot change these ΔG MAEs.

## Pair-set effect

The historic 2.56265 meV value belongs to a pair containing C2221, whereas E1 deliberately excludes C2221 because its representative structure was not comparable. It also uses a different historical gauge-C run, rather than either of the polynomial/T-log checkpoint families. Thus the values answer different experimental questions and must not appear together in a common ΔG-MAE comparison.

## Required reporting convention

Use: ‘equal-weight mean of pair-level ΔG MAE, full 851–2499 K grid, three-phase subset, with no-crossing crist–trid reported separately (and excluded from the primary aggregate).’ Never label the isolated C2221-pair number as ID aggregate MAE.

## Deviations from design

No new model was trained; this is a read-only re-evaluation of frozen checkpoint curves and a read-only inspection of the archived historical metrics file.
