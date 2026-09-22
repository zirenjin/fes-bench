# E1 skill scores against the training-mean constant ΔG floor

For each SiO2 pair, `constant_delta_g` is the mean reference ΔG over frozen temp-extrap training indices 0–591 (851–1442 K), then evaluated on the full 851–2499 K grid. `skill = 1 − MAE_model / MAE_constant`. The all-pair aggregate applies that formula to the mean pair MAE; it is not the arithmetic mean of pairwise skills.

| Pair | training mean ΔG (eV/atom) | constant ΔG MAE (eV/atom) |
|---|---:|---:|
| crist–quartz | 0.00287686 | 0.00434745 |
| crist–trid | -0.00019687 | 0.00007268 |
| quartz–trid | -0.00307373 | 0.00427486 |

## Per-checkpoint scores

`crist–trid` is a no-reference-crossing pair with a 0.0727 meV/atom floor; it is shown, but its scale makes individual-pair skills highly unstable. **Every aggregate score below is negative**, meaning the checkpoint has larger full-grid ΔG MAE than the corresponding constant floor.

| basis | seed | all pairs: model / constant MAE (meV), skill | excluding crist–trid: model / constant MAE (meV), skill | crist–quartz skill | crist–trid skill | quartz–trid skill |
|---|---:|---|---|---:|---:|---:|
| polynomial | 11 | 12.789 / 2.898, **-3.412** | 17.540 / 4.311, **-3.068** | -2.673 | -44.212 | -3.471 |
| polynomial | 23 | 12.843 / 2.898, **-3.431** | 17.625 / 4.311, **-3.088** | -2.693 | -44.120 | -3.490 |
| polynomial | 37 | 12.868 / 2.898, **-3.440** | 17.661 / 4.311, **-3.097** | -2.701 | -44.163 | -3.499 |
| polynomial | 51 | 12.795 / 2.898, **-3.415** | 17.549 / 4.311, **-3.071** | -2.675 | -44.226 | -3.473 |
| polynomial | 67 | 12.835 / 2.898, **-3.428** | 17.607 / 4.311, **-3.084** | -2.688 | -44.277 | -3.487 |
| tlog_polynomial | 11 | 6.525 / 2.898, **-1.251** | 8.692 / 4.311, **-1.016** | -0.754 | -29.140 | -1.283 |
| tlog_polynomial | 23 | 6.508 / 2.898, **-1.245** | 8.667 / 4.311, **-1.010** | -0.749 | -29.104 | -1.277 |
| tlog_polynomial | 37 | 6.516 / 2.898, **-1.248** | 8.680 / 4.311, **-1.013** | -0.752 | -29.094 | -1.280 |
| tlog_polynomial | 51 | 6.517 / 2.898, **-1.249** | 8.681 / 4.311, **-1.014** | -0.752 | -29.121 | -1.280 |
| tlog_polynomial | 67 | 6.507 / 2.898, **-1.245** | 8.666 / 4.311, **-1.010** | -0.748 | -29.111 | -1.277 |

## Aggregation sensitivity

| basis | full-grid model MAE, all pairs (meV/atom) | excluding crist–trid (meV/atom) |
|---|---:|---:|
| polynomial | 12.8260095 ± 0.0300687 | 17.5963761 ± 0.0459907 |
| tlog_polynomial | 6.5145488 ± 0.0067310 | 8.6774183 ± 0.0097252 |

Recommendation: report the no-reference-crossing-pair aggregate as the primary ΔG MAE, with the all-pair value in a footnote/supplement. The crist–trid pair remains useful for false-crossing diagnostics, but its 0.059 meV/atom reference standard deviation (and 0.073 meV/atom constant MAE) otherwise gives it disproportionate leverage.

## Deviations from design

No new model was trained. The E1 historical checkpoints are evaluated on their original full grid; the only fitted quantity is the explicitly requested pairwise training-mean constant.
