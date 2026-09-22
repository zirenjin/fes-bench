# E1 — SiO2 Tc re-evaluation on the three-phase subset

Historic four-phase checkpoints are evaluated only on quartz_beta, cristobalite_beta, and tridymite_p63mmc. C2221 is excluded. Every model was evaluated on the original 851–2499 K, 1 K grid; the predictor-level gauge uses the frozen temp-extrap train window (851–1442 K).

Each seed cell is `predicted Tc (Tc error; false/missed)`. The reference slope is `|dΔG/dT|`; the slope conversion is reported as the evaluator's crossing-local `|ΔG error|/|dΔG/dT|`.

## Table 6-style per-seed Tc results

| predictor | pair | reference Tc (K) | |dΔG/dT| (eV/atom/K) | seed 11 | seed 23 | seed 37 | seed 51 | seed 67 | mean Tc error ± seed std (K) | ΔG→Tc (K) |
|---|---|---:|---:|---|---|---|---|---|---:|---:|
| polynomial | sio2:cristobalite_beta_minus_quartz_beta#1 | 1542.71 | 0.00000755 | 2213.30 (670.59; 0/0) | 2211.22 (668.52; 0/0) | 2210.69 (667.98; 0/0) | 2213.15 (670.44; 0/0) | 2212.37 (669.67; 0/0) | 669.44 ± 1.04 | 2833.00 |
| polynomial | sio2:cristobalite_beta_minus_tridymite_p63mmc#1 | — | — | 2187.77 (—; 1/0) | 2184.04 (—; 1/0) | 2183.76 (—; 1/0) | 2187.70 (—; 1/0) | 2186.66 (—; 1/0) | — ± — | — |
| polynomial | sio2:quartz_beta_minus_tridymite_p63mmc#1 | 1582.32 | 0.00000745 | 2209.22 (626.90; 0/0) | 2206.91 (624.59; 0/0) | 2206.42 (624.10; 0/0) | 2209.09 (626.76; 0/0) | 2208.28 (625.96; 0/0) | 625.66 ± 1.13 | 3325.67 |
| tlog_polynomial | sio2:cristobalite_beta_minus_quartz_beta#1 | 1542.71 | 0.00000755 | 2102.16 (559.45; 0/0) | 2101.62 (558.91; 0/0) | 2102.09 (559.38; 0/0) | 2101.92 (559.21; 0/0) | 2101.55 (558.85; 0/0) | 559.16 ± 0.24 | 1319.24 |
| tlog_polynomial | sio2:cristobalite_beta_minus_tridymite_p63mmc#1 | — | — | 2259.19 (—; 1/0) | 2259.13 (—; 1/0) | 2257.59 (—; 1/0) | 2258.78 (—; 1/0) | 2258.93 (—; 1/0) | — ± — | — |
| tlog_polynomial | sio2:quartz_beta_minus_tridymite_p63mmc#1 | 1582.32 | 0.00000745 | 2129.78 (547.46; 0/0) | 2129.35 (547.03; 0/0) | 2129.48 (547.16; 0/0) | 2129.53 (547.21; 0/0) | 2129.28 (546.96; 0/0) | 547.16 ± 0.17 | 1638.57 |

## Aggregate curve metrics

| predictor | G MAE (eV/atom) | ΔG MAE (eV/atom), mean across pairs/seeds | sign accuracy, mean |
|---|---:|---:|---:|
| polynomial | 4.13651 ± 0.012 | 0.012826 ± 0.0069 | 0.674752 ± 0.096 |
| tlog_polynomial | 1.36809 ± 0.0073 | 0.00651455 ± 0.0032 | 0.727673 ± 0.089 |
