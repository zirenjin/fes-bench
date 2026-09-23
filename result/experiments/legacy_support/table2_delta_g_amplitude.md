# Table 2 candidate — reference ΔG amplitude

Amplitude is the population standard deviation of pair ΔG over the full common reference grid (eV/atom; multiply by 1000 for meV/atom). ΔG is defined as `G(left) − G(right)`.

| system | pair | grid | ΔG std (meV/atom) | min–max (meV/atom) |
|---|---|---|---:|---:|
| SiO₂ | cristobalite_beta − quartz_beta | 851–2499 K, 1649 points | 3.607 | −7.494 to +4.884 |
| SiO₂ | cristobalite_beta − tridymite_p63mmc | 851–2499 K, 1649 points | 0.059 | −0.317 to −0.109 |
| SiO₂ | quartz_beta − tridymite_p63mmc | 851–2499 K, 1649 points | 3.576 | −4.993 to +7.304 |
| Hf | bcc − hcp | 1076–2327 K, 1252 points | 11.987 | −9.300 to +32.000 |

The cristobalite–tridymite pair has no reference crossing and a much smaller amplitude; this is why it is retained for diagnostics but excluded from the recommended primary aggregate.

Source: `results/reference_delta_g_stats.json`; `std` uses population normalization (`ddof=0`).

## Deviations from design

No deviation. Values are a direct aggregation of the persisted reference-statistics artifact; no new interpolation or metric definition was introduced.
