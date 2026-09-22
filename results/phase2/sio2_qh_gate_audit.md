# SiO₂ QH gate audit — code and cell geometry

## F_QH imaginary-frequency handling (source inspection)

In `fes_bench/qh/run.py`, `_thermal_for_volume` calls:

```python
phonon.run_thermal_properties(
    t_min=min(temperatures), t_max=max(temperatures), t_step=1.0
)
```

No `cutoff_frequency` or `pretend_real` argument is passed. The installed
phonopy 2.29.0 source defines `cutoff_frequency=None` as a zero cutoff in
`Phonopy.run_thermal_properties`: modes with frequencies below zero are
excluded from the thermal free-energy sum; frequencies are not made real by
absolute value. The configuration key
`imaginary_frequency_cutoff_THz: -0.05` is therefore metadata only in this
run and was not forwarded to phonopy. The raw mesh minimum is still recorded
in each `*_fqh.csv`, and the separate diagnose pass measures the negative-mode
fractions and BZ locations.

## Actual `[1,1,1]` supercell edge lengths

The QH input structures have these lattice-vector lengths; with
`supercell_matrix: [1,1,1]`, the phonopy supercell is exactly the input cell.

| phase | lattice-vector lengths (Å) | shortest edge (Å) | volume (Å³) |
|---|---|---:|---:|
| quartz_beta | 14.517843, 14.517843, 15.975435 | **14.517843** | 2916.000001 |
| cristobalite_beta | 13.998051, 13.998051, 13.998051 | **13.998051** | 2742.854400 |
| tridymite_p63mmc | 14.726499, 14.726499, 16.100973 | **14.726499** | 3023.999999 |

Evidence: `data/sio2/*/structure.extxyz` in the isolated thu-GenSi workspace;
values were computed from the `Lattice` vectors with `numpy.linalg.norm` and
`abs(det(cell))`.

## Gate decision

Until this source-level cutoff clarification and the actual cell geometry are
explicitly recorded, the SiO₂ QH gate is **pending**. This audit does not
change the existing QH curves or the frozen protocol.
