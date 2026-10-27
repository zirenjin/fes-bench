# Carbon multipressure harmonic-QHA layer (v2)

## Scope and decision

This scoped extension releases paired graphite/diamond harmonic-QHA labels at 0, 10, and 15 GPa for 300, 600, 900, and 1200 K. The previously admitted 5 GPa v1 layer is preserved unchanged. The three added graphite volumes are V/V0 = 0.92, 1.10, and 1.12; diamond uses the already accepted seven-volume set. No unique atomic structure is claimed for each G(T,P) value: G and ΔG are thermodynamic ensemble quantities, while the phase-resolved structures and volume sets are archived as calculation provenance.

## Method and acceptance

VASP 5.4.4 PAW-PBE+D3(BJ), ENCUT=600 eV, static SIGMA=0.05 eV with Te=0 K; harmonic quantum nuclei; 2×2×2 phonon supercells; 16³ Gamma-centered q mesh. Added graphite volumes have paired 32-atom force calculations and 12³/16³/20³ mesh checks. At the three new volumes, 16³ minimum non-acoustic frequencies are 0.920 V/V0: 0.226544 THz, 1.100 V/V0: 0.096981 THz, 1.120 V/V0: 0.087852 THz. All q meshes have no non-acoustic imaginary modes.

Independent reconstruction used a separate standard-library q-weighted oscillator sum and polynomial stationary-point search. Result: 15 new Fvib rows and 96 phase fit cases passed; all requested phase minima are interior at polynomial degrees 2, 3, and 4.

## Paired ΔG candidates and conservative numerical bounds

ΔG is G(graphite) − G(diamond). Bounds below sum degree-2/4 fit sensitivity and q-mesh sensitivity by phase without assuming cancellation. These are numerical sensitivities only; they do not include functional, supercell, anharmonic, or electronic-temperature systematics.

| P (GPa) | T (K) | ΔG (meV/atom) | Numerical bound (meV/atom) |
|---:|---:|---:|---:|
| 0 | 300 | -108.968 | ±10.069 |
| 0 | 600 | -118.156 | ±10.251 |
| 0 | 900 | -127.999 | ±10.457 |
| 0 | 1200 | -137.526 | ±10.647 |
| 10 | 300 | 39.175 | ±11.824 |
| 10 | 600 | 29.416 | ±11.790 |
| 10 | 900 | 18.561 | ±11.759 |
| 10 | 1200 | 7.809 | ±11.720 |
| 15 | 300 | 114.286 | ±15.124 |
| 15 | 600 | 104.272 | ±14.987 |
| 15 | 900 | 92.958 | ±14.772 |
| 15 | 1200 | 81.648 | ±14.596 |

The cubic fits place the 1200 K, 10 GPa ΔG within its estimated numerical sensitivity; the sign is not resolved there. These values are scoped harmonic-QHA candidates, not complete physical free-energy truth.

## Reproducibility files

- `EXTENDED_FVT_CANDIDATE.csv`: 70 prior P5 points plus 15 new graphite F(V,T) points.
- `EXTENDED_MULTIPRESSURE_G_DELTA_CANDIDATE.csv`: degree 2/3/4 phase fits and paired ΔG candidates.
- `EXTENDED_VOLUME_HARMONIC_GATES.csv` and `EXTENDED_PHONON_QMESH_CONVERGENCE.csv`: stability and q-mesh evidence.
- `INDEPENDENT_RECONSTRUCTION_CHECK.json`: independent reconstruction result.
- `STATIC_BOUNDARY_AUDIT.csv`, `COMPUTE_LEDGER.csv`, and archived per-volume inputs/outputs: provenance.
