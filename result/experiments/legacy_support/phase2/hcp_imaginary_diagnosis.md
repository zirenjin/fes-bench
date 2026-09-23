# Hf hcp imaginary-mode diagnosis

This finite-displacement diagnostic used the 2×2×2 supercell, a 12×12×12
Gamma-centered mesh, and the unmodified force constants at volume scale 1.0.
It is not itself a production QH result.

- Γ neighborhood: 5% of the conservative BZ radius (`0.0311751 Å⁻¹`), one
  sampled Γ point.
- Negative threshold: `-0.05 THz`.
- ASR applied: `False`; force-constant symmetry applied: `False`.

| Head | min freq (THz) | negative modes | negative fraction | negative q-points | Γ-neighborhood fraction | Γ acoustic contribution |
|---|---:|---:|---:|---:|---:|---:|
| Domains_SSE_PBE | -2.495846 | 162/798 | 0.203008 | 101/133 | 0/162 | 0/162 |
| Domains_Alloy | 0.000777913 | 0/798 | 0 | 0/133 | 0/0 | 0/0 |

## Decision

The SSE-PBE negatives are scattered across reciprocal space: 101 of 133
sampled q-points are negative, while the 5%-radius Γ neighborhood contains
0 of 162 negative modes and 0 negative acoustic entries.  This is decision
**(b)**, not a Γ-local ASR artifact.  I therefore selected the metallic-domain
`Domains_Alloy` head, did not apply ASR, and did not exclude a Γ neighborhood.
The selected head has no mode below the threshold.  A single-point sanity
check and a fresh hcp QH rerun are recorded separately before the hcp gate is
used.

The bcc gate remains false.
