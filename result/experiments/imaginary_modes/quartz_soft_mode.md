# Hf hcp imaginary-mode diagnosis

This is a finite-displacement phonon diagnostic at the configured volume and mesh; it is not a production QH result.

- Γ neighborhood: 0.05 of the conservative BZ radius (0.0281043 Å⁻¹).
- Negative threshold: -0.05 THz.
- ASR applied: `False`; fc-symmetry applied: `False`.

| Head | min freq (THz) | negative modes | negative fraction | negative q-points | Γ-neighborhood fraction | Γ acoustic contribution |
|---|---:|---:|---:|---:|---:|---:|
| Domains_SSE_PBE | -1.1534 | 20 | 0.00274348 | 19 | 0.05 | 0.05 |

## Decision

The measured SSE-PBE negatives are scattered across the Brillouin-zone mesh rather than Γ-local: the Γ-neighborhood contains 1 of 20 negative modes and 1 negative acoustic entries. Therefore this is decision (b): do not apply an ASR or Γ-neighborhood exclusion; switch the representative hcp calculation to the metallic-domain `Domains_Alloy` head. The selected head is supported by the same finite-displacement test showing zero modes below the threshold. A single-point sanity check and a fresh hcp QH rerun are required before unlocking the hcp gate.

## Deviations from design

- The Γ radius is defined as 5% of a conservative Brillouin-zone radius (half the shortest reciprocal-vector norm), recorded above.
- The first pass applies neither ASR nor force-constant symmetry; both flags are recorded explicitly so a permitted follow-up can be compared against the same raw diagnosis.
- Acoustic contribution is counted as negative modes among the three lowest-frequency branches at q-points in the Γ neighborhood.
