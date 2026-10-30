# SiO₂ information-intervention results

## S2: what extra target information buys

On a common 1,647-point fair scoring support, the zero-shot held-out-phase pair-gap MAE is 0.4398 meV/atom for C-simple and 1.1818 meV/atom for C-DPA. Their stable-envelope accuracy is 97.57%; both predict the one reference quartz/cristobalite boundary at 1582.5 K, 39.8 K above the reference.

The **oracle-gap** intervention uses the target gap at 1,675 K to remove the local anchor error. Its pair-gap MAE is 0.4137 meV/atom for C-simple and 1.1815 meV/atom for C-DPA; the quartz/cristobalite boundary errors are 14.82 K and 35.93 K. It predicts a false cristobalite–tridymite crossing in both methods; this crossing is already present in zero-shot. This is a one-label oracle diagnostic, not the predeclared few-shot arm.

The **one-label few-shot** arm instead consumes the target label at 1,263 K. On its scored support, pair-gap MAE is 0.4961 and 1.4471 meV/atom (C-simple and C-DPA), and boundary errors are 26.58 K and 70.68 K; it also predicts a false cristobalite–tridymite crossing in both methods; the crossing is already present in zero-shot. Thus the 0.4961/1.4471 MAEs belong to the 1,263 K few-shot arm, while 14.82/35.93 K are the boundary errors for the 1,675 K oracle-gap arm.

Adding an oracle local slope strongly improves interpolation around the phase boundary: pair-gap MAE falls to 0.0790 and 0.1159 meV/atom; boundary errors fall to 0.87 and 1.26 K. This slope is estimated from 11 target temperatures (1,670–1,680 K), so it is not a one-label intervention.

Two target points at 1,263 and 2,087 K give 0.0582 and 0.0862 meV/atom, with boundary errors 7.4 and 10.8 K; neither predicts a false cristobalite–tridymite crossing. One-point correction is less reliable away from its calibration point and can move the phase boundary in the wrong direction. Pair-derived curves retain antisymmetry/cycle closure. Full-support zero-shot C-simple/C-DPA MAE and envelope scores exactly reproduce the archived V2 outputs within floating-point precision; all scalar gap comparisons differ by at most 3.47×10⁻¹⁷ eV/atom.

The improvement is diagnostic of the information needed for this one held-out phase and this reference curve. It does not establish performance across phases or materials. The oracle uses reference target information and is not a deployable prediction.

## S1: rank of the saved source readouts

With one independent supervised source-phase contrast, both archived readouts and their predicted pair-coefficient contrasts have effective rank 1 at a relative singular-value threshold of 10⁻⁴. C-DPA has default numerical matrix rank 4 for W because its remaining singular values are small but nonzero; the audit preserves those values and does not call its matrix exactly rank 1. The post-hoc reference coefficient contrasts have effective rank 2, but were fit from held-out reference curves after readouts were frozen. They are diagnostic only and cannot show that adding C2221 would cure the source-rank limit. C2221 remains inadmissible as a structure-linked source phase pending rowwise identity and source-contract closure. See [rank SVD figure](figures/S1_RANK_SVD.pdf), [source values](data/S1_RANK_SVD_FIGURE_DATA.json), and [C2221 admission audit](data/S1_RANK_SVD_AND_C2221_ADMISSION.json).

## S1 and remaining experiments

C2221 remains ineligible for a structure-linked source contrast because the table-to-geometry mapping and pressure/size-error contract are unresolved. S1 existing C readouts are available as historical evidence; new C2221 source-addition/rank experiments are blocked. S3 reuses the same C LOPO outputs; corrected T2 and true T3 LOPO were not trained, and the SiO₂ source-only control is absent. S4 and S5 are blocked by the reference decomposition and phase-specific vibrational applicability gates. Detailed status and provenance are in `S3_S4_S5_AUDIT.md` and `data/EXPERIMENT_STATUS_DETAILED.csv`.

## Figure quality assurance

The figure source is Python/Matplotlib. The four S2 PDFs passed static preflight (21/21), PDF text audit (minimum 6.2 pt), and panel alignment. The S1 rank figure also passed 21/21 preflight checks, text-size audit (minimum 5.25 pt), and panel alignment; a supplementary MuPDF text-line check found no line-box intersections. The official full collision audit remains unavailable because PyMuPDF is missing; see [Figure QA report](FIGURE_QA_REPORT.md) and `figures/official_collision_audit_status.json`. Browser screenshot QA for the updated page is also blocked by the available runtime; the manifest records this limitation.
