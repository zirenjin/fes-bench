# Figure QA report

The Python/Matplotlib source and rendered figures were visually inspected at full output resolution.

| Figure | Static preflight | PDF text audit | Panel alignment | Supplementary text-line overlap | Official collision audit |
|---|---:|---:|---:|---:|---|
| S1_RANK_SVD | 21/21 pass | min 5.25 pt; 0 below 5 pt | PASS; 1.5 pt tolerance | 0 pairs | Blocked: PyMuPDF unavailable |
| S2_INFORMATION_LADDER | recorded prior: 21/21 pass | min 6.5 pt; 0 below 5 pt | PASS | 0 pairs | Blocked: PyMuPDF unavailable |
| S2_PAIR_CURVES_C_simple | recorded prior: 21/21 pass | min 6.2 pt; 0 below 5 pt | PASS | 0 pairs | Blocked: PyMuPDF unavailable |
| S2_PAIR_CURVES_C_DPA | recorded prior: 21/21 pass | min 6.2 pt; 0 below 5 pt | PASS | 0 pairs | Blocked: PyMuPDF unavailable |
| S2_ENVELOPE_AND_BOUNDARIES | recorded prior: 21/21 pass | min 6.2 pt; 0 below 5 pt | PASS | 0 pairs | Blocked: PyMuPDF unavailable |

The supplementary overlap check uses MuPDF structured-text line boxes. It detects intersecting text lines only; it does not inspect text crossed by vector paths, text clipping, or text over filled shapes. It is not a substitute for the required rendered collision audit. The official collision script exited 2 for all five PDFs because `fitz` is not installed. A network package installation had failed under DNS/network restrictions. Do not claim full collision QA passed.

The new S1 source-rank figure was adjusted after its first supplementary check found one overlap between exact-zero labels and the threshold note; the final export has no text-line intersections.

Browser screenshots of the updated public page also remain unavailable: the installed Node is 18, Playwright requires Node 20+, and the bundled Windows Node executable could not start from this WSL session. Static page asset links, internal paths, release hashes, and frozen source-snapshot hash were checked separately.
