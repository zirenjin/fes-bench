# Run status

## Completed

- Froze the S2 protocol before scoring.
- Scored the C-simple/C-DPA information ladder on common fair support and crosschecked the original V2 zero-shot result.
- Completed the C2221 admission and rank audit; pressure rounding impact is quantified.
- Reused and rescored the archived C readouts under the same LOPO scorer for S3.
- Built the public review / gh-pages release with the S2 page, data, plots and portable generator scripts. Installed the update into the existing `gh-pages` working tree without committing or pushing.

## Blocked or unavailable

- C2221 cannot yet serve as a structure-linked source phase: rowwise source geometry, pressure and size-error contracts remain unresolved.
- Corrected T2 and phase-LOPO T3 need retraining. Existing T2 double-counts its static term; available T3 evidence is temperature-tail only. The local DeepMD runtime and remote Expanse connection were unavailable. The SiO2 source-only LOPO readout is absent.
- S4 lacks matched same-branch static energies and geometry identity. S5 lacks a common vibrationally valid phase subset.
- Browser screenshot QA could not run because bundled Playwright requires Node 20+ and the Windows Node executable could not launch from this WSL session. PDF collision QA is unresolved because PyMuPDF is unavailable and package installation failed under DNS/network restrictions.

The gh-pages worktree now contains local modifications on branch `gh-pages`; it remains uncommitted and unpublished to the remote. The original frozen review snapshot hash is preserved. Detailed status is in `EXPERIMENT_STATUS.csv` and `data/EXPERIMENT_STATUS_DETAILED.csv`.
