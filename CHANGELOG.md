# Changelog

## 2026-09-30

- Corrected `data/processed/splits_v2/system_loso.json` to the metals-only
  folds `{hf, ti, zr}`. Each fold trains on the other two metals; SiO₂ is
  excluded. The Hf fold now records the `overlap_T` subset through 1527 K
  (the maximum Ti/Zr training temperature).
- Previous SHA256: `e8c891c1fa52dd776aaf9541e74c866ff8bc49f1c62dd3242894eb6ef16311a4`.
- Corrected SHA256: `30b680759b3ded88ccd1459413bf0914cf5a2707f43c50ad1c623a8322b08a43`.
