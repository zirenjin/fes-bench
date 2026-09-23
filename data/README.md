# Public data manifest

The benchmark keeps public archives out of git. `data/download.sh` downloads
the currently configured SiO₂ and Hf archives and verifies their SHA-256
values; the complete per-file provenance is in
`data/processed/sio2/download.json` and `data/processed/hf/download.json`.
Use `bash data/download.sh --dry-run` to inspect the
entry point without downloading.

| dataset | DOI / record | local landing area | status |
|---|---|---|---|
| SiO₂ direct upsampling | [10.18419/DARUS-4999](https://doi.org/10.18419/DARUS-4999) | `data/raw/sio2/`, processed files under `data/processed/sio2/` | active |
| Ti/Zr/Hf thermodynamic tables | [10.18419/DARUS-3582](https://doi.org/10.18419/DARUS-3582) | `data/raw/{ti,zr,hf}/`, processed files under `data/processed/{ti,zr,hf}/` | active where manifest entries exist |
| CaSiO₃ DP-TI | [Zenodo 10460440](https://doi.org/10.5281/zenodo.10460440) | `data/raw/casio3/`, processed files under `data/processed/casio3/` | hold; not downloaded |

Every downloaded archive must be listed with URL, DOI, destination, download
time, action, and SHA-256 in the corresponding `download.json`. Raw archives
are ignored by git; normalized processed inputs are versioned.
