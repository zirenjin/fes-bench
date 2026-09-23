# thu-GenSi environment lock

Captured 2026-09-22 UTC in the isolated workspace
`/share/jzr/codex_fes_overnight_20260922_1790052610710037915`.

| field | value |
|---|---|
| source tree commit at capture | `1260f19a8811eaf7c13e7566df97411e7b9d616b` |
| execution environment | `/root/miniconda3/envs/mat-agent` |
| Python | 3.11.15 |
| conda export | `env_lock_conda.yml` in the isolated workspace; SHA-256 `8c2eb1a3dc702990cb9348626ce00086d37f8682a6173282d890293b7e152d62` |
| pip freeze | `env_lock_pip_freeze.txt` in the isolated workspace; 183 lines; SHA-256 `79f437521671a3fca301f3d09f40120d53b18ba46132720b508dac5612133a2e` |
| remote git metadata | absent in the isolated workspace; the source commit above is the local repository commit used for synchronization |

The exports were generated with:

```text
/root/miniconda3/bin/conda env export -n mat-agent --no-builds > env_lock_conda.yml
/root/miniconda3/envs/mat-agent/bin/python -m pip freeze > env_lock_pip_freeze.txt
```

Key runtime packages: `deepmd-kit==0.1.dev4172+g4889cf879`, `phonopy==2.29.0`,
`torch==2.10.0`, `numpy==1.26.4`, `scipy==1.17.1`, `pymatgen==2026.5.4`,
and `matplotlib==3.10.9`.

## Resource status

thu-GenSi is the only environment used for Phase 6 readiness. The V100 is
listed separately as a backup resource; its balance is pending recharge and it
is not a readiness gate.
