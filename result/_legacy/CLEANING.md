# Legacy cleaning record

The six legacy delivery records retain their original wording, conclusions,
acceptance results, and deviation statements. Only machine-specific paths and
internal host identifiers were replaced before the records were added to Git.

| Original form | Clean form | Applies to |
|---|---|---|
| `/mnt/c/Users/ziren/Documents/dpa-adapt/fes-bench/data` | `<repo-root>/data` | phase0 |
| `/GenSIvePFS/users/zirenj/fes-bench` | `<remote-repo>/fes-bench` | phase1 |
| `/GenSIvePFS/users/zirenj/deepmd-kit` | `<remote-repo>/deepmd-kit` | phase4 |
| `/share/jzr/fes-bench/runs/...` | `<remote-run>/...` | phase3, phase5 |
| `/share/jzr/deepmd-kit` | `<remote-repo>/deepmd-kit` | phase4 |
| `/root/fes-bench-env/bin/python` | `<remote-env>/bin/python` | phase3, phase5 |
| `/root/miniconda3/envs/mat-agent/bin/python` | `<remote-env>/bin/python` | phase4 |
| `/tmp/phase4_*.py` | `<temporary-script>/phase4_*.py` | phase4 |
| `thu-GenSi` | `remote host` | phase1, phase2, phase3, phase4 |
| host-specific audit/log path names containing `thu` | corresponding `remote` path names | phase1, phase3 |

No numerical value, date, test result, scientific conclusion, or deviation
statement was changed.
