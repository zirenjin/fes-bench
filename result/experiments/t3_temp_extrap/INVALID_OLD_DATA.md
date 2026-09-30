# Invalidated pre-v2 SiO₂ T3 runs

The pre-v2 runs in the isolated Thu-GenSi directory
`/share/jzr/codex_sse_policy_20260930/fes-bench/result/experiments/`
were inspected before restarting training.  Their DeepMD inputs contain
243, 192, and 216 atoms for quartz, cristobalite, and tridymite respectively;
the canonical representatives contain 9, 24, and 12 atoms.  Their provenance
also records the old `Domains_Alloy` branch, whereas the adopted SiO₂ policy is
`Domains_SSE_PBE`.  They therefore used both the wrong structures and wrong
head, and are invalid for the benchmark.  They were moved to the remote
`result/experiments/legacy_support/t3_invalid_20260930/` archive and must not
be copied into `result/tables` or the main predictor comparison.

Canonical structure file SHA-256 values are recorded by
`src/experiments/data_prep/build_fes_training_data.py` in each generated
fold's `provenance.json`.
