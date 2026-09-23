# Phase 4 delivery — physics-baseline interface (scoped validation complete)

## Implemented scoped changes

The `feat/fes-head` checkout at
`<remote-repo>/deepmd-kit` now has a scoped change in:

- `deepmd/pt/model/task/free_energy.py`
- `deepmd/pt/model/model/free_energy_model.py`
- `deepmd/utils/argcheck.py`
- `source/tests/pt/test_fes.py`

`FreeEnergyFittingNet` accepts `physics_baseline_column`.  It extracts that
non-temperature supplied fparam column as a per-atom physical baseline,
removes it before the correction network and its input statistics, adds the
new `fes_physics_baseline` output, and serves
`fes_baseline + fes_physics_baseline + fes_correction` in additive mode.
The option is serialized and rejected for the incompatible RSTA path.

## Actual validation

```bash
cd <remote-repo>/deepmd-kit
PYTHONPATH=$PWD <remote-env>/bin/python <temporary-script>/phase4_smoke.py
```

Result: `phase4 physics-baseline smoke: PASS`.

The focused smoke constructs a three-column `[T, P, F_QH]` FES model and
verifies that changing F_QH changes only `fes_physics_baseline`, leaves the
correction invariant, and gives exact three-term total energy decomposition.
`py_compile` and `git diff --check` also passed for all four scoped files.

The frozen-export check was additionally run as an isolated temporary model:

```bash
cd <remote-repo>/deepmd-kit
PYTHONPATH=$PWD <remote-env>/bin/python <temporary-script>/phase4_deepeval_smoke.py
```

Result: `phase4 physics-baseline DeepEval smoke: PASS`.  It scripts the same
three-column model, opens it through `DeepEval`, and verifies the exported
global free energy against eager output.  The FES model is dispatched through
DeepEval's `DeepTensor` interface, so the public API exposes its declared
free-energy tensor rather than the auxiliary decomposition outputs.  New
repository tests cover both serialization preservation of
`physics_baseline_column` and the TorchScript auxiliary output.

The five new regression functions plus two pre-existing statistics regressions
were also invoked directly (using a
temporary import-only pytest shim because the environment does not ship
pytest); all passed.  The command and complete output are preserved in
`results/phase4_direct_tests_20260919.log`.

The final added regression confirms that the benchmark's selected
`continuous_tlog_polynomial` path computes finite correction statistics after
removing F_QH, while preserving the three-column external `[T, P, F_QH]`
input width.

`argcheck` now exposes `physics_baseline_column` as a nullable, zero-based
external-fparam column with explicit `[T, P, F_QH]` semantics and the same
non-temperature-column guard documented by the model.  Its direct schema
registration check printed `physics_baseline_column None`.

The scoped upstream FES test file was subsequently run on the V100 checkout
`<remote-repo>/deepmd-kit` (`feat/fes-head`):

```bash
cd <remote-repo>/deepmd-kit
PYTHONPATH=$PWD pytest -q source/tests/pt/test_fes.py
```

It returned `33 passed, 72 warnings in 11.77s`; the warnings are the existing
unknown `timeout` option and TorchScript deprecations. The complete excerpt is
preserved in `results/phase4/v100_upstream_fes_20260921.log`.

The three paper ablations are declared under `configs/models/` with one shared
`model_base.yaml` and strict loading validation in
`fes_bench/models/ablation_config.py`:

- `qh_only`: head-free `E_DPA + F_QH + c_system` plan;
- `qh_residual`: `[T, P, F_QH]`, F_QH at column 2, additive continuous-tlog
  residual head;
- `repr_only`: `[T, P]`, no F_QH channel, same continuous-tlog/additive head.

All plans retain the same five seeds and train-only one-constant-per-system
calibration policy. V100 r16 returned `15 passed in 6.20s` and directly loaded
all three plans; the command excerpt is `results/phase4/v100_r16.log`.

## Deviations from design

- `mat-agent` on the remote host still has no pytest, so its focused smoke and direct
  regressions remain the remote evidence. The V100 has now run the complete
  scoped `source/tests/pt/test_fes.py` file, but not the entire deepmd-kit
  repository test suite; no broader claim is made.
