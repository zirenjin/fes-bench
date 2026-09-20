from pathlib import Path

from fes_bench.models.ablation_config import load_ablation


def test_all_declared_ablations_keep_the_fixed_fes_basis_and_calibration_scope():
    root = Path(__file__).resolve().parents[1] / "configs" / "models"
    qh_only = load_ablation("qh_only", root)
    residual = load_ablation("qh_residual", root)
    repr_only = load_ablation("repr_only", root)

    assert qh_only["uses_fes_head"] is False
    assert residual["fitting_net"]["physics_baseline_column"] == 2
    assert repr_only["fitting_net"]["physics_baseline_column"] is None
    for plan in (qh_only, residual, repr_only):
        assert plan["calibration"]["fit_rows"] == "frozen_train_only"
