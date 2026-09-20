from fes_bench.synthetic.e3_recovery import recover


def test_e3_recovery_uses_distinct_bases_and_stays_within_factor_two():
    result = recover()
    assert result["n_phases"] == 4
    assert result["acceptance_MAE_within_factor_two"]
    assert result["basis"] == {
        "qh_residual_fqh_zero": "continuous_tlog_polynomial",
        "repr_only": "mlp",
    }
    assert result["MAE_ratio_qh_residual_over_repr_only"] != 1.0
    assert result["fits"]["repr_only"]["coefficients"] is None
