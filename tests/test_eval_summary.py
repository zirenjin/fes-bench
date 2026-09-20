from fes_bench.eval.run import _summary_markdown


def test_summary_includes_pairwise_table_and_fold_context():
    text = _summary_markdown(
        {
            "n_test_frames": 2,
            "G_MAE_eV_per_atom": 0.0,
            "folds": {
                "toy:alpha": {
                    "pairs": {
                        "toy:alpha_minus_beta": {
                            "pair_support": "test_and_train_partner",
                            "delta_G_MAE_eV_per_atom": 0.0,
                            "delta_G_RMSE_eV_per_atom": 0.0,
                            "sign_accuracy": 1.0,
                            "reference_Tc_K": [5.0],
                            "predicted_Tc_K": [5.0],
                            "Tc_error_K": [0.0],
                            "Tc_err_from_dG_K": [0.0],
                            "false_crossings": 0,
                            "missed_crossings": 0,
                        }
                    }
                }
            },
        },
        "toy",
    )
    assert "test_and_train_partner" in text
    assert "Tc from ΔG" in text
    assert "toy:alpha_minus_beta" in text
