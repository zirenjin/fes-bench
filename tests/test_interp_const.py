import numpy as np

from fes_bench.baselines.interp_const import InterpConst


def test_interpolation_and_constant_extrapolation():
    model = InterpConst()
    model.fit([
        {"system": "x", "phase": "a", "T_K": 300.0, "G_eV_per_atom": -1.0},
        {"system": "x", "phase": "a", "T_K": 500.0, "G_eV_per_atom": -0.8},
    ])
    result = model.predict("x", "a", np.array([200.0, 400.0, 700.0]))
    np.testing.assert_allclose(result, [-1.0, -0.9, -0.8])
