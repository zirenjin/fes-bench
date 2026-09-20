import numpy as np

from fes_bench.synthetic.e3 import generate


def test_e3_is_deterministic_and_four_phase():
    first, second = generate(), generate()
    assert first["clean_G_eV_per_atom"].shape[0] == 4
    assert first["T_K"][0] == 1000.0
    assert first["T_K"][-1] == 2000.0
    assert np.array_equal(first["noisy_G_eV_per_atom"], second["noisy_G_eV_per_atom"])
