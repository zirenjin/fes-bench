import numpy as np

from fes_bench.baselines.bartel2018 import predict


def test_same_composition_difference_is_volume_only():
    temperature = np.array([1000.0, 1500.0])
    left = predict(temperature, -5.0, 10.0, 20.0)
    right = predict(temperature, -5.0, 12.0, 20.0)
    assert not np.allclose(left, right)
