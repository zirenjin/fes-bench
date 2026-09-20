import numpy as np

from fes_bench.baselines.phase_id_mlp import PhaseIdMLP


def test_phase_id_mlp_fits_small_temperature_function():
    temperature = np.tile(np.linspace(300.0, 900.0, 10), 2)
    phases = ["a"] * 10 + ["b"] * 10
    energy = np.concatenate((-1.0 + temperature[:10] * 1e-4, -0.7 + temperature[10:] * 1e-4))
    model = PhaseIdMLP(["a", "b"], steps=600, learning_rate=1e-2)
    model.fit(temperature, phases, energy)
    predicted = model.predict(temperature, phases)
    assert predicted.shape == (20,)
    assert np.mean(np.abs(predicted - energy)) < 0.01
