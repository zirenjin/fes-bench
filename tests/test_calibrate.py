from __future__ import annotations

from fes_bench.models.calibrate import apply_system_constants, fit_system_constants


def test_system_calibration_cannot_absorb_phase_offset() -> None:
    train = [
        {"system": "A", "phase": "alpha", "T_index": 0},
        {"system": "A", "phase": "beta", "T_index": 0},
    ]
    reference = {("A", "alpha", 0): -1.0, ("A", "beta", 0): -0.7}
    baseline = {("A", "alpha", 0): -1.2, ("A", "beta", 0): -1.1}

    constants = fit_system_constants(train, reference, baseline)
    assert abs(constants["A"] - 0.3) < 1.0e-12
    predicted = apply_system_constants(train, baseline, constants)
    assert abs(predicted[("A", "alpha", 0)] + 0.9) < 1.0e-12
    assert abs(predicted[("A", "beta", 0)] + 0.8) < 1.0e-12
    # A shared system gauge leaves the phase difference changed by the
    # underlying baseline error; it cannot make both labels exact.
    assert abs(predicted[("A", "alpha", 0)] - predicted[("A", "beta", 0)] + 0.1) < 1.0e-12
