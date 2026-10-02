"""Training-only thermo-form ΔG predictor.

The predictor fits ``a + b T + c T log(T)`` to shared training-region
reference ΔG values and evaluates the fitted curve at the requested
temperatures.  It never reads test-region labels.
"""

from __future__ import annotations

import numpy as np


def design(temperature: np.ndarray) -> np.ndarray:
    """Return the three-column thermo-form design matrix."""

    temperature = np.asarray(temperature, dtype=float)
    if np.any(temperature <= 0.0):
        raise ValueError("thermo-form temperatures must be positive")
    return np.column_stack((np.ones_like(temperature), temperature, temperature * np.log(temperature)))


def fit(temperature: np.ndarray, delta_g: np.ndarray) -> np.ndarray:
    """Fit ``a + b T + c T log(T)`` and return its coefficients."""

    temperature = np.asarray(temperature, dtype=float)
    delta_g = np.asarray(delta_g, dtype=float)
    if temperature.shape != delta_g.shape or temperature.ndim != 1 or len(temperature) < 3:
        raise ValueError("thermo-form fitting requires three or more aligned points")
    coefficients, *_ = np.linalg.lstsq(design(temperature), delta_g, rcond=None)
    return coefficients


def predict(temperature: np.ndarray, coefficients: np.ndarray) -> np.ndarray:
    """Evaluate a fitted thermo-form curve."""

    return design(np.asarray(temperature, dtype=float)) @ np.asarray(coefficients, dtype=float)


def fit_predict(train_temperature: np.ndarray, train_delta_g: np.ndarray, temperature: np.ndarray) -> np.ndarray:
    """Fit on training points and predict at ``temperature``."""

    return predict(temperature, fit(train_temperature, train_delta_g))
