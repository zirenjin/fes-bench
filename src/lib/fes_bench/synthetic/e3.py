"""Deterministic four-phase synthetic free-energy fixture for E3."""

from __future__ import annotations

import numpy as np


def generate(seed: int = 20260919, noise_eV: float = 0.001) -> dict[str, np.ndarray]:
    """Generate four SiO2-scale curves ``a + bT + cT ln(T)`` with noise."""

    temperature = np.arange(1000.0, 2001.0)
    coefficients = np.array([
        [-8.010, -1.00e-4, 1.10e-5],
        [-8.002, -1.06e-4, 1.14e-5],
        [-8.006, -0.98e-4, 1.08e-5],
        [-7.998, -1.10e-4, 1.16e-5],
    ])
    design = np.column_stack((np.ones_like(temperature), temperature, temperature * np.log(temperature)))
    clean = coefficients @ design.T
    noisy = clean + np.random.default_rng(seed).normal(0.0, noise_eV, size=clean.shape)
    return {"T_K": temperature, "clean_G_eV_per_atom": clean, "noisy_G_eV_per_atom": noisy, "coefficients": coefficients}
