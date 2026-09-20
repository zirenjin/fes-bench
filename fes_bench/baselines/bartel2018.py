"""Bartel et al. (2018) temperature-dependent free-energy baseline.

For a phase with 0-K energy ``E0`` (eV/atom), volume ``V`` (A^3/atom), and
reduced atomic mass ``m`` (amu), the implemented correction is
``(-2.48e-4 ln(V) - 8.94e-5 m/V) T + 0.181 ln(T) - 0.882`` eV/atom.

The coefficients, units, and fixed-composition polymorph limitation were
checked against Eq. 4 of Bartel *et al.*, Nat. Commun. 9, 4168 (2018),
https://doi.org/10.1038/s41467-018-06682-4 (open-access article and its
linked supplementary material, checked 2026-09-19).

Within a fixed-composition polymorph family ``m`` is identical, leaving volume
as the only structural discriminator.  This limitation is intentional and is
reported alongside benchmark results.
"""

from __future__ import annotations

import numpy as np


def gibbs_correction(T_K: np.ndarray | float, volume_A3_per_atom: float, reduced_mass_amu: float) -> np.ndarray:
    """Return the Bartel free-energy correction in eV/atom."""

    temperature = np.asarray(T_K, dtype=float)
    if np.any(temperature <= 0) or volume_A3_per_atom <= 0 or reduced_mass_amu <= 0:
        raise ValueError("temperature, volume, and reduced mass must be positive")
    slope = -2.48e-4 * np.log(volume_A3_per_atom) - 8.94e-5 * reduced_mass_amu / volume_A3_per_atom
    return slope * temperature + 0.181 * np.log(temperature) - 0.882


def predict(T_K: np.ndarray | float, E0_eV_per_atom: float, volume_A3_per_atom: float, reduced_mass_amu: float) -> np.ndarray:
    """Return ``E0 + G_delta`` in eV/atom."""

    return float(E0_eV_per_atom) + gibbs_correction(T_K, volume_A3_per_atom, reduced_mass_amu)
