"""Per-phase training-range interpolation with constant extrapolation."""

from __future__ import annotations

import numpy as np


class InterpConst:
    def __init__(self) -> None:
        self._curves: dict[tuple[str, str], tuple[np.ndarray, np.ndarray]] = {}

    def fit(self, frames: list[dict[str, object]]) -> None:
        grouped: dict[tuple[str, str], list[tuple[float, float]]] = {}
        for frame in frames:
            key = (str(frame["system"]), str(frame["phase"]))
            grouped.setdefault(key, []).append((float(frame["T_K"]), float(frame["G_eV_per_atom"])))
        self._curves = {}
        for key, values in grouped.items():
            values.sort()
            temperature, energy = map(np.asarray, zip(*values, strict=True))
            if len(np.unique(temperature)) != len(temperature):
                raise ValueError(f"duplicate training temperatures for {key}")
            self._curves[key] = temperature, energy

    def predict(self, system: str, phase: str, temperature_K: np.ndarray) -> np.ndarray:
        try:
            temperature, energy = self._curves[(system, phase)]
        except KeyError as exc:
            raise ValueError(f"no training curve for {(system, phase)}") from exc
        return np.interp(np.asarray(temperature_K, dtype=float), temperature, energy, left=energy[0], right=energy[-1])
