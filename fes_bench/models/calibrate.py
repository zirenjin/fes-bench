"""System-level energy-gauge calibration for benchmark predictors.

The calibration is deliberately one constant per chemical system, never one
constant per phase.  Fitting it only from frozen training triples prevents a
predictor from absorbing the target relative phase free energies.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping


class CalibrationError(ValueError):
    """Raised for incomplete or inconsistent calibration inputs."""


def fit_system_constants(
    train_rows: Iterable[Mapping[str, object]],
    reference_eV_per_atom: Mapping[tuple[str, str, int], float],
    baseline_eV_per_atom: Mapping[tuple[str, str, int], float],
) -> dict[str, float]:
    """Fit ``mean(reference - baseline)`` separately for each system.

    ``train_rows`` must be the frozen split's training triples.  The lookup
    mappings are explicit to keep this routine model-agnostic and make it
    impossible for a hidden test frame to participate in the fit.
    """

    residuals: dict[str, list[float]] = defaultdict(list)
    for row in train_rows:
        try:
            key = (str(row["system"]), str(row["phase"]), int(row["T_index"]))
            reference = reference_eV_per_atom[key]
            baseline = baseline_eV_per_atom[key]
        except (KeyError, TypeError, ValueError) as exc:
            raise CalibrationError(f"missing or invalid calibration row: {row!r}") from exc
        residuals[key[0]].append(float(reference) - float(baseline))
    if not residuals:
        raise CalibrationError("cannot calibrate an empty training set")
    return {system: sum(values) / len(values) for system, values in residuals.items()}


def apply_system_constants(
    rows: Iterable[Mapping[str, object]],
    baseline_eV_per_atom: Mapping[tuple[str, str, int], float],
    constants_eV_per_atom: Mapping[str, float],
) -> dict[tuple[str, str, int], float]:
    """Apply previously fitted constants without inspecting reference labels."""

    predicted: dict[tuple[str, str, int], float] = {}
    for row in rows:
        key = (str(row["system"]), str(row["phase"]), int(row["T_index"]))
        if key[0] not in constants_eV_per_atom:
            raise CalibrationError(f"missing system calibration constant for {key[0]!r}")
        try:
            predicted[key] = float(baseline_eV_per_atom[key]) + float(constants_eV_per_atom[key[0]])
        except KeyError as exc:
            raise CalibrationError(f"missing baseline for {key!r}") from exc
    return predicted
