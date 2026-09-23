"""Run the E3 synthetic tlog-polynomial coefficient-recovery check."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from fes_bench.config import load_mapping
from fes_bench.baselines.phase_id_mlp import PhaseIdMLP
from fes_bench.synthetic.e3 import generate


def _design(temperature_K: np.ndarray) -> np.ndarray:
    temperature = np.asarray(temperature_K, dtype=float)
    return np.column_stack((np.ones_like(temperature), temperature, temperature * np.log(temperature)))


def _fit_tlog(temperature_K: np.ndarray, energy_eV_per_atom: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    design = _design(temperature_K)
    # Fixture arrays are [phase, temperature], while ``lstsq`` accepts one
    # target curve per output column.
    coefficients, _, _, _ = np.linalg.lstsq(design, np.asarray(energy_eV_per_atom).T, rcond=None)
    return coefficients.T, (design @ coefficients).T


def _fit_phase_mlp(
    temperature_K: np.ndarray,
    energy_eV_per_atom: np.ndarray,
    seed: int,
    steps: int,
    learning_rate: float,
) -> np.ndarray:
    phases = [f"phase_{index}" for index in range(np.asarray(energy_eV_per_atom).shape[0])]
    temperatures = np.tile(np.asarray(temperature_K, dtype=float), len(phases))
    phase_labels = [phase for phase in phases for _ in temperature_K]
    targets = np.asarray(energy_eV_per_atom, dtype=float).reshape(-1)
    model = PhaseIdMLP(phases, seed=seed, steps=steps, learning_rate=learning_rate)
    model.fit(temperatures, phase_labels, targets)
    return model.predict(temperatures, phase_labels).reshape(len(phases), -1)


def recover(
    seed: int = 20260919,
    noise_eV: float = 0.001,
    qh_residual_basis: str = "continuous_tlog_polynomial",
    repr_only_basis: str = "mlp",
    repr_only_steps: int = 3000,
    repr_only_learning_rate: float = 0.01,
) -> dict[str, object]:
    """Recover E3 curves with distinct residual and representation bases.

    The residual path is the constrained continuous ``T log(T)`` basis, while
    the representation-only path is the unconstrained phase-ID MLP comparator.
    E3 explicitly sets ``F_QH = 0``; it therefore tests basis behavior rather
    than claiming an atomistic QH or DPA training result.
    """

    if qh_residual_basis != "continuous_tlog_polynomial":
        raise ValueError("qh_residual_basis must be continuous_tlog_polynomial")
    if repr_only_basis != "mlp":
        raise ValueError("repr_only_basis must be mlp")

    fixture = generate(seed=seed, noise_eV=noise_eV)
    temperature = fixture["T_K"]
    truth = fixture["coefficients"]
    noisy = fixture["noisy_G_eV_per_atom"]
    fitted: dict[str, dict[str, object]] = {}
    coefficients, prediction = _fit_tlog(temperature, noisy)
    fitted["qh_residual_fqh_zero"] = {
        "basis": qh_residual_basis,
        "coefficients": coefficients.tolist(),
        "coefficient_relative_error": (
            np.abs(coefficients - truth) / np.maximum(np.abs(truth), 1.0e-12)
        ).tolist(),
        "MAE_eV_per_atom": float(np.mean(np.abs(prediction - noisy))),
    }
    mlp_prediction = _fit_phase_mlp(
        temperature,
        noisy,
        seed=seed,
        steps=repr_only_steps,
        learning_rate=repr_only_learning_rate,
    )
    fitted["repr_only"] = {
        "basis": repr_only_basis,
        "coefficients": None,
        "coefficient_relative_error": None,
        "MAE_eV_per_atom": float(np.mean(np.abs(mlp_prediction - noisy))),
    }
    ratio = fitted["qh_residual_fqh_zero"]["MAE_eV_per_atom"] / fitted["repr_only"]["MAE_eV_per_atom"]
    return {
        "seed": seed,
        "noise_eV_per_atom": noise_eV,
        "n_phases": int(noisy.shape[0]),
        "truth_coefficients": truth.tolist(),
        "fits": fitted,
        "basis": {
            "qh_residual_fqh_zero": qh_residual_basis,
            "repr_only": repr_only_basis,
        },
        "MAE_ratio_qh_residual_over_repr_only": float(ratio),
        "acceptance_MAE_within_factor_two": bool(0.5 <= ratio <= 2.0),
    }


def run(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    args = parser.parse_args(argv)
    config = load_mapping(args.config)
    result = recover(
        seed=int(config.get("seed", 20260919)),
        noise_eV=float(config.get("noise_eV", 0.001)),
        qh_residual_basis=str(config.get("qh_residual_basis", "continuous_tlog_polynomial")),
        repr_only_basis=str(config.get("repr_only_basis", "mlp")),
        repr_only_steps=int(config.get("repr_only_steps", 3000)),
        repr_only_learning_rate=float(config.get("repr_only_learning_rate", 0.01)),
    )
    output = Path(str(config["output"])).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
