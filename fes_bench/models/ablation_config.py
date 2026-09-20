"""Load and validate the fixed Phase-4 predictor-ablation plans."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from fes_bench.config import ConfigError, load_mapping


_VARIANTS = ("qh_only", "qh_residual", "repr_only")


def _merge(base: dict[str, Any], overlay: dict[str, Any]) -> dict[str, Any]:
    merged = dict(base)
    for key, value in overlay.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def load_ablation(name: str, config_root: Path) -> dict[str, Any]:
    """Load one declared predictor plan and enforce its non-negotiable fields.

    These are experiment-plan documents, not training launchers.  Keeping the
    validation here makes changes to the only permitted FES ablation fields
    explicit before any future runner can consume them.
    """

    if name not in _VARIANTS:
        raise ConfigError(f"unknown ablation {name!r}; choose from {_VARIANTS}")
    base = load_mapping(config_root / "model_base.yaml")
    variant = load_mapping(config_root / f"{name}.yaml")
    if variant.pop("base", None) != "model_base.yaml":
        raise ConfigError(f"{name} must explicitly inherit model_base.yaml")
    plan = _merge(base, variant)
    fitting = plan.get("fitting_net")
    if not isinstance(fitting, dict):
        raise ConfigError(f"{name} must declare fitting_net")
    if fitting.get("baseline_mode") != "additive":
        raise ConfigError(f"{name} must use additive baseline mode")
    if name == "qh_only":
        if plan.get("predictor_kind") != "qh_only" or plan.get("uses_fes_head") is not False:
            raise ConfigError("qh_only must be declared head-free")
    elif name == "qh_residual":
        if fitting.get("temperature_basis") != "continuous_tlog_polynomial":
            raise ConfigError("qh_residual must use continuous_tlog_polynomial")
        if fitting.get("physics_baseline_column") != 2 or fitting.get("numb_state_fparam") != 3:
            raise ConfigError("qh_residual must expose [T, P, F_QH] with F_QH at column 2")
    else:
        if fitting.get("temperature_basis") != "mlp":
            raise ConfigError("repr_only must use the unconstrained mlp basis")
        if fitting.get("physics_baseline_column") is not None or fitting.get("numb_state_fparam") != 2:
            raise ConfigError("repr_only must use only [T, P] and no F_QH channel")
    return plan
