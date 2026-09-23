"""Audit persisted benchmark artifacts without starting training.

The audit deliberately reports ``pass_with_holds`` when all infrastructure
checks pass while source or physical limitations remain.  It never treats a
missing representative, invalid QH curve, or unavailable checkpoint as a
passing substitute.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

from fes_bench.config import ConfigError, load_mapping
from fes_bench.models.ablation_config import load_ablation
from fes_bench.splits.verify import verify


FQH_FIELDS = (
    "T_K",
    "P_GPa",
    "V_min_A3_per_atom",
    "F_vib_eV_per_atom",
    "F_QH_eV_per_atom",
)


def _resolve(config_file: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else (config_file.parent / path).resolve()


def _require(path: Path, checks: dict[str, str], name: str) -> None:
    if not path.is_file():
        raise ConfigError(f"missing {name}: {path}")
    checks[name] = str(path)


def _pair_views(metrics: dict[str, Any]):
    folds = metrics.get("folds")
    return folds.values() if isinstance(folds, dict) else (metrics,)


def _check_metric(path: Path, predictor: str, checks: dict[str, str]) -> None:
    _require(path, checks, f"{predictor}:{path.parent.name}:metrics")
    try:
        metrics = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ConfigError(f"invalid metrics JSON: {path}") from exc
    if metrics.get("predictor") != predictor:
        raise ConfigError(f"{path} predictor does not match {predictor!r}")
    if predictor == "reference" and float(metrics.get("G_MAE_eV_per_atom")) != 0.0:
        raise ConfigError(f"{path} reference G MAE is not zero")
    for view in _pair_views(metrics):
        if not isinstance(view, dict):
            continue
        pairs = view.get("pairs", {})
        for pair_name, pair in pairs.items():
            if not isinstance(pair, dict):
                raise ConfigError(f"{path} pair {pair_name} is not an object")
            if predictor == "reference":
                if any(float(value) != 0.0 for value in pair.get("Tc_error_K", [])):
                    raise ConfigError(f"{path} reference Tc error is not zero")
            if pair.get("false_crossings", 0) or pair.get("missed_crossings", 0):
                raise ConfigError(f"{path} pair {pair_name} has false or missed crossings")
            if predictor.startswith("reference_noise:") and "delta_G_MAE_at_crossing_eV_per_atom" not in pair:
                raise ConfigError(f"{path} lacks crossing-local MAE")


def audit(config_path: str | Path) -> dict[str, Any]:
    config_file = Path(config_path).expanduser().resolve()
    config = load_mapping(config_file)
    checks: dict[str, str] = {}
    data_root = _resolve(config_file, str(config.get("data_root", "../data")))
    split_root = _resolve(config_file, str(config.get("split_root", "../splits")))
    results_root = _resolve(config_file, str(config.get("results_root", "../results")))
    model_root = _resolve(config_file, str(config.get("model_root", "../configs/models")))

    try:
        split_report = verify(split_root)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        raise ConfigError(f"frozen split audit failed: {exc}") from exc
    checks["frozen_splits"] = json.dumps(split_report, sort_keys=True)

    for split_name in config.get("splits", ["temp_extrap", "phase_lopo", "system_loso"]):
        _check_metric(results_root / "reference" / split_name / "metrics.json", "reference", checks)
        _require(results_root / "reference" / split_name / "summary.md", checks, f"reference:{split_name}:summary")
        _require(results_root / "reference" / split_name / "delta_g_curves.png", checks, f"reference:{split_name}:curve_png")
        _require(results_root / "reference" / split_name / "delta_g_curves.pdf", checks, f"reference:{split_name}:curve_pdf")
        _check_metric(results_root / "reference_noise_0.005" / split_name / "metrics.json", "reference_noise:0.005", checks)
        _require(results_root / "reference_noise_0.005" / split_name / "summary.md", checks, f"noise:{split_name}:summary")
        _require(results_root / "reference_noise_0.005" / "plots" / f"{split_name}.png", checks, f"noise:{split_name}:curve_png")
        _require(results_root / "reference_noise_0.005" / "plots" / f"{split_name}.pdf", checks, f"noise:{split_name}:curve_pdf")

    for phase in ("hcp", "bcc"):
        phase_root = data_root / "hf" / phase
        fqh = phase_root / "fqh.csv"
        _require(fqh, checks, f"qh:{phase}:fqh")
        with fqh.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            if tuple(reader.fieldnames or ()) != FQH_FIELDS or not next(reader, None):
                raise ConfigError(f"{fqh} does not have the canonical non-empty FQH schema")
        report = phase_root / "phonon_report.json"
        _require(report, checks, f"qh:{phase}:phonon_report")
        payload = json.loads(report.read_text(encoding="utf-8"))
        expected_reliable = {"hcp": True, "bcc": False}[phase]
        if payload.get("qh_reliable") is not expected_reliable or payload.get("diagnostic_only") is not True:
            raise ConfigError(f"{report} must record qh_reliable={expected_reliable} and diagnostic_only=true")

    for variant in ("qh_only", "qh_residual", "repr_only"):
        load_ablation(variant, model_root)
        checks[f"ablation:{variant}"] = str(model_root / f"{variant}.yaml")

    for label, setting, marker in (
        ("phase3:v100_noise_log", "phase3_noise_log", "14 passed"),
        ("phase4:v100_fes_log", "phase4_fes_log", "33 passed"),
        ("phase5:v100_recheck_log", "phase5_recheck_log", "14 passed"),
    ):
        path = _resolve(config_file, str(config.get(setting)))
        _require(path, checks, label)
        if marker not in path.read_text(encoding="utf-8"):
            raise ConfigError(f"{path} does not contain required validation marker {marker!r}")

    e3 = results_root / "phase5" / "e3_recovery.json"
    _require(e3, checks, "phase5:e3_recovery")
    e3_payload = json.loads(e3.read_text(encoding="utf-8"))
    if not e3_payload.get("acceptance_MAE_within_factor_two") or float(e3_payload.get("MAE_ratio_qh_residual_over_repr_only", 99.0)) > 2.0:
        raise ConfigError("E3 recovery acceptance failed")

    for hold in config.get("required_holds", []):
        _require(_resolve(config_file, str(hold)), checks, f"hold:{hold}")

    output = _resolve(config_file, str(config.get("output", "../result/experiments/legacy_support/completion_audit_20260921.json")))
    report = {
        "status": "pass_with_holds",
        "checks": checks,
        "holds": [str(item) for item in config.get("required_holds", [])],
        "note": "Infrastructure checks pass; listed source/physical holds remain intentionally unresolved.",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(output)
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    args = parser.parse_args(argv)
    try:
        audit(args.config)
    except (ConfigError, OSError, ValueError, KeyError) as exc:
        print(f"status audit failed: {exc}")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
