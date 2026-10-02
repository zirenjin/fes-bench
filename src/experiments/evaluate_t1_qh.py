"""Produce the T1 E + F_QH controlled-ablation predictor — plan experiment T1.

The script reads canonical representative energies and QH curves, constructs
``G_T1(T) = E_DPA + F_QH(T)``, and delegates every metric (including frozen
split filtering and crossing diagnostics) to the shared evaluator.  It never
trains a model and writes a provenance-bearing raw run for each split.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from fes_bench.eval.run import evaluate


SYSTEM_PHASES = {
    "sio2": ("quartz_beta", "cristobalite_beta", "tridymite_p63mmc"),
    "hf": ("hcp", "bcc"),
    "ti": ("hcp", "bcc"),
    "zr": ("hcp", "bcc"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_commit(root: Path) -> str:
    try:
        return subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def _json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"{path} must contain an object")
    return payload


def _energy(meta: dict[str, Any], path: Path) -> float:
    value = meta.get("formal_energy_eV_per_atom")
    if value is None:
        value = meta.get("representative", {}).get("energy_eV_per_atom", {}).get("after")
    if not isinstance(value, (int, float)):
        raise ValueError(f"{path} has no canonical E_DPA energy")
    return float(value)


def _fqh(path: Path) -> dict[float, float]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or "T_K" not in rows[0] or "F_QH_eV_per_atom" not in rows[0]:
        raise ValueError(f"{path} is missing the canonical F_QH_eV_per_atom column")
    values: dict[float, float] = {}
    for row in rows:
        temperature, value = float(row["T_K"]), float(row["F_QH_eV_per_atom"])
        if temperature in values:
            raise ValueError(f"{path} has duplicate temperature {temperature:g}")
        values[temperature] = value
    return values


def build_curves(root: Path) -> tuple[dict[str, Any], list[Path], dict[str, dict[str, Any]]]:
    predictions: dict[str, Any] = {}
    inputs: list[Path] = []
    phase_meta: dict[str, dict[str, Any]] = {}
    for system, phases in SYSTEM_PHASES.items():
        for phase in phases:
            phase_dir = root / "data/processed" / system / phase
            meta_path = phase_dir / "meta.json"
            fqh_path = phase_dir / "fqh.csv"
            meta = _json(meta_path)
            fqh = _fqh(fqh_path)
            e_dpa = _energy(meta, meta_path)
            report_path = phase_dir / "phonon_report.json"
            report = _json(report_path) if report_path.exists() else {}
            predictions[f"{system}:{phase}"] = {
                "T_K": [temperature for temperature in sorted(fqh)],
                "G_eV_per_atom": [e_dpa + fqh[temperature] for temperature in sorted(fqh)],
            }
            phase_meta[f"{system}:{phase}"] = {
                "qh_reliable": report.get("qh_reliable", meta.get("qh_reliable")),
                "e_dpa_eV_per_atom": e_dpa,
                "energy_head": meta.get("energy_head_policy") or meta.get("representative", {}).get("head"),
                "meta": str(meta_path.relative_to(root)),
                "fqh": str(fqh_path.relative_to(root)),
            }
            inputs.extend([meta_path, fqh_path])
            if report_path.exists():
                inputs.append(report_path)
    return predictions, inputs, phase_meta


def _with_overlap(split: dict[str, Any], curves_path: Path, data_root: Path) -> dict[str, Any]:
    """Evaluate a split and attach the Hf overlap_T subset when present."""
    result = evaluate(f"precomputed:{curves_path}", split, data_root, [0])
    if isinstance(split.get("folds"), dict) and isinstance(result.get("folds"), dict):
        for fold_name, source_fold in split["folds"].items():
            if not isinstance(source_fold, dict) or not isinstance(source_fold.get("overlap_T"), list):
                continue
            subset = {"test": source_fold["overlap_T"], "train": source_fold.get("train", [])}
            result["folds"][str(fold_name)]["overlap_T"] = evaluate(
                f"precomputed:{curves_path}", subset, data_root, [0]
            )
    return result


def _json_or_marker(value: Any) -> str:
    return json.dumps(value, separators=(",", ":"), ensure_ascii=False) if isinstance(value, list) else str(value)


def write_pair_details(root: Path, split_name: str, metrics: dict[str, Any], phase_meta: dict[str, dict[str, Any]]) -> Path:
    path = root / "result/experiments/t1_qh" / f"pair_details_{split_name}.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "split", "fold", "eval_subset", "pair", "left_qh_reliable", "right_qh_reliable",
        "delta_G_MAE_eV_per_atom", "sign_accuracy", "Tc_error_K", "Tc_err_from_dG_K",
        "false_crossings", "missed_crossings",
    ]
    rows: list[dict[str, Any]] = []
    def add(fold_name: str, subset: str, view: dict[str, Any]) -> None:
        for pair, record in (view.get("pairs", {}) if isinstance(view.get("pairs"), dict) else {}).items():
            if ":" not in pair:
                continue
            system, pair_name = pair.split(":", 1)
            if "_minus_" in pair_name:
                left, right = pair_name.split("_minus_", 1)
            else:
                left, right = pair_name, ""
            left_meta = phase_meta.get(f"{system}:{left}", {})
            right_meta = phase_meta.get(f"{system}:{right}", {})
            rows.append({
                "split": split_name,
                "fold": fold_name,
                "eval_subset": subset,
                "pair": pair,
                "left_qh_reliable": left_meta.get("qh_reliable", ""),
                "right_qh_reliable": right_meta.get("qh_reliable", ""),
                "delta_G_MAE_eV_per_atom": record.get("delta_G_MAE_eV_per_atom", ""),
                "sign_accuracy": record.get("sign_accuracy", ""),
                "Tc_error_K": _json_or_marker(record.get("Tc_error_K", "")),
                "Tc_err_from_dG_K": _json_or_marker(record.get("Tc_err_from_dG_K", "")),
                "false_crossings": record.get("false_crossings", ""),
                "missed_crossings": record.get("missed_crossings", ""),
            })
    if isinstance(metrics.get("folds"), dict):
        for fold_name, view in metrics["folds"].items():
            if isinstance(view, dict):
                add(str(fold_name), "full", view)
                if isinstance(view.get("overlap_T"), dict):
                    add(str(fold_name), "overlap_T", view["overlap_T"])
    else:
        add("all", "full", metrics)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return path


def run(root: Path, split_names: list[str]) -> int:
    root = root.resolve()
    data_root = root / "data/processed"
    curves, curve_inputs, phase_meta = build_curves(root)
    curves_path = root / "result/experiments/t1_qh/curves.json"
    curves_path.parent.mkdir(parents=True, exist_ok=True)
    curves_path.write_text(json.dumps({"predictions": curves}, indent=2) + "\n", encoding="utf-8")
    head_policy = root / "configs/models/head_policy.yaml"
    predictor_config = root / "configs/predictors/qh_only.yaml"
    checkpoint_sha = ""
    summary = root / "result/experiments/quasi_harmonic_sio2_sse_pbe/raw_runs/sio2/qh_summary.json"
    if summary.exists():
        checkpoint_sha = str(_json(summary).get("checkpoint_sha256", ""))
    for split_name in split_names:
        split_path = root / "data/processed/splits_v2" / f"{split_name}.json"
        split = _json(split_path)
        metrics = _with_overlap(split, curves_path, data_root)
        metrics["method"] = "qh_only"
        raw_path = root / "result/experiments/external_baselines" / f"raw_{split_name}" / "qh_only" / "seed_none" / "metrics.json"
        raw_path.parent.mkdir(parents=True, exist_ok=True)
        raw_payload = {
            "provenance": {
                "git_commit": git_commit(root),
                "split": split_name,
                "split_sha256": sha256(split_path),
                "predictor_config": str(predictor_config.relative_to(root)),
                "predictor_config_sha256": sha256(predictor_config),
                "checkpoint_path": "external/checkpoints/DPA-3.1-3M.pt",
                "checkpoint_sha256": checkpoint_sha,
                "head_policy": str(head_policy.relative_to(root)),
                "head_policy_sha256": sha256(head_policy),
                "evaluator_version": "fes_bench.eval.run precomputed curves",
                "timestamp_utc": datetime.now(timezone.utc).isoformat(),
                "seed_policy": "n/a:no_seed",
                "formula": "G_T1(T) = E_DPA + F_QH(T)",
                "fqh_imaginary_modes": "excluded by adopted QH producer",
                "phase_meta": phase_meta,
                "input_sha256": [{"path": str(path.relative_to(root)), "sha256": sha256(path)} for path in sorted(set(curve_inputs))],
            },
            "metrics": metrics,
        }
        raw_path.write_text(json.dumps(raw_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        detail_path = write_pair_details(root, split_name, metrics, phase_meta)
        print(f"T1 {split_name}: {raw_path}")
        print(f"pair details: {detail_path}")
    curves_path.with_suffix(".meta.json").write_text(
        json.dumps({"git_commit": git_commit(root), "inputs": [{"path": str(path.relative_to(root)), "sha256": sha256(path)} for path in sorted(set(curve_inputs))], "formula": "G_T1(T) = E_DPA + F_QH(T)", "checkpoint_sha256": checkpoint_sha, "head_policy": str(head_policy.relative_to(root)), "qh_imaginary_modes": "excluded by adopted QH producer"}, indent=2) + "\n",
        encoding="utf-8",
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--splits", default="temp_extrap,phase_lopo,system_loso")
    args = parser.parse_args(argv)
    return run(args.repo_root, [item.strip() for item in args.splits.split(",") if item.strip()])


if __name__ == "__main__":
    raise SystemExit(main())
