"""Normalize imported metrics into result/experiments raw outputs — migration support."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def commit(root: Path) -> str:
    try:
        return subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def sanitize(value: Any) -> Any:
    """Keep archived payloads repository-relative without changing numeric values."""
    if isinstance(value, dict):
        return {str(key): sanitize(item) for key, item in value.items()}
    if isinstance(value, list):
        return [sanitize(item) for item in value]
    if isinstance(value, str) and (value.startswith("/") or ":/" in value):
        prefix = ""
        if ":/" in value and not value.startswith("/"):
            prefix, value = value.split(":/", 1)
            prefix += ":"
        return f"{prefix}external/{Path(value).name}"
    return value


def provenance(root: Path, source: Path, split: str, predictor: str, seed: str, missing: list[dict[str, str]]) -> dict[str, Any]:
    config = root / "configs" / "predictors" / f"{predictor}.yaml"
    return {
        "git_commit": commit(root),
        "split": split,
        "split_sha256": sha256(root / "data/processed/splits" / f"{split}.json") if (root / "data/processed/splits" / f"{split}.json").exists() else None,
        "predictor_config": str(config.relative_to(root)) if config.exists() else None,
        "predictor_config_sha256": sha256(config) if config.exists() else None,
        "checkpoint_path": None,
        "checkpoint_sha256": None,
        "evaluator_version": "fes-bench historical result normalizer v1",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "input_metrics": str(source.relative_to(root)),
        "input_metrics_sha256": sha256(source),
        "missing": missing,
    }


def write_run(root: Path, split: str, predictor: str, seed: str, payload: Any, source: Path, missing: list[dict[str, str]]) -> None:
    if split == "e1_full_grid":
        out = root / "result/experiments/crossing_reevaluation/raw_runs" / predictor / f"seed_{seed}" / "metrics.json"
    elif split == "e2_calibration_window":
        out = root / "result/experiments/calibration_window/raw_runs" / predictor / f"seed_{seed}" / "metrics.json"
    elif split == "e3_synthetic_recovery":
        out = root / "result/experiments/synthetic_recovery/raw_runs" / predictor / f"seed_{seed}" / "metrics.json"
    elif split == "reference_statistics":
        out = root / "result/experiments/reference_statistics/raw_runs" / predictor / f"seed_{seed}" / "metrics.json"
    else:
        out = root / "result/experiments/external_baselines" / f"raw_{split}" / predictor / f"seed_{seed}" / "metrics.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    body = {"provenance": provenance(root, source, split, predictor, seed, missing), "metrics": sanitize(payload)}
    out.write_text(json.dumps(body, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def normalize(root: Path, manifest: dict[str, Any]) -> int:
    for spec in manifest["runs"]:
        source = root / spec["source"]
        if not source.exists():
            raise FileNotFoundError(source)
        payload = json.loads(source.read_text(encoding="utf-8"))
        mode = spec["mode"]
        if mode == "nested":
            for split in sorted(payload):
                if spec.get("split") not in (None, "*") and split != spec["split"]:
                    continue
                for predictor in sorted(payload[split]):
                    if spec.get("predictor") not in (None, "*") and predictor != spec["predictor"]:
                        continue
                    write_run(root, split, predictor, str(spec.get("seed", "none")), payload[split][predictor], source, [])
        elif mode == "flat":
            write_run(root, spec["split"], spec["predictor"], str(spec.get("seed", "none")), payload, source, [])
        elif mode == "records":
            for record in payload["checkpoint_records"]:
                predictor = str(record[spec.get("predictor_field", "basis")])
                seed = str(record.get("seed", "record"))
                missing = [{"field": "checkpoint_path", "reason": "historical external checkpoint is not in this repository"}]
                write_run(root, spec["split"], predictor, seed, record, source, missing)
        else:
            raise ValueError(f"Unknown mode: {mode}")
    build_inventory(root)
    return 0


def build_inventory(root: Path) -> None:
    systems: dict[str, Any] = {}
    phases: dict[str, Any] = {}
    for system_path in sorted((root / "data/processed").glob("*/system.json")):
        system = system_path.parent.name
        system_data = json.loads(system_path.read_text(encoding="utf-8"))
        systems[system] = sanitize(system_data)
        for phase in system_data.get("phases", []):
            phase_dir = root / "data" / "processed" / system / phase
            meta_path = phase_dir / "meta.json"
            phase_data = json.loads(meta_path.read_text(encoding="utf-8")) if meta_path.exists() else {}
            table_path = phase_dir / "reference_G.csv"
            n_points = sum(1 for _ in table_path.open(encoding="utf-8")) - 1 if table_path.exists() else None
            edge_lengths: list[float] = []
            structure_path = phase_dir / "structure.extxyz"
            atom_count = None
            if structure_path.exists():
                first_lines = structure_path.read_text(encoding="utf-8").splitlines()[:2]
                if first_lines:
                    try:
                        atom_count = int(first_lines[0].strip())
                    except ValueError:
                        atom_count = None
                match = re.search(r'Lattice="([^"]+)"', " ".join(first_lines))
                if match:
                    values = [float(value) for value in match.group(1).split()]
                    edge_lengths = [sum(values[row * 3 + col] ** 2 for col in range(3)) ** 0.5 for row in range(3)]
            phases[f"{system}:{phase}"] = {
                "system": system,
                "phase": phase,
                "meta": sanitize(phase_data),
                "reference_points": n_points,
                "structure_path": str((phase_dir / "structure.extxyz").relative_to(root)) if (phase_dir / "structure.extxyz").exists() else None,
                "lattice_edge_lengths_A": edge_lengths,
                "shortest_edge_A": min(edge_lengths) if edge_lengths else None,
                "atom_count": atom_count,
            }
    splits: dict[str, Any] = {}
    for split_path in sorted((root / "data/processed/splits").glob("*.json")):
        split = json.loads(split_path.read_text(encoding="utf-8"))
        if "folds" in split:
            folds = {name: {subset: len(value.get(subset, [])) for subset in ("train", "test")} for name, value in split["folds"].items()}
        else:
            folds = {"all": {subset: len(split.get(subset, [])) for subset in ("train", "test")}}
        splits[split_path.stem] = {"sha256": sha256(split_path), "git_commit": split.get("git_commit"), "generation_parameters": split.get("generation_parameters", {}), "folds": folds}
    qh_summaries: dict[str, Any] = {}
    qh_paths = sorted((root / "result/experiments/quasi_harmonic/raw_runs").glob("*/seed_none/qh_summary.json"))
    qh_paths += sorted((root / "result/experiments/legacy_support/phase2").glob("*/qh_summary.json"))
    qh_paths += sorted((root / "result/experiments/legacy_support/phase2_domains_alloy").glob("*/qh_summary.json"))
    qh_paths += sorted((root / "data/processed").glob("*/qh_remote/qh_summary.json"))
    for qh_path in qh_paths:
        summary = json.loads(qh_path.read_text(encoding="utf-8")); qh_summaries[str(summary.get("system", qh_path.parent.parent.name))] = summary
    qh_diagnostics: dict[tuple[str, str], dict[str, Any]] = {}
    for diagnosis_path in sorted((root / "result/experiments/quasi_harmonic/raw_runs").glob("*/seed_none/*_imaginary_diagnosis.json")):
        diagnosis = json.loads(diagnosis_path.read_text(encoding="utf-8")); key = (str(diagnosis.get("system")), str(diagnosis.get("phase"))); qh_diagnostics[key] = diagnosis
    for phase_key, phase_data in phases.items():
        system, phase = phase_key.split(":", 1)
        summary = qh_summaries.get(system, {})
        qh_phase = summary.get("phases", {}).get(phase)
        if qh_phase is not None:
            minima = qh_phase.get("minimum_frequency_THz_by_volume", {})
            diagnosis = qh_diagnostics.get((system, phase), {})
            fraction = diagnosis.get("negative_mode_fraction")
            phase_data["qh"] = {"supercell": summary.get("supercell_matrix"), "mesh": summary.get("mesh"), "imaginary_fraction": fraction, "minimum_frequency_THz": min(minima.values()) if minima else None, "qh_reliable": None if fraction is None else bool(float(fraction) <= 0.01 and diagnosis.get("negative_modes_are_gamma_local", False))}
    payload = {
        "systems": systems,
        "phases": phases,
        "splits": splits,
        "metric_definitions": "fes_bench.eval.run: evaluator owns G/ΔG MAE, RMSE, sign accuracy, crossings, slopes, and 2σ coverage definitions",
    }
    source = root / "configs" / "raw_runs.json"
    out = root / "result/experiments/data_prep/raw_inventory" / "inventory" / "seed_none" / "metrics.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"provenance": provenance(root, source, "_inventory", "inventory", "none", []), "metrics": payload}, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--manifest", type=Path, default=Path("configs/raw_runs.json"))
    args = parser.parse_args()
    root = args.repo_root.resolve()
    manifest_path = args.manifest if args.manifest.is_absolute() else root / args.manifest
    return normalize(root, json.loads(manifest_path.read_text(encoding="utf-8")))


if __name__ == "__main__":
    raise SystemExit(main())
