"""Build leakage-checked DeepMD FES training data from canonical processed inputs.

The builder is deliberately the only path from ``data/processed`` into T3
training directories.  It reads frozen split indices (never creates a split),
repeats each canonical representative structure at the requested reference
temperatures, and records hashes and head policy in a provenance file.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


SYSTEM_PHASES = {
    "sio2": ("quartz_beta", "cristobalite_beta", "tridymite_p63mmc"),
    "hf": ("hcp", "bcc"),
    "ti": ("hcp", "bcc"),
    "zr": ("hcp", "bcc"),
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_config(path: Path) -> dict[str, Any]:
    # Repository configs use JSON-compatible YAML.  Keep the script usable in
    # the minimal builder environment; PyYAML is optional.
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        try:
            import yaml  # type: ignore
        except ImportError as exc:  # pragma: no cover - only for non-JSON YAML
            raise RuntimeError(f"cannot parse YAML without PyYAML: {path}") from exc
        return yaml.safe_load(path.read_text(encoding="utf-8"))


def parse_extxyz(path: Path) -> tuple[list[str], np.ndarray, np.ndarray]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or int(lines[0].strip()) <= 0:
        raise ValueError(f"invalid extxyz atom count: {path}")
    natoms = int(lines[0].strip())
    if len(lines) < natoms + 2:
        raise ValueError(f"truncated extxyz: {path}")
    lattice = re.search(r'Lattice="([^"]+)"', lines[1])
    if lattice is None:
        raise ValueError(f"extxyz has no Lattice attribute: {path}")
    cell = np.asarray([float(x) for x in lattice.group(1).split()], dtype=float)
    if cell.size != 9:
        raise ValueError(f"expected 9 lattice values: {path}")
    species: list[str] = []
    coords: list[list[float]] = []
    for line in lines[2 : 2 + natoms]:
        fields = line.split()
        if len(fields) < 4:
            raise ValueError(f"invalid atom row in {path}: {line!r}")
        species.append(fields[0])
        coords.append([float(value) for value in fields[1:4]])
    if len(species) != natoms:
        raise ValueError(f"atom count mismatch in {path}")
    return species, np.asarray(coords, dtype=float), cell.reshape(3, 3)


def reference_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    required = {"T_K", "P_GPa", "G_eV_per_atom"}
    if not rows or not required <= set(rows[0]):
        raise ValueError(f"reference_G.csv missing required columns: {path}")
    return rows


def fold_indices(split: dict[str, Any]) -> dict[str, dict[str, list[dict[str, Any]]]]:
    if "folds" not in split:
        return {"all": {"train": list(split["train"]), "test": list(split["test"])}}
    return {
        str(name): {"train": list(value["train"]), "test": list(value["test"])}
        for name, value in split["folds"].items()
    }


def model_type_map(species: list[str], system: str) -> list[str]:
    unique = set(species)
    if system == "sio2":
        expected = {"Si", "O"}
        if unique != expected:
            raise AssertionError(f"SiO2 representative must contain Si/O, got {sorted(unique)}")
        # E1 checkpoints use this order.  system.json historically lists O,Si;
        # type.raw is generated against the checkpoint order, not that metadata.
        return ["Si", "O"]
    if len(unique) != 1:
        raise AssertionError(f"metal phase contains unexpected species: {sorted(unique)}")
    return sorted(unique)


def build_phase(root: Path, output: Path, system: str, phase: str,
                selected: list[dict[str, Any]], head: str) -> dict[str, Any]:
    processed = root / "data/processed" / system / phase
    structure = processed / "structure.extxyz"
    ref_path = processed / "reference_G.csv"
    species, coords, cell = parse_extxyz(structure)
    type_map = model_type_map(species, system)
    type_ids = np.asarray([type_map.index(symbol) for symbol in species], dtype=np.int32)
    rows = reference_rows(ref_path)
    by_index = {index: row for index, row in enumerate(rows)}
    indexes = sorted({int(row["T_index"]) for row in selected})
    missing = [index for index in indexes if index not in by_index]
    if missing:
        raise IndexError(f"{system}/{phase}: T_index out of range: {missing[:5]}")
    temperatures = np.asarray([float(by_index[index]["T_K"]) for index in indexes])
    pressures = np.asarray([float(by_index[index]["P_GPa"]) for index in indexes])
    energies = np.asarray([float(by_index[index]["G_eV_per_atom"]) for index in indexes])
    output.mkdir(parents=True, exist_ok=True)
    (output / "type.raw").write_text("\n".join(str(int(value)) for value in type_ids) + "\n", encoding="utf-8")
    set_path = output / "set.000"
    set_path.mkdir(parents=True, exist_ok=True)
    nframes = len(indexes)
    np.save(set_path / "coord.npy", np.tile(coords.reshape(1, -1), (nframes, 1)))
    np.save(set_path / "box.npy", np.tile(cell.reshape(1, -1), (nframes, 1)))
    np.save(set_path / "fparam.npy", np.column_stack([temperatures, pressures]))
    np.save(set_path / "free_energy.npy", energies.reshape(-1, 1))
    return {
        "system": system,
        "phase": phase,
        "n_atoms": len(species),
        "n_frames": nframes,
        "train_T_indices": indexes,
        "structure_path": str(structure.relative_to(root)),
        "structure_sha256": sha256(structure),
        "reference_path": str(ref_path.relative_to(root)),
        "reference_sha256": sha256(ref_path),
        "type_map": type_map,
        "head": head,
    }


def build_fold(root: Path, split_path: Path, split_name: str, fold: str,
               data: dict[str, Any], output_root: Path, systems: list[str],
               policy: dict[str, Any]) -> Path:
    folds = fold_indices(data)
    if fold not in folds:
        raise KeyError(f"unknown {split_name} fold {fold}; available={sorted(folds)}")
    train = folds[fold]["train"]
    test = folds[fold]["test"]
    train_keys = {(r["system"], r["phase"], int(r["T_index"])) for r in train}
    test_keys = {(r["system"], r["phase"], int(r["T_index"])) for r in test}
    overlap = train_keys & test_keys
    if overlap:
        raise AssertionError(f"{split_name}/{fold}: train/test overlap: {sorted(overlap)[:3]}")
    destination = output_root / split_name / fold
    phases: list[dict[str, Any]] = []
    for system in systems:
        head = policy["systems"][system]["energy_head"]
        selected = [row for row in train if row["system"] == system]
        by_phase: dict[str, list[dict[str, Any]]] = {}
        for row in selected:
            by_phase.setdefault(str(row["phase"]), []).append(row)
        for phase, rows in sorted(by_phase.items()):
            info = build_phase(root, destination / system / phase, system, phase, rows, head)
            phase_keys = {(system, phase, index) for index in info["train_T_indices"]}
            if phase_keys & test_keys:
                raise AssertionError(f"{split_name}/{fold}/{system}/{phase}: test frame emitted")
            phases.append(info)
    provenance = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "split": split_name,
        "fold": fold,
        "split_path": str(split_path.relative_to(root)),
        "split_sha256": sha256(split_path),
        "head_policy": "configs/models/head_policy.yaml",
        "head_policy_sha256": sha256(root / "configs/models/head_policy.yaml"),
        "systems": systems,
        "train_frame_count": len(train),
        "test_frame_count": len(test),
        "train_test_overlap": [],
        "phases": phases,
    }
    (destination / "provenance.json").parent.mkdir(parents=True, exist_ok=True)
    (destination / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    return destination


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--split", type=Path, action="append", required=True,
                        help="frozen split JSON; may be repeated")
    parser.add_argument("--output-root", type=Path, default=Path("result/experiments/t3_training_data_v2"))
    parser.add_argument("--systems", nargs="+", choices=sorted(SYSTEM_PHASES), default=sorted(SYSTEM_PHASES))
    args = parser.parse_args()
    root = args.repo_root.resolve()
    policy = read_config(root / "configs/models/head_policy.yaml")
    output_root = args.output_root if args.output_root.is_absolute() else root / args.output_root
    for split_arg in args.split:
        split_path = split_arg if split_arg.is_absolute() else root / split_arg
        data = read_config(split_path)
        split_name = str(data.get("name") or split_path.stem)
        for fold in fold_indices(data):
            build_fold(root, split_path, split_name, fold, data, output_root, args.systems, policy)
    print(output_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
