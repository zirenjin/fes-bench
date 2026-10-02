"""Compute DPA-vs-reference-DFT snapshot bias without training a model."""

from __future__ import annotations

import argparse
import csv
import gc
import hashlib
import io
import json
import subprocess
import sys
import tarfile
import tempfile
import zipfile
from pathlib import Path

import numpy as np


PHASES = {
    "hf": ("hcp", "bcc"),
    "ti": ("hcp", "bcc"),
    "zr": ("hcp", "bcc"),
    "sio2": ("quartz_beta", "cristobalite_beta", "tridymite_p63mmc"),
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def parse_outcar(raw: bytes, name: str):
    from ase.io import read

    handle = io.StringIO(raw.decode("utf-8", errors="replace"))
    handle.name = name  # type: ignore[attr-defined]
    try:
        atoms = read(handle, format="vasp-out", index=-1)
    except (StopIteration, ValueError, IndexError):
        return None
    return atoms, float(atoms.get_potential_energy())


def metal_snapshots(path: Path):
    with tarfile.open(path, "r:gz") as archive:
        members = [
            member
            for member in archive.getmembers()
            if member.isfile() and member.name.endswith("OUTCAR") and "/c.training_set_by_high_DFT/" in member.name
        ]
        for member in sorted(members, key=lambda item: item.name):
            handle = archive.extractfile(member)
            if handle is not None:
                record = parse_outcar(handle.read(), f"{path}::{member.name}")
                if record is not None:
                    yield record


def sio2_snapshots(path: Path, phase: str):
    tgz_name = f"direct_upsampling_data_rungs_1-3/moment_tensor_potentials/training_sets/{phase}_pbe_d3bj.tgz"
    with zipfile.ZipFile(path) as outer:
        payload = outer.read(tgz_name)
    with tarfile.open(fileobj=io.BytesIO(payload), mode="r:gz") as archive:
        members = [member for member in archive.getmembers() if member.isfile() and member.name.endswith("OUTCAR")]
        for member in sorted(members, key=lambda item: item.name):
            handle = archive.extractfile(member)
            if handle is not None:
                record = parse_outcar(handle.read(), f"{path}::{tgz_name}::{member.name}")
                if record is not None:
                    yield record


def model_inputs(atoms):
    symbols = atoms.get_chemical_symbols()
    cell = np.asarray(atoms.cell.array, dtype=np.float64)
    if not np.any(cell):
        raise ValueError("snapshot has no periodic cell")
    return np.asarray(atoms.positions, dtype=np.float64).reshape(-1), cell.reshape(-1), symbols


def evaluate_batched(checkpoint: Path, head: str, records, batch_size: int, device: str) -> np.ndarray:
    """Evaluate batches in short-lived workers to bound DeepMD inference memory."""
    buffers: dict[tuple[int, tuple[str, ...]], tuple[list[np.ndarray], list[np.ndarray], list[float]]] = {}
    errors: list[float] = []

    def flush(key: tuple[int, tuple[str, ...]]) -> None:
        coords, cells, dft_values = buffers[key]
        with tempfile.TemporaryDirectory(prefix="dft-static-batch-") as temp_dir:
            batch_file = Path(temp_dir) / "batch.npz"
            np.savez(
                batch_file,
                coords=np.asarray(coords),
                cells=np.asarray(cells),
                dft_values=np.asarray(dft_values),
            )
            command = [
                sys.executable,
                str(Path(__file__).resolve()),
                "--worker",
                "--batch-file",
                str(batch_file),
                "--checkpoint",
                str(checkpoint.resolve()),
                "--head",
                head,
                "--device",
                device,
                "--symbols",
                json.dumps(key[1]),
            ]
            completed = subprocess.run(command, check=False, capture_output=True, text=True)
            if completed.returncode:
                raise RuntimeError(
                    "DPA worker failed with exit code "
                    f"{completed.returncode}: {completed.stderr.strip()}"
                )
            result = json.loads(completed.stdout.strip().splitlines()[-1])
            errors.extend(result["errors_meV_per_atom"])
        buffers[key] = ([], [], [])

    for atoms, dft_energy in records:
        coord, cell, symbols = model_inputs(atoms)
        key = (len(symbols), tuple(symbols))
        coords, cells, dft_values = buffers.setdefault(key, ([], [], []))
        coords.append(coord)
        cells.append(cell)
        dft_values.append(float(dft_energy) / len(symbols))
        if len(coords) >= batch_size:
            flush(key)
    for key, value in list(buffers.items()):
        if value[0]:
            flush(key)
    return np.asarray(errors, dtype=float)


def evaluate_persistent_batched(checkpoint: Path, head: str, records, batch_size: int, device: str) -> np.ndarray:
    """Evaluate batches with one persistent DeepMD model in the worker."""
    from deepmd.infer import DeepPot

    model = DeepPot(str(checkpoint.resolve()), head=head, device=device)
    type_map = list(model.get_type_map())
    buffers: dict[tuple[int, tuple[str, ...]], tuple[list[np.ndarray], list[np.ndarray], list[float]]] = {}
    errors: list[float] = []

    def flush(key: tuple[int, tuple[str, ...]]) -> None:
        coords, cells, dft_values = buffers[key]
        type_ids = np.asarray([type_map.index(symbol) for symbol in key[1]], dtype=np.int32)
        energy, _, _ = model.eval(np.asarray(coords), np.asarray(cells), type_ids)
        n_atoms = key[0]
        errors.extend(
            ((np.asarray(energy).reshape(-1) / n_atoms - np.asarray(dft_values)) * 1000.0).tolist()
        )
        buffers[key] = ([], [], [])

    for atoms, dft_energy in records:
        coord, cell, symbols = model_inputs(atoms)
        key = (len(symbols), tuple(symbols))
        coords, cells, dft_values = buffers.setdefault(key, ([], [], []))
        coords.append(coord)
        cells.append(cell)
        dft_values.append(float(dft_energy) / len(symbols))
        if len(coords) >= batch_size:
            flush(key)
    for key, value in list(buffers.items()):
        if value[0]:
            flush(key)
    return np.asarray(errors, dtype=float)


def evaluate_inprocess(model, records, type_map: list[str]) -> np.ndarray:
    """Evaluate one snapshot at a time with one model per phase."""
    errors: list[float] = []
    for atoms, dft_energy in records:
        coord, cell, symbols = model_inputs(atoms)
        type_ids = np.asarray([type_map.index(symbol) for symbol in symbols], dtype=np.int32)
        energy, _, _ = model.eval(np.asarray([coord]), np.asarray([cell]), type_ids)
        errors.append(float((np.asarray(energy).reshape(-1)[0] / len(symbols) - dft_energy / len(symbols)) * 1000.0))
    return np.asarray(errors, dtype=float)


def worker(args: argparse.Namespace) -> int:
    from deepmd.infer import DeepPot

    payload = np.load(args.batch_file)
    coords = payload["coords"]
    cells = payload["cells"]
    dft_values = payload["dft_values"]
    symbols = json.loads(args.symbols)
    model = DeepPot(str(args.checkpoint.resolve()), head=args.head, device=args.device)
    type_map = list(model.get_type_map())
    type_ids = np.asarray([type_map.index(symbol) for symbol in symbols], dtype=np.int32)
    energy, _, _ = model.eval(coords, cells, type_ids)
    n_atoms = len(symbols)
    errors = ((np.asarray(energy).reshape(-1) / n_atoms - dft_values) * 1000.0).tolist()
    print(json.dumps({"errors_meV_per_atom": errors}))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument(
        "--checkpoint-map",
        type=Path,
        help="JSON object mapping each system to its checkpoint; overrides --checkpoint.",
    )
    parser.add_argument("--system", choices=sorted(PHASES))
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--batch-file", type=Path)
    parser.add_argument("--head")
    parser.add_argument("--symbols")
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--output", type=Path, default=Path("result/experiments/dft_static/dpa_vs_dft_snapshots.csv"))
    parser.add_argument("--head-policy", type=Path, default=Path("configs/models/head_policy.yaml"))
    args = parser.parse_args()
    if args.worker:
        return worker(args)
    root = args.repo_root.resolve()
    checkpoint_map: dict[str, Path] = {}
    if args.checkpoint_map:
        raw_map = json.loads(args.checkpoint_map.read_text(encoding="utf-8"))
        checkpoint_map = {
            system: (Path(path) if Path(path).is_absolute() else root / path)
            for system, path in raw_map.items()
        }
    if not args.checkpoint and not checkpoint_map:
        raise SystemExit("one of --checkpoint or --checkpoint-map is required")
    from deepmd.infer import DeepPot

    policy = json.loads((root / args.head_policy).read_text(encoding="utf-8"))
    phase_rows: list[dict[str, object]] = []
    systems = [args.system] if args.system else list(PHASES)
    for system in systems:
        phases = PHASES[system]
        checkpoint = checkpoint_map.get(system, args.checkpoint)
        if checkpoint is None:
            raise SystemExit(f"no checkpoint supplied for {system}")
        checkpoint = checkpoint.resolve()
        head = policy["systems"][system]["energy_head"]
        for phase in phases:
            if system == "sio2":
                archive = root / "data/raw/sio2/direct_upsampling_data_rungs_1-3.zip"
                records = sio2_snapshots(archive, {"tridymite_p63mmc": "tridymite_beta"}.get(phase, phase)) if archive.is_file() else iter(())
                source = str(archive.relative_to(root)) if archive.is_file() else "missing"
            else:
                archive = root / f"data/raw/{system}/{system.capitalize()}_{phase}_PBE.tar.gz"
                records = metal_snapshots(archive) if archive.is_file() else iter(())
                source = str(archive.relative_to(root)) if archive.is_file() else "missing"
            if args.batch_size == 0:
                model = DeepPot(str(checkpoint), head=head, device=args.device)
                type_map = list(model.get_type_map())
                error = evaluate_inprocess(model, records, type_map)
                del model
            else:
                error = evaluate_persistent_batched(checkpoint, head, records, args.batch_size, args.device)
            gc.collect()
            if not len(error):
                phase_rows.append({"row_type": "phase", "system": system, "phase": phase, "phase_low": "", "phase_high": "", "bias_meV_per_atom": "", "std_meV_per_atom": "", "n_snapshots": 0, "status": "no_available_snapshots", "head": head, "source": source})
                continue
            phase_rows.append({"row_type": "phase", "system": system, "phase": phase, "phase_low": "", "phase_high": "", "bias_meV_per_atom": float(np.mean(error)), "std_meV_per_atom": float(np.std(error, ddof=1)) if len(error) > 1 else 0.0, "n_snapshots": len(error), "status": "ok", "head": head, "source": source})
    pair_rows: list[dict[str, object]] = []
    for system in systems:
        phases = PHASES[system]
        if system == "sio2":
            pairs = (("quartz_beta", "cristobalite_beta"), ("quartz_beta", "tridymite_p63mmc"))
        else:
            pairs = (("hcp", "bcc"),)
        by_phase = {row["phase"]: row for row in phase_rows if row["system"] == system and row["row_type"] == "phase"}
        for low, high in pairs:
            low_row, high_row = by_phase[low], by_phase[high]
            available = low_row["status"] == "ok" and high_row["status"] == "ok"
            pair_rows.append({"row_type": "pair", "system": system, "phase": "", "phase_low": low, "phase_high": high, "bias_meV_per_atom": (float(high_row["bias_meV_per_atom"]) - float(low_row["bias_meV_per_atom"])) if available else "", "std_meV_per_atom": "", "n_snapshots": "", "status": "ok" if available else "no_available_snapshots", "head": high_row["head"], "source": "phase bias difference"})
    rows = phase_rows + pair_rows
    args.output = root / args.output
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = ["row_type", "system", "phase", "phase_low", "phase_high", "bias_meV_per_atom", "std_meV_per_atom", "n_snapshots", "status", "head", "source"]
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    meta = {
        "checkpoint": str(args.checkpoint) if args.checkpoint else None,
        "checkpoint_sha256": sha256(args.checkpoint) if args.checkpoint else None,
        "checkpoint_map": {
            system: {"path": str(path), "sha256": sha256(path)}
            for system, path in checkpoint_map.items()
        },
        "head_policy": str(args.head_policy),
        "energy_definition": "ASE vasp-out get_potential_energy, normally the VASP energy without entropy",
        "systems": systems,
        "status": "complete",
    }
    (args.output.parent / "dpa_vs_dft_snapshots.meta.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(args.output), "phase_rows": len(phase_rows), "pair_rows": len(pair_rows)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
