"""Resume-safe DPA inference for one metal snapshot phase."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from analyze_snapshots import PHASES, metal_snapshots, model_inputs


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--system", choices=("ti", "zr"), required=True)
    parser.add_argument("--phase", choices=("hcp", "bcc"), required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--chunk-size", type=int, default=16)
    parser.add_argument(
        "--start-index",
        type=int,
        default=0,
        help="Skip this many valid snapshots before evaluating new chunks.",
    )
    parser.add_argument(
        "--chunk-index-start",
        type=int,
        default=0,
        help="Index for the first newly written chunk.",
    )
    parser.add_argument("--device", default="cpu")
    args = parser.parse_args()
    if args.chunk_size < 1:
        raise SystemExit("--chunk-size must be positive")

    archive = args.repo_root / f"data/raw/{args.system}/{args.system.capitalize()}_{args.phase}_PBE.tar.gz"
    args.output_dir.mkdir(parents=True, exist_ok=True)
    records = metal_snapshots(archive)
    model = None
    type_map = None
    chunk_index = args.chunk_index_start
    record_index = 0
    pending_coords: list[np.ndarray] = []
    pending_cells: list[np.ndarray] = []
    pending_dft: list[float] = []
    pending_symbols: tuple[str, ...] | None = None

    def flush() -> None:
        nonlocal model, type_map, chunk_index, pending_coords, pending_cells, pending_dft, pending_symbols
        if not pending_coords:
            return
        target = args.output_dir / f"chunk_{chunk_index:05d}.json"
        if not target.is_file():
            if model is None:
                from deepmd.infer import DeepPot

                model = DeepPot(str(args.checkpoint.resolve()), head="Domains_Alloy", device=args.device)
                type_map = list(model.get_type_map())
            type_ids = np.asarray([type_map.index(symbol) for symbol in pending_symbols], dtype=np.int32)
            energy, _, _ = model.eval(
                np.asarray(pending_coords), np.asarray(pending_cells), type_ids
            )
            n_atoms = len(pending_symbols)
            errors = (
                np.asarray(energy).reshape(-1) / n_atoms - np.asarray(pending_dft)
            ) * 1000.0
            target.write_text(
                json.dumps({
                    "system": args.system,
                    "phase": args.phase,
                    "chunk": chunk_index,
                    "n_snapshots": len(errors),
                    "errors_meV_per_atom": errors.tolist(),
                }) + "\n",
                encoding="utf-8",
            )
        chunk_index += 1
        pending_coords, pending_cells, pending_dft = [], [], []
        pending_symbols = None

    for atoms, dft_energy in records:
        if record_index < args.start_index:
            record_index += 1
            continue
        record_index += 1
        coord, cell, symbols = model_inputs(atoms)
        symbol_tuple = tuple(symbols)
        if pending_symbols is not None and symbol_tuple != pending_symbols:
            flush()
        pending_symbols = symbol_tuple
        pending_coords.append(coord)
        pending_cells.append(cell)
        pending_dft.append(float(dft_energy) / len(symbols))
        if len(pending_coords) >= args.chunk_size:
            flush()
    flush()
    all_errors: list[float] = []
    for path in sorted(args.output_dir.glob("chunk_*.json")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        all_errors.extend(float(value) for value in payload["errors_meV_per_atom"])
    if not all_errors:
        raise SystemExit("no valid snapshots found")
    result = {
        "system": args.system,
        "phase": args.phase,
        "status": "ok",
        "bias_meV_per_atom": float(np.mean(all_errors)),
        "std_meV_per_atom": float(np.std(all_errors, ddof=1)) if len(all_errors) > 1 else 0.0,
        "n_snapshots": len(all_errors),
        "head": "Domains_Alloy",
        "source": str(archive.relative_to(args.repo_root)),
    }
    (args.output_dir / "phase_result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
