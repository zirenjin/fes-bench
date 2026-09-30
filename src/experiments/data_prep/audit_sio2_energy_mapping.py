"""Audit SiO2 composition, bonds, type maps, and same-head energies — plan Table 2b."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--head", default="Domains_SSE_PBE")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    from ase.io import read
    from ase.neighborlist import neighbor_list
    from deepmd.calculator import DP
    import numpy as np

    root = args.repo_root.resolve()
    calculator = DP(model=args.checkpoint, head=args.head)
    rows = []
    for phase in ("quartz_beta", "cristobalite_beta", "tridymite_p63mmc"):
        old = root / "data/processed/sio2/_archive_md_snapshot" / phase / "structure.extxyz"
        new = root / "data/processed/sio2" / phase / "structure.extxyz"
        for label, path in (("old", old), ("new", new)):
            atoms = read(path)
            atoms.calc = calculator
            energy = float(atoms.get_potential_energy())
            distances = neighbor_list("ijd", atoms, cutoff=2.2, self_interaction=False)
            sio = [float(distance) for i, j, distance in zip(*distances) if sorted((atoms[i].symbol, atoms[j].symbol)) == ["O", "Si"]]
            rows.append({
                "phase": phase,
                "version": label,
                "path": str(path.relative_to(root)),
                "n_atoms": len(atoms),
                "species_counts": dict(Counter(atoms.get_chemical_symbols())),
                "si_o_ratio": f"{Counter(atoms.get_chemical_symbols()).get('Si', 0)}:{Counter(atoms.get_chemical_symbols()).get('O', 0)}",
                "nearest_SiO_A": min(sio) if sio else None,
                "mean_four_shortest_SiO_A": float(np.mean(sorted(sio)[:4])) if sio else None,
                "total_energy_eV": energy,
                "energy_eV_per_atom": energy / len(atoms),
            })
    payload = {
        "provenance": {
            "checkpoint_path": "<external-checkpoint>/DPA-3.1-3M.pt",
            "checkpoint_runtime_path": args.checkpoint,
            "head": args.head,
            "checkpoint_type_map_observed": "periodic-table universal map; O and Si identified by symbols (DPA-3.1-3M.pt)",
            "system_type_map": ["O", "Si"],
            "mapping_rule": "DeepMD maps chemical symbols to the checkpoint type map; system.json order is not used to relabel extxyz symbols.",
        },
        "rows": rows,
    }
    output = root / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
