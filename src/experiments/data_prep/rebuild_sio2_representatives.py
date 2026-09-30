"""Build and symmetry-constrained-relax ideal SiO2 representatives — plan Table 2b.

The generator reads the AFLOW prototype CIFs recorded in the rebuild config,
writes candidate extxyz files and provenance, and (when ``--relax`` is used)
applies ASE FixSymmetry with the configured Domains_SSE_PBE head.  The default
relaxation is hydrostatic variable-cell, so the DPA equilibrium volume is not
confused with a finite-temperature DaRUS snapshot volume.  It never trains a
model and refuses to promote a relaxation whose detected space group does not
equal the requested target.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def space_group(atoms) -> str | None:
    import spglib

    dataset = spglib.get_symmetry_dataset((atoms.cell.array, atoms.get_scaled_positions(), atoms.numbers), symprec=1e-3)
    if dataset is None:
        return None
    return str(getattr(dataset, "international", dataset["international"]))


def run(config_path: Path, *, relax: bool, checkpoint: str | None = None) -> int:
    from ase.io import read, write

    cfg = json.loads(config_path.read_text(encoding="utf-8"))
    root = config_path.parents[2]
    source_root = (root / cfg["source_root"]).resolve()
    output_root = (root / cfg["output_root"]).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    for phase, item in cfg["phases"].items():
        source = source_root / item["source"]
        out = output_root / phase
        out.mkdir(parents=True, exist_ok=True)
        atoms = read(source)
        before = space_group(atoms)
        if before != item["target_space_group"]:
            raise RuntimeError(f"{phase}: source space group {before!r} != target {item['target_space_group']!r}")
        write(out / "structure.source.extxyz", atoms, format="extxyz")
        meta = {
            "status": "ideal_source_generated",
            "source": str(source.relative_to(root)),
            "source_sha256": digest(source),
            "source_url": item["source_url"],
            "prototype": item["prototype"],
            "icsd": item.get("icsd"),
            "darus_phase": item.get("darus_phase"),
            "reference_model_audit": item.get("reference_model"),
            "reference_paper": "Forslund et al., npj Computational Materials 12, 14 (2026), DOI:10.1038/s41524-025-01874-1; DaRUS-4999",
            "space_group": {"before": before},
            "target_space_group": item["target_space_group"],
            "calculator": "deepmd.calculator.DP" if relax else None,
            "head": cfg["head"] if relax else None,
        }
        if relax:
            from ase.constraints import FixSymmetry
            from ase.filters import FrechetCellFilter
            from ase.optimize import FIRE
            from deepmd.calculator import DP

            atoms.set_constraint(FixSymmetry(atoms, symprec=1e-3))
            atoms.calc = DP(model=checkpoint or cfg["checkpoint"], head=cfg["head"])
            before_energy = atoms.get_potential_energy() / len(atoms)
            before_positions = atoms.get_positions().copy()
            dynamics = FrechetCellFilter(atoms, scalar_pressure=0.0, hydrostatic_strain=True) if bool(cfg.get("relax_cell", True)) else atoms
            opt = FIRE(dynamics, logfile=str(out / "relaxation.log"))
            opt.run(fmax=float(cfg["fmax_eV_per_A"]), steps=int(cfg["max_steps"]))
            after = space_group(atoms)
            if after != item["target_space_group"]:
                raise RuntimeError(f"{phase}: relaxed space group {after!r} != target {item['target_space_group']!r}")
            write(out / "structure.extxyz", atoms, format="extxyz")
            meta.update({"status": "relaxed", "converged": bool(opt.converged()), "relaxation_mode": "hydrostatic_variable_cell" if bool(cfg.get("relax_cell", True)) else "fixed_cell", "space_group": {"before": before, "after": after}, "energy_eV_per_atom": {"before": before_energy, "after": atoms.get_potential_energy() / len(atoms)}, "max_displacement_A": float((((atoms.get_positions() - before_positions) ** 2).sum(axis=1) ** 0.5).max()), "fmax_eV_per_A": float(cfg["fmax_eV_per_A"]), "max_steps": int(cfg["max_steps"])})
        (out / "meta.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
        print(phase, meta["status"], before, meta.get("space_group", {}).get("after"))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("configs/representatives/sio2_ideal_rebuild.yaml"))
    parser.add_argument("--relax", action="store_true")
    parser.add_argument("--checkpoint")
    args = parser.parse_args()
    return run(args.config.resolve(), relax=args.relax, checkpoint=args.checkpoint)


if __name__ == "__main__":
    raise SystemExit(main())
