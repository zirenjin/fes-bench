"""Promote validated ideal SiO2 representatives and archive MD snapshots — plan Table 2b."""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    args = parser.parse_args()
    root = args.repo_root.resolve(); sio2 = root / "data/processed/sio2"
    candidates = sio2 / "_rebuild_candidates"; archive = sio2 / "_archive_md_snapshot"
    phases = ["quartz_beta", "cristobalite_beta", "tridymite_p63mmc"]
    for phase in phases:
        current = sio2 / phase; source = candidates / phase
        target_archive = archive / phase; target_archive.mkdir(parents=True, exist_ok=True)
        for name in ("structure.source.extxyz", "structure.extxyz", "meta.json", "representative_relaxation.log"):
            old = current / name
            if old.exists() and not (target_archive / name).exists():
                shutil.copy2(old, target_archive / name)
        for name in ("structure.source.extxyz", "structure.extxyz"):
            shutil.copy2(source / name, current / name)
        old_meta = json.loads((current / "meta.json").read_text(encoding="utf-8"))
        new_meta = json.loads((source / "meta.json").read_text(encoding="utf-8"))
        for label, path in (("source_atoms", source / "structure.source.extxyz"), ("relaxed_atoms", source / "structure.extxyz")):
            try:
                new_meta[label] = int(path.read_text(encoding="utf-8").splitlines()[0].strip())
            except (OSError, ValueError, IndexError):
                new_meta[label] = None
        old_meta["representative"] = new_meta
        old_meta["representative"]["archive_previous"] = str((target_archive / "meta.json").relative_to(root))
        old_meta["downstream_recompute_required"] = [
            "result/experiments/external_baselines (Bartel and phase-level baselines)",
            "result/experiments/quasi_harmonic",
            "result/experiments/imaginary_modes",
            "result/experiments/reference_statistics",
            "result/experiments/representative_structures (E4)",
        ]
        (current / "meta.json").write_text(json.dumps(old_meta, indent=2) + "\n", encoding="utf-8")
        print(f"promoted {phase}; archived {target_archive}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
