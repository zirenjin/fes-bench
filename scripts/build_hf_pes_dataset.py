"""Extract high-DFT Hf OUTCAR snapshots into DeepMD raw PES systems.

The source archives are read-only DaRUS copies.  This script writes only to
the caller-selected isolated run directory and records that virials are not
present in the source OUTCARs (ISIF=0).
"""

from __future__ import annotations

import argparse
import io
import json
import tarfile
from pathlib import Path

import numpy as np


def read_archive(path: Path, phase: str):
    from ase.io import read

    rows = []
    with tarfile.open(path, "r:gz") as archive:
        members = [m for m in archive.getmembers() if m.isfile() and m.name.endswith("OUTCAR") and "/c.training_set_by_high_DFT/" in m.name]
        for member in sorted(members, key=lambda m: m.name):
            handle = archive.extractfile(member)
            if handle is None:
                continue
            raw = io.StringIO(handle.read().decode("utf-8", errors="replace"))
            # ASE's OUTCAR parser uses fd.name for chunk diagnostics.
            raw.name = member.name  # type: ignore[attr-defined]
            atoms = read(raw, format="vasp-out", index=-1)
            rows.append({
                "phase": phase,
                "source": member.name,
                "energy": float(atoms.get_potential_energy()),
                "force": np.asarray(atoms.get_forces(), dtype=np.float64),
                "coord": np.asarray(atoms.get_positions(), dtype=np.float64),
                "box": np.asarray(atoms.cell.array, dtype=np.float64),
            })
    return rows


def write_system(root: Path, name: str, rows: list[dict], valid_stride: int = 5) -> dict:
    train = [row for i, row in enumerate(rows) if i % valid_stride]
    valid = [row for i, row in enumerate(rows) if i % valid_stride == 0]
    phase_root = root / name
    type_dir = phase_root / "train"
    valid_dir = phase_root / "valid"
    type_dir.mkdir(parents=True, exist_ok=True)
    valid_dir.mkdir(parents=True, exist_ok=True)
    (phase_root / "type.raw").write_text("0\n" * len(rows[0]["force"]), encoding="utf-8")
    (phase_root / "type_map.raw").write_text("Hf\n", encoding="utf-8")
    for folder, subset in ((type_dir, train), (valid_dir, valid)):
        set_dir = folder / "set.000"
        set_dir.mkdir(parents=True, exist_ok=True)
        np.save(set_dir / "energy.npy", np.asarray([r["energy"] for r in subset]))
        np.save(set_dir / "force.npy", np.asarray([r["force"] for r in subset]).reshape(len(subset), -1))
        np.save(set_dir / "coord.npy", np.asarray([r["coord"] for r in subset]).reshape(len(subset), -1))
        np.save(set_dir / "box.npy", np.asarray([r["box"] for r in subset]).reshape(len(subset), -1))
    return {"phase": name, "n_total": len(rows), "n_train": len(train), "n_valid": len(valid), "n_atoms": int(len(rows[0]["force"])), "virial_present": False, "functional": "PBE", "source": "DaRUS Hf *_PBE.tar.gz high-DFT OUTCAR", "source_members": [r["source"] for r in rows]}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hcp-archive", required=True)
    parser.add_argument("--bcc-archive", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    root = Path(args.output)
    root.mkdir(parents=True, exist_ok=True)
    hcp = read_archive(Path(args.hcp_archive), "hcp")
    bcc = read_archive(Path(args.bcc_archive), "bcc")
    manifest = {"protocol": "read final ASE frame from each high-DFT OUTCAR; deterministic 4:1 valid split", "phases": [write_system(root, "hcp", hcp), write_system(root, "bcc", bcc)], "virial_note": "Source OUTCARs show ISIF=0; no stress/virial arrays were synthesized."}
    (root / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"hcp": len(hcp), "bcc": len(bcc), "output": str(root)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
