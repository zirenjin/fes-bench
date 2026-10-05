"""Prepare PES fine-tuning manifests from DFT snapshots — plan PES step 2.

This preparation path is intentionally free-energy blind: it scans only DFT
archive members (OUTCAR/POSCAR/CONTCAR) and never opens ``reference_G.csv`` or
any other thermodynamic table.  The optional extraction stage is kept separate
from the manifest scan so a large archive can be staged on the GPU host first.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import tarfile
from datetime import datetime, timezone
from pathlib import Path


ARCHIVES = {
    "sio2": {"functional": "PBE-D3(BJ)", "source": "data/raw/sio2/direct_upsampling_data_rungs_1-3.zip"},
    "hf": {"functional": "PBE", "source_glob": "data/raw/hf/*.tar.gz"},
    "ti": {"functional": "PBE", "source_glob": "data/raw/ti/*.tar.gz"},
    "zr": {"functional": "PBE", "source_glob": "data/raw/zr/*.tar.gz"},
}

FULL_TYPE_MAP = (
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"
).split()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def split_for_trajectory(key: str) -> str:
    bucket = int(hashlib.sha256(key.encode()).hexdigest()[:2], 16) % 10
    return "train" if bucket < 8 else "validation" if bucket == 8 else "test"


def trajectory_key(member: str) -> str:
    path = member.rsplit("/", 1)[0]
    path = re.sub(r"/vasp_step_[^/]+", "/trajectory", path)
    path = re.sub(r"/step[_-]?\d+", "/trajectory", path)
    return path


def scan_archive(root: Path, system: str, archive: Path) -> dict[str, object]:
    if "reference_G.csv" in str(archive):
        raise AssertionError("PES preparation cannot consume free-energy tables")
    with tarfile.open(archive, "r:gz") as handle:
        members = [member for member in handle.getmembers() if member.isfile() and member.name.endswith("/OUTCAR")]
    trajectories = {}
    for member in members:
        key = trajectory_key(member.name)
        trajectories.setdefault(key, {"trajectory": key, "split": split_for_trajectory(key), "n_outcar": 0, "members": []})
        trajectories[key]["n_outcar"] += 1
        trajectories[key]["members"].append(member.name)
    counts = {name: sum(1 for value in trajectories.values() if value["split"] == name) for name in ("train", "validation", "test")}
    return {
        "system": system,
        "archive": str(archive.relative_to(root)),
        "archive_sha256": sha256(archive),
        "n_outcar": len(members),
        "n_trajectory_groups": len(trajectories),
        "trajectory_group_counts": counts,
        "trajectories": list(trajectories.values()),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--output-root", type=Path, default=Path("data/processed/pes_finetune_v1"))
    parser.add_argument("--manifest-only", action="store_true", help="scan archives without extracting DFT frames")
    args = parser.parse_args()
    root = args.repo_root.resolve()
    output = (root / args.output_root).resolve()
    output.mkdir(parents=True, exist_ok=True)
    archives = []
    for system, spec in ARCHIVES.items():
        if "source_glob" in spec:
            paths = sorted(root.glob(spec["source_glob"]))
        else:
            paths = [root / spec["source"]]
        for archive in paths:
            if not archive.exists():
                archives.append({"system": system, "archive": str(archive.relative_to(root)), "status": "missing"})
                continue
            if archive.suffix == ".zip":
                archives.append({"system": system, "archive": str(archive.relative_to(root)), "status": "zip_requires_darus_extractor"})
            else:
                record = scan_archive(root, system, archive)
                record["functional"] = spec["functional"]
                record["status"] = "scanned_manifest_only" if args.manifest_only else "scanned_manifest_only_extraction_not_enabled"
                archives.append(record)
    manifest = {
        "schema": "fes_bench.pes_finetune_manifest.v1",
        "status": "manifest_only_extraction_pending_remote_staging",
        "free_energy_blind": True,
        "forbidden_inputs": ["data/processed/*/*/reference_G.csv", "data/raw/*/*thermodynamic*.tab", "data/raw/*/thermodynamic_properties/*"],
        "functional_by_system": {system: spec["functional"] for system, spec in ARCHIVES.items()},
        "type_map": FULL_TYPE_MAP,
        "split_rule": "trajectory-level deterministic SHA-256 buckets: 80% train, 10% validation, 10% test; no frame-level random split",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "archives": archives,
        "remote_staging": {"host": "thu-GenSi", "required_path": "/share/jzr/pes_finetune_v1", "status": "blocked_by_full_filesystem_at_last_check"},
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    (output / "README.md").write_text(
        "# PES fine-tuning data manifest\n\n"
        "This directory is a manifest-only checkpoint. Source archives are DFT energy/force/virial data; no free-energy table is opened. Full DeepMD extraction is to be run after the archives are staged in the isolated thu-GenSi directory.\n",
        encoding="utf-8",
    )
    print(output / "manifest.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
