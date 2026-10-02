"""Generate provenance-bearing VASP inputs for the dft_static experiment."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

from ase.io import read, write


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


def read_settings(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    return {row["system"]: row for row in rows}


def variants(row: dict[str, str]) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in row["potcar_titles"].split(";"):
        element, title = item.split("=", 1)
        match = re.search(r"PAW_PBE\s+([^\s]+)", title)
        if not match:
            raise ValueError(f"cannot derive POTCAR directory from {title!r}")
        result[element.strip()] = match.group(1)
    return result


def write_incar(path: Path, row: dict[str, str], mode: str, encut: float, mesh: tuple[int, int, int]) -> None:
    lines = [
        "SYSTEM = fes-bench dft_static",
        "GGA = PE",
        f"ENCUT = {encut:.8f}",
        f"PREC = {row['prec']}",
        "EDIFF = 1E-8",
        f"ISMEAR = {row['ismear']}",
        f"SIGMA = {row['sigma_eV']}",
        f"ISPIN = {row['ispin']}",
        f"LASPH = {row['lasph'].upper()}",
        f"ADDGRID = {row['addgrid'].upper()}",
        "ISYM = 2",
    ]
    if row["ivdw"]:
        lines.append(f"IVDW = {row['ivdw']}")
    if mode == "relax":
        lines.extend(["IBRION = 2", "NSW = 200", "ISIF = 3", "EDIFFG = -0.001"])
    else:
        lines.extend(["IBRION = -1", "NSW = 0", "ISIF = 2"])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_kpoints(path: Path, mesh: tuple[int, int, int]) -> None:
    path.write_text(
        "fes-bench dft_static\n0\nGamma\n%s %s %s\n0 0 0\n" % mesh,
        encoding="utf-8",
    )


def write_potcar(path: Path, atoms, potcar_dir: Path, row: dict[str, str]) -> list[dict[str, str]]:
    mapping = variants(row)
    # ASE's VASP writer with sort=True uses alphabetical species order.
    # Use the same order for POTCAR, including structures whose source atom
    # order is interleaved (for example O-Si-O).
    ordered = sorted(set(atoms.get_chemical_symbols()))
    records = []
    with path.open("wb") as output:
        for symbol in ordered:
            if symbol not in mapping:
                raise ValueError(f"no reference POTCAR title for {symbol} in {row['system']}")
            source = potcar_dir / mapping[symbol] / "POTCAR"
            if not source.is_file():
                raise FileNotFoundError(source)
            title = next((x.strip() for x in row["potcar_titles"].split(";") if x.strip().startswith(symbol + "=")), "")
            source_title = title.split("=", 1)[1]
            first_titel = next(
                (line.strip() for line in source.read_text(errors="replace").splitlines() if "TITEL" in line),
                "",
            )
            if source_title not in first_titel:
                raise ValueError(f"POTCAR TITEL mismatch for {symbol}: expected {source_title!r}, got {first_titel!r}")
            output.write(source.read_bytes())
            records.append({"element": symbol, "source": str(source), "expected_titel": source_title, "source_titel": first_titel})
    return records


def generate(root: Path, output: Path, potcar_dir: Path, system: str, phase: str, mode: str, mesh: tuple[int, int, int], encut_factor: float) -> dict:
    settings = read_settings(root / "result/experiments/dft_static/reference_settings.csv")[system]
    atoms = read(root / "data/processed" / system / phase / "structure.extxyz")
    output.mkdir(parents=True, exist_ok=True)
    # POSCAR groups atoms by species; keep this consistent with the concatenated POTCAR.
    write(output / "POSCAR", atoms, format="vasp", direct=True, sort=True, vasp5=True)
    encut = float(settings["encut_eV"]) * encut_factor
    write_incar(output / "INCAR", settings, mode, encut, mesh)
    write_kpoints(output / "KPOINTS", mesh)
    potcar_records = write_potcar(output / "POTCAR", atoms, potcar_dir, settings)
    manifest = {
        "system": system,
        "phase": phase,
        "mode": mode,
        "mesh": list(mesh),
        "encut_eV": encut,
        "reference_settings": settings,
        "structure": str((root / "data/processed" / system / phase / "structure.extxyz").relative_to(root)),
        "structure_sha256": sha256(root / "data/processed" / system / phase / "structure.extxyz"),
        "potcar": potcar_records,
        "status": "inputs_generated",
    }
    (output / "input_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--potcar-dir", type=Path, required=True)
    parser.add_argument("--system", choices=sorted(PHASES), action="append")
    parser.add_argument("--phase", action="append")
    parser.add_argument("--mode", choices=("single_point", "relax", "final_static"), default="single_point")
    parser.add_argument("--mesh", nargs=3, type=int, default=(1, 1, 1))
    parser.add_argument("--encut-factor", type=float, default=1.0)
    args = parser.parse_args()
    root = args.repo_root.resolve()
    systems = args.system or sorted(PHASES)
    phases = args.phase or [phase for system in systems for phase in PHASES[system]]
    selected = [(system, phase) for system in systems for phase in PHASES[system] if phase in phases]
    rows = []
    for system, phase in selected:
        target = args.output_root / system / phase / args.mode / ("k%s%s%s_encut%.3f" % (*args.mesh, args.encut_factor))
        rows.append(generate(root, target, args.potcar_dir.resolve(), system, phase, args.mode, tuple(args.mesh), args.encut_factor))
    (args.output_root / "input_manifest.json").parent.mkdir(parents=True, exist_ok=True)
    (args.output_root / "input_manifest.json").write_text(json.dumps({"rows": rows}, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"generated": len(rows), "output_root": str(args.output_root)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
