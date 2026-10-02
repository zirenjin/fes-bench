"""Parse VASP outputs into compact static-energy evidence."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def last_float(lines: list[str], pattern: str) -> float | None:
    values = []
    regex = re.compile(pattern)
    for line in lines:
        match = regex.search(line)
        if match:
            values.append(float(match.group(1)))
    return values[-1] if values else None


def symmetry_record(path: Path | None, expected: str | None) -> dict[str, object]:
    result = {
        "structure": str(path) if path else "",
        "structure_sha256": sha256(path) if path and path.is_file() else "",
        "spacegroup_number": None,
        "spacegroup_symbol": "",
        "expected_spacegroup": expected or "",
        "spacegroup_matches": "not_checked",
    }
    if path is None or not path.is_file():
        return result
    from ase.io import read
    import spglib

    atoms = read(path, format="vasp", index=-1)
    dataset = spglib.get_symmetry_dataset(
        (atoms.cell.array, atoms.get_scaled_positions(), atoms.numbers),
        symprec=1.0e-3,
    )
    if dataset is None:
        result["spacegroup_matches"] = "symmetry_not_found"
        return result
    result["spacegroup_number"] = int(dataset["number"])
    result["spacegroup_symbol"] = str(dataset["international"])
    result["spacegroup_matches"] = (
        "yes" if expected and result["spacegroup_symbol"] == expected else
        "no" if expected else "not_checked"
    )
    return result


def parse_outcar(path: Path, structure: Path | None = None, expected: str | None = None) -> dict[str, object]:
    text = path.read_text(errors="replace")
    lines = text.splitlines()
    energy = last_float(lines, r"free\s+energy\s+TOTEN\s*=\s*([-+0-9.Ee]+)")
    energy_zero = last_float(lines, r"energy\s+without\s+entropy\s*=\s*([-+0-9.Ee]+)")
    nions = last_float(lines, r"NIONS\s*=\s*([0-9]+)")
    if energy is None and energy_zero is None:
        raise ValueError(f"no VASP energy found in {path}")
    result = {
        "outcar": str(path),
        "outcar_sha256": sha256(path),
        "energy_free_eV": energy,
        "energy_without_entropy_eV": energy_zero,
        "energy_eV": energy_zero if energy_zero is not None else energy,
        "n_atoms": int(nions) if nions is not None else None,
        "converged": "aborting loop because EDIFF is reached" in text or "reached required accuracy" in text,
        "vasp_version": next((line.strip() for line in lines if line.startswith(" vasp.")), "unknown"),
    }
    result.update(symmetry_record(structure, expected))
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outcar", type=Path, action="append", required=True)
    parser.add_argument("--structure", type=Path, action="append")
    parser.add_argument("--expected-spacegroup", action="append")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    structures = args.structure or []
    expected = args.expected_spacegroup or []
    rows = [
        parse_outcar(
            path,
            structures[index] if index < len(structures) else None,
            expected[index] if index < len(expected) else None,
        )
        for index, path in enumerate(args.outcar)
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({"parsed": len(rows), "output": str(args.output)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
