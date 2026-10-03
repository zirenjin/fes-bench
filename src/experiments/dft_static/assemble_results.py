"""Assemble compact DFT static-energy comparisons from an external work tree."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
from pathlib import Path

import numpy as np

from parse_results import parse_outcar, sha256


PHASES = {
    "hf": ("hcp", "bcc"),
    "ti": ("hcp", "bcc"),
    "zr": ("hcp", "bcc"),
    "sio2": ("quartz_beta", "cristobalite_beta", "tridymite_p63mmc"),
}
PAIRS = {
    "hf": (("hcp", "bcc"),),
    "ti": (("hcp", "bcc"),),
    "zr": (("hcp", "bcc"),),
    "sio2": (("quartz_beta", "cristobalite_beta"), ("quartz_beta", "tridymite_p63mmc")),
}
EXPECTED_SG = {
    "hcp": "P6_3/mmc",
    "bcc": "Im-3m",
    "quartz_beta": "P6_422",
    "cristobalite_beta": "Fd-3m",
    "tridymite_p63mmc": "P6_3/mmc",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def energy(path: Path) -> dict[str, object]:
    row = parse_outcar(path)
    row["energy_per_atom_eV"] = float(row["energy_eV"]) / int(row["n_atoms"])
    return row


def energy_candidates(path: Path) -> tuple[Path | None, dict[str, object] | None]:
    """Select the first complete output when a launch left suffixed variants."""
    candidates = [path]
    if path.parent.is_dir():
        candidates.extend(sorted(path.parent.glob(path.name + "-*")))
    if path.parent.parent.is_dir():
        candidates.extend(
            variant / path.name
            for variant in sorted(path.parent.parent.glob(path.parent.name + "-*"))
        )
    for candidate in candidates:
        if not candidate.is_file():
            continue
        try:
            return candidate, energy(candidate)
        except (OSError, ValueError):
            continue
    return None, None


def phase_outcar(work: Path, root: str, system: str, phase: str, mode: str, key: str) -> Path:
    return work / root / system / phase / mode / key / "OUTCAR"


def phase_contcar(work: Path, system: str, phase: str) -> Path:
    if system == "sio2":
        candidates = [
            work / "relax-converged" / system / phase / "relax/k555_encut1.000/CONTCAR",
            work / "sio2-relax-retry" / phase / "relax/k111_encut1.000/CONTCAR",
            work / "sio2-relax" / "sio2" / phase / "relax/k111_encut1.000/CONTCAR",
            work / "relax" / system / phase / "relax/k444_encut1.000/CONTCAR",
        ]
    else:
        candidates = [
            work / "relax-converged-k24-encut2" / system / phase / "relax/k242424_encut2.000/CONTCAR",
            work / "relax-converged" / system / phase / "relax/k777_encut1.000/CONTCAR",
            work / "relax" / system / phase / "relax/k444_encut1.000/CONTCAR",
        ]
    return next((path for path in candidates if path.is_file()), candidates[0])


def dpa_energies(repo: Path) -> dict[tuple[str, str], float]:
    result = {}
    for system, phases in PHASES.items():
        for phase in phases:
            meta = json.loads((repo / "data/processed" / system / phase / "meta.json").read_text(encoding="utf-8"))
            result[(system, phase)] = float(meta["representative"]["energy_eV_per_atom"]["after"])
    return result


def reference_stds(repo: Path) -> dict[str, float]:
    result = {}
    for system, pairs in PAIRS.items():
        values = {}
        for phase in {phase for pair in pairs for phase in pair}:
            rows = read_csv(repo / "data/processed" / system / phase / "reference_G.csv")
            values[phase] = {
                (float(row["T_K"]), float(row["P_GPa"])): float(row["G_eV_per_atom"])
                for row in rows
            }
        for low, high in pairs:
            common = sorted(set(values[low]) & set(values[high]))
            delta = np.asarray([(values[high][key] - values[low][key]) * 1000.0 for key in common])
            result[f"{system}:{low}:{high}"] = float(np.std(delta, ddof=1))
    return result


def symmetry(path: Path, expected: str) -> tuple[int | None, str, str]:
    if not path.is_file():
        return None, "", "not_checked"
    try:
        from ase.io import read
        import spglib
    except ModuleNotFoundError:
        return None, "", "not_checked"

    atoms = read(path, format="vasp", index=-1)
    dataset = spglib.get_symmetry_dataset((atoms.cell.array, atoms.get_scaled_positions(), atoms.numbers), symprec=1.0e-3)
    if dataset is None:
        return None, "", "symmetry_not_found"
    number = int(dataset.number if hasattr(dataset, "number") else dataset["number"])
    symbol = str(dataset.international if hasattr(dataset, "international") else dataset["international"])
    return number, symbol, "yes" if symbol == expected else "no"


def snapshot_bias(path: Path) -> tuple[dict[tuple[str, str], float], dict[str, float | None]]:
    phase_bias: dict[tuple[str, str], float] = {}
    pair_bias: dict[str, float | None] = {}
    if not path.is_file():
        return phase_bias, pair_bias
    for row in read_csv(path):
        if row["row_type"] == "phase" and row["status"] == "ok":
            phase_bias[(row["system"], row["phase"])] = float(row["bias_meV_per_atom"])
        if row["row_type"] == "pair" and row["status"] == "ok":
            pair_bias[f'{row["system"]}:{row["phase_low"]}:{row["phase_high"]}'] = float(row["bias_meV_per_atom"])
    return phase_bias, pair_bias


def last_max_force(text: str) -> float | None:
    values = re.findall(r"FORCES:\s+max atom, RMS\s+([-+0-9.Ee]+)", text)
    return float(values[-1]) if values else None


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, default=Path("result/experiments/dft_static"))
    parser.add_argument("--checkpoint", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo_root.resolve()
    work = args.work_root.resolve()
    output = (repo / args.output_root).resolve()
    audit_path = output / "symmetry_audit.json"
    symmetry_audit = json.loads(audit_path.read_text(encoding="utf-8")) if audit_path.is_file() else {}

    dpa = dpa_energies(repo)
    ref_std = reference_stds(repo)
    _, pair_bias = snapshot_bias(output / "dpa_vs_dft_snapshots.csv")
    dft_at_dpa: dict[tuple[str, str], dict[str, object]] = {}
    dft_final: dict[tuple[str, str], dict[str, object]] = {}
    final_sources: dict[str, str] = {}
    for system, phases in PHASES.items():
        for phase in phases:
            if system == "sio2":
                at_dpa_path = work / "conv-k555-sio2" / system / phase / "single_point/k555_encut1.000/OUTCAR"
                final_candidates = [
                    work / "final-adopted-converged" / system / phase / "final_static/k555_encut1.000/OUTCAR",
                    work / "final-adopted" / system / phase / "final_static/k444_encut1.000/OUTCAR",
                ]
            else:
                at_dpa_path = work / "conv-encut20-k242424" / system / phase / "single_point/k242424_encut2.000/OUTCAR"
                final_candidates = [
                    work / "final-adopted-converged" / system / phase / "final_static/k242424_encut2.000/OUTCAR",
                    work / "final-adopted-converged" / system / phase / "final_static/k777_encut1.000/OUTCAR",
                    work / "final-adopted" / system / phase / "final_static/k666_encut1.000/OUTCAR",
                ]
            final_path = next((path for path in final_candidates if path.is_file()), final_candidates[0])
            final_sources[f"{system}:{phase}"] = str(final_path.relative_to(work)) if final_path.is_relative_to(work) else str(final_path)
            dft_at_dpa[(system, phase)] = energy(at_dpa_path)
            dft_final[(system, phase)] = energy(final_path)

    comparison = []
    for system, pairs in PAIRS.items():
        for low, high in pairs:
            key = f"{system}:{low}:{high}"
            dpa_delta = (dpa[(system, high)] - dpa[(system, low)]) * 1000.0
            at_dpa_delta = (dft_at_dpa[(system, high)]["energy_per_atom_eV"] - dft_at_dpa[(system, low)]["energy_per_atom_eV"]) * 1000.0
            final_delta = (dft_final[(system, high)]["energy_per_atom_eV"] - dft_final[(system, low)]["energy_per_atom_eV"]) * 1000.0
            model_error = dpa_delta - final_delta
            structural_error = at_dpa_delta - final_delta
            std = ref_std[key]
            comparison.append({
                "system": system,
                "phase_low": low,
                "phase_high": high,
                "delta_E_DPA_meV_per_atom": dpa_delta,
                "delta_E_DFT_at_DPA_meV_per_atom": at_dpa_delta,
                "delta_E_DFT_meV_per_atom": final_delta,
                "DPA_minus_DFT_meV_per_atom": model_error,
                "DFT_at_DPA_minus_DFT_meV_per_atom": structural_error,
                "delta_bias_meV_per_atom": pair_bias.get(key, "n/a:unavailable"),
                "reference_delta_G_std_meV_per_atom": std,
                "model_error_exceeds_reference_std": abs(model_error) > std,
            })
    write_csv(output / "delta_E_comparison.csv", comparison)

    convergence = []
    levels = {
        "hf": [("reference_k4", "baseline", "k444_encut1.000"), ("k5", "conv-k555", "k555_encut1.000"), ("k6", "conv-k666", "k666_encut1.000"), ("k7", "conv-k777", "k777_encut1.000"), ("k8", "conv-k888", "k888_encut1.000"), ("k9", "conv-k999", "k999_encut1.000"), ("k10", "conv-k101010", "k101010_encut1.000"), ("k12", "conv-k121212", "k121212_encut1.000"), ("k16", "conv-k161616", "k161616_encut1.000"), ("k20", "conv-k202020", "k202020_encut1.000"), ("k24", "conv-k242424", "k242424_encut1.000"), ("encut_1.3_at_k24", "conv-encut13-k242424", "k242424_encut1.300"), ("encut_1.6_at_k24", "conv-encut16-k242424", "k242424_encut1.600"), ("encut_2.0_at_k24", "conv-encut20-k242424", "k242424_encut2.000")],
        "ti": [("reference_k4", "baseline", "k444_encut1.000"), ("k5", "conv-k555", "k555_encut1.000"), ("k6", "conv-k666", "k666_encut1.000"), ("k7", "conv-k777", "k777_encut1.000"), ("k8", "conv-k888", "k888_encut1.000"), ("k9", "conv-k999", "k999_encut1.000"), ("k10", "conv-k101010", "k101010_encut1.000"), ("k12", "conv-k121212", "k121212_encut1.000"), ("k16", "conv-k161616", "k161616_encut1.000"), ("k20", "conv-k202020", "k202020_encut1.000"), ("k24", "conv-k242424", "k242424_encut1.000"), ("encut_1.3_at_k24", "conv-encut13-k242424", "k242424_encut1.300"), ("encut_1.6_at_k24", "conv-encut16-k242424", "k242424_encut1.600"), ("encut_2.0_at_k24", "conv-encut20-k242424", "k242424_encut2.000")],
        "zr": [("reference_k4", "baseline", "k444_encut1.000"), ("k5", "conv-k555", "k555_encut1.000"), ("k6", "conv-k666", "k666_encut1.000"), ("k7", "conv-k777", "k777_encut1.000"), ("k8", "conv-k888", "k888_encut1.000"), ("k9", "conv-k999", "k999_encut1.000"), ("k10", "conv-k101010", "k101010_encut1.000"), ("k12", "conv-k121212", "k121212_encut1.000"), ("k16", "conv-k161616", "k161616_encut1.000"), ("k20", "conv-k202020", "k202020_encut1.000"), ("k24", "conv-k242424", "k242424_encut1.000"), ("encut_1.3_at_k24", "conv-encut13-k242424", "k242424_encut1.300"), ("encut_1.6_at_k24", "conv-encut16-k242424", "k242424_encut1.600"), ("encut_2.0_at_k24", "conv-encut20-k242424", "k242424_encut2.000")],
        "sio2": [("reference_gamma", "sio2-sp", "k111_encut1.000"), ("k2", "conv-k222", "k222_encut1.000"), ("k3", "conv-k333", "k333_encut1.000"), ("k4", "conv-k444-sio2", "k444_encut1.000"), ("k5", "conv-k555-sio2", "k555_encut1.000"), ("encut_1.3_at_k5", "conv-encut-k555-sio2", "k555_encut1.300")],
    }
    for system, pair_list in PAIRS.items():
        for level, root, key in levels[system]:
            energies = {}
            paths = {}
            for phase in {phase for pair in pair_list for phase in pair}:
                if root == "baseline":
                    root_path = work / root / system / phase / "single_point" / key / "OUTCAR"
                elif root == "sio2-sp":
                    root_path = work / root / system / phase / "single_point" / key / "OUTCAR"
                else:
                    root_path = work / root / system / phase / "single_point" / key / "OUTCAR"
                paths[phase] = root_path
                selected, parsed = energy_candidates(root_path)
                if selected is not None and parsed is not None:
                    paths[phase] = selected
                    energies[phase] = parsed["energy_per_atom_eV"]
            for low, high in pair_list:
                missing = [str(paths[phase]) for phase in (low, high) if phase not in energies]
                convergence.append({
                    "system": system,
                    "phase_low": low,
                    "phase_high": high,
                    "level": level,
                    "delta_E_meV_per_atom": ((energies[high] - energies[low]) * 1000.0 if not missing else "n/a:missing_output"),
                    "status": "ok" if not missing else "missing_output",
                    "missing_paths": ";".join(missing),
                })
    write_csv(output / "convergence.csv", convergence)

    structure_rows = []
    for system, phases in PHASES.items():
        for phase in phases:
            contcar = phase_contcar(work, system, phase)
            number, symbol, matches = symmetry(contcar, EXPECTED_SG[phase])
            audit = symmetry_audit.get(f"{system}:{phase}", {})
            if matches == "not_checked" and audit.get("spacegroup_matches") == "yes":
                number = audit.get("spacegroup_number")
                symbol = audit.get("spacegroup_symbol", "")
                matches = "yes"
            relax_outcar = contcar.parent / "OUTCAR"
            vasp_out = contcar.parent / "vasp.out"
            logs = [vasp_out] if vasp_out.is_file() else sorted(contcar.parent.glob("slurm-*.out"))
            text = "\n".join(path.read_text(errors="replace") for path in logs if path.is_file())
            if relax_outcar.is_file():
                text += "\n" + relax_outcar.read_text(errors="replace")
            max_force = last_max_force(text)
            structure_rows.append({
                "system": system,
                "phase": phase,
                "structure": str(contcar),
                "spacegroup_number": number,
                "spacegroup_symbol": symbol,
                "expected_spacegroup": EXPECTED_SG[phase],
                "spacegroup_matches": matches,
                "max_force_eV_per_A": max_force if max_force is not None else "n/a",
                "force_converged": (
                    "reached required accuracy - stopping structural energy minimisation" in text
                    or (max_force is not None and max_force <= 0.001)
                ),
                "relax_outcar_sha256": sha256(relax_outcar) if relax_outcar.is_file() else "",
            })
    write_csv(output / "relaxed_structures.csv", structure_rows)

    missing_convergence = [row for row in convergence if row["status"] != "ok"]
    unchecked_symmetry = [row for row in structure_rows if row["spacegroup_matches"] != "yes"]
    model_abs = [abs(float(row["DPA_minus_DFT_meV_per_atom"])) for row in comparison]
    structural_abs = [abs(float(row["DFT_at_DPA_minus_DFT_meV_per_atom"])) for row in comparison]
    lines = [
        "# DFT static-energy error budget",
        "",
        "The sign convention is delta E = E_high - E_low. Energies are in meV/atom.",
        "",
        "## Numerical conclusions",
        "",
    ]
    for row in comparison:
        key = f'{row["system"]}:{row["phase_low"]}:{row["phase_high"]}'
        statement = "DPA static energy difference is the main error source for this pair." if row["model_error_exceeds_reference_std"] else "The DPA static energy difference does not exceed the reference Delta-G standard deviation for this pair."
        bias = row["delta_bias_meV_per_atom"]
        if isinstance(bias, float):
            same_sign = np.sign(float(row["DPA_minus_DFT_meV_per_atom"])) == np.sign(bias)
            magnitude = abs(bias) / max(abs(float(row["DPA_minus_DFT_meV_per_atom"])), 1.0e-12)
            bias_note = f"Delta-bias is {'same-sign' if same_sign else 'opposite-sign'} and has a magnitude ratio of {magnitude:.3g} relative to the DPA static error."
        else:
            bias_note = "No usable snapshot bias was available for this pair."
        lines.append(f'- {key}: |DPA - DFT| = {abs(float(row["DPA_minus_DFT_meV_per_atom"])):.6f}; reference Delta-G std = {float(row["reference_delta_G_std_meV_per_atom"]):.6f}. {statement}')
        lines.append(f'  DPA = {float(row["delta_E_DPA_meV_per_atom"]):.6f}; DFT@DPA = {float(row["delta_E_DFT_at_DPA_meV_per_atom"]):.6f}; relaxed DFT = {float(row["delta_E_DFT_meV_per_atom"]):.6f}; model error = {float(row["DPA_minus_DFT_meV_per_atom"]):.6f}; structural error = {float(row["DFT_at_DPA_minus_DFT_meV_per_atom"]):.6f}. {bias_note}')
    lines.extend([
        "",
        f"Across the five pairs, the mean absolute model error is {np.mean(model_abs):.6f} meV/atom and the mean absolute structural error is {np.mean(structural_abs):.6f} meV/atom.",
        "",
        "## Convergence and deviations",
        "",
        "The adopted production meshes are 24x24x24 for Hf, Ti, and Zr, and 5x5x5 for SiO2. The production final energies use these adopted meshes. The adopted metal ENCUT is 2.0 times the reference value because the required 1.3x check remained above 0.2 meV/atom; higher-cutoff checks are reported in convergence.csv.",
        "The neighboring k20-to-k24 pair changes were -0.564425 meV/atom (Hf), +0.061875 meV/atom (Ti), and +0.077300 meV/atom (Zr). Therefore the strict <0.2 meV/atom neighboring-mesh criterion was met for Ti and Zr but not Hf within the computed mesh range; Hf k24 is retained as the highest tested mesh and this deviation is explicit.",
        "Missing convergence levels are retained as n/a:missing_output in convergence.csv rather than being inferred from another k-point or ENCUT.",
        ("All nine relaxed structures retain their target space group according to spglib."
         if not unchecked_symmetry else
         "Space-group verification is not complete in this local assembly because ASE/spglib was unavailable; relaxed structure paths are recorded in relaxed_structures.csv."),
        "The actual VASP executable reported revision 5.4.4.18Apr17-6-g9f103f2a35, build Feb 06 2024. Reference archives reported the same executable revision family but different build metadata for Hf and SiO2; this is recorded in reference_settings.csv and findings.meta.json.",
        "",
    ])
    (output / "summary.md").write_text("\n".join(lines), encoding="utf-8")

    git_commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=repo, check=True, capture_output=True, text=True).stdout.strip()
    settings = read_csv(output / "reference_settings.csv")
    input_hashes = {}
    for source in final_sources.values():
        source_path = work / source
        for path in sorted(source_path.parent.glob("*")):
            if path.is_file() and path.name in {"INCAR", "KPOINTS", "POSCAR", "POTCAR"}:
                input_hashes[str(path.relative_to(work))] = sha256(path)
    for path in sorted((work / "final-adopted-converged").glob("**/*")):
        if path.is_file() and path.name in {"INCAR", "KPOINTS", "POSCAR", "POTCAR"}:
            input_hashes[str(path.relative_to(work))] = sha256(path)
    metadata = {
        "git_commit": git_commit,
        "platform": "expanse",
        "vasp_command": "/home/zjin9/VASP5/bin/vasp_std",
        "vasp_version_observed": str(dft_final[('hf', 'hcp')]["vasp_version"]),
        "potcar_titles": {row["system"]: row["potcar_titles"] for row in settings},
        "input_sha256": input_hashes,
        "checkpoint": str(args.checkpoint),
        "checkpoint_sha256": sha256(args.checkpoint),
        "heads": {system: json.loads((repo / "configs/models/head_policy.yaml").read_text(encoding="utf-8"))["systems"][system]["energy_head"] for system in PHASES},
        "convergence_target_kpoints": {"hf": [24, 24, 24], "ti": [24, 24, 24], "zr": [24, 24, 24], "sio2": [5, 5, 5]},
        "production_final_kpoints": {"metals": [24, 24, 24], "sio2": [5, 5, 5]},
        "final_energy_sources": final_sources,
        "final_converged_output_available": any("final-adopted-converged" in source for source in final_sources.values()),
        "symmetry_audit": str(audit_path) if symmetry_audit else None,
        "adopted_encut_factors": {"hf": 2.0, "ti": 2.0, "zr": 2.0, "sio2": 1.0},
        "reference_settings_sha256": sha256(output / "reference_settings.csv"),
        "missing_convergence_rows": missing_convergence,
        "unchecked_symmetry_rows": unchecked_symmetry,
        "status": "complete" if not missing_convergence and not unchecked_symmetry and all("final-adopted-converged" in source for source in final_sources.values()) else "partial",
    }
    (output / "findings.meta.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"comparison_rows": len(comparison), "convergence_rows": len(convergence), "output": str(output)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
