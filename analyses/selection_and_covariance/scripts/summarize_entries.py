"""Per-entry reference summary for the Alexandria PBE phonon files (run in the deepmd-kit WSL env: phonopy, spglib).

For each file: formula, space group, static energy and volume per atom, and the stored thermal free energy
(kJ/mol per phonopy primitive cell -> meV/atom) at the five stored temperatures. It also rebuilds the phonon from the stored
force constants (stored primitive matrix, 30x30x30-equivalent mesh by density) and checks the stored free energy is
reproduced, and records the lowest mesh frequency (dynamical-stability re-check).
Usage: summarize_entries.py <raw_dir> <out.csv> [mp-id ...]
Per-entry results are cached in <out.csv>.cache/ so reruns only process new files.
"""
import bz2, csv, json, sys
from pathlib import Path
import numpy as np, yaml
LOADER = getattr(yaml, "CSafeLoader", yaml.SafeLoader)  # libyaml when available; same safe semantics
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

KJMOL_TO_MEV = 1000 / 96.48533212        # 1 kJ/mol = 10.3643 meV
raw, out = Path(sys.argv[1]), Path(sys.argv[2])
ids = sys.argv[3:] or sorted(p.name.split(".")[0] for p in raw.glob("mp-*.yaml.bz2"))


def formula(symbols):
    from collections import Counter
    from math import gcd
    from functools import reduce
    c = Counter(symbols); g = reduce(gcd, c.values())
    return "".join(f"{el}{'' if n // g == 1 else n // g}" for el, n in sorted(c.items()))


cache = Path(str(out) + ".cache"); cache.mkdir(exist_ok=True)
rows = []
for mid in ids:
    cf = cache / f"{mid}.json"
    if cf.exists():
        rows.append(json.loads(cf.read_text())); continue
    d = yaml.load(bz2.open(raw / f"{mid}.yaml.bz2"), Loader=LOADER)
    uc = d["unit_cell"]; pts = uc["points"]; sym = [p["symbol"] for p in pts]; n = len(pts)
    cell = PhonopyAtoms(symbols=sym, cell=np.array(uc["lattice"]), scaled_positions=np.array([p["coordinates"] for p in pts]),
                        masses=[p["mass"] for p in pts])
    ph = Phonopy(cell, supercell_matrix=np.array(d["supercell_matrix"]), primitive_matrix=np.array(d.get("primitive_matrix", np.eye(3))))
    fc = np.array(d["force_constants"]["elements"]).reshape(d["force_constants"]["shape"] + [3, 3])
    ph.force_constants = fc
    ph.run_mesh(100.0, with_eigenvectors=False)
    fmin = float(ph.mesh.frequencies.min())
    T = np.array(d["temperatures"], dtype=float)
    ph.run_thermal_properties(temperatures=T)
    # stored free_e, like phonopy's thermal properties, is per primitive cell (checked: BaTe with a primitive matrix)
    n_prim = len(ph.primitive)
    f_re = np.array(ph.get_thermal_properties_dict()["free_energy"])  # dict API works in phonopy 2.x and 4.x
    f_st = np.array(d["free_e"], dtype=float)
    rows.append(dict(mp_id=mid, formula=formula(sym), spacegroup=d["space_group"]["number"], natoms_unit=n,
                     natoms_supercell=len(d["supercell"]["points"]), n_displacements=len(d["displacements"]),
                     energy_eV_atom=d["energy"] / n, volume_A3_atom=d["volume"] / n, min_mesh_freq_THz=fmin,
                     natoms_primitive=n_prim, **{f"Fvib_{int(t)}K_meV_atom": f * KJMOL_TO_MEV / n_prim for t, f in zip(T, f_st)},
                     Fvib_rebuild_max_abs_diff_meV_atom=float(abs(f_re - f_st).max() * KJMOL_TO_MEV / n_prim),
                     phonopy_version_file=d["phonopy"]["version"]))
    cf.write_text(json.dumps(rows[-1]))
    print(mid, rows[-1]["formula"], rows[-1]["spacegroup"], f"fmin {fmin:.3f}", f"rebuild diff {rows[-1]['Fvib_rebuild_max_abs_diff_meV_atom']:.3g} meV/atom", flush=True)
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
