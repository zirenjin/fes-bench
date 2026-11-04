"""Exploratory (post-hoc): frequency-band decomposition of the F0 thermal free-energy error, per entry.

Same phonopy path as mlip_vs_dft_phonons.py (stored displaced supercells; DFT forces from file, MLIP forces from the
model; mesh density 100). Saves, per entry, the mesh frequencies (THz) of both sides and the q-point weights, so that
mode contributions can be matched band-by-band (sorted per q-point) in the analysis.
Usage: freqband_decomp.py <ckpt> <head> <out_dir> <entry.yaml.bz2> [more ...]   (model loaded once; resume-safe)
"""
import bz2, sys, traceback
from pathlib import Path
import numpy as np, torch, yaml
LOADER = getattr(yaml, "CSafeLoader", yaml.SafeLoader)
torch.serialization.add_safe_globals([slice])
from deepmd.infer import DeepPot  # noqa: E402
from phonopy import Phonopy  # noqa: E402
from phonopy.structure.atoms import PhonopyAtoms  # noqa: E402

ckpt, head, out_dir = sys.argv[1:4]; out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
dp = DeepPot(ckpt, head=head); tmap = dp.get_type_map()


def run(src):
    d = yaml.load(bz2.open(src), Loader=LOADER)
    uc = d["unit_cell"]; pts = uc["points"]
    cell = PhonopyAtoms(symbols=[p["symbol"] for p in pts], cell=np.array(uc["lattice"]), scaled_positions=np.array([p["coordinates"] for p in pts]),
                        masses=[p["mass"] for p in pts])
    pm = np.array(d.get("primitive_matrix", np.eye(3)))

    def mesh(forces_list):
        ph = Phonopy(cell, supercell_matrix=np.array(d["supercell_matrix"]), primitive_matrix=pm)
        ph.dataset = {"natom": len(ph.supercell), "first_atoms": [dict(number=x["atom"] - 1, displacement=np.array(x["displacement"]), forces=np.array(f))
                                                                for x, f in zip(d["displacements"], forces_list)]}
        ph.produce_force_constants(); ph.run_mesh(100.0, with_eigenvectors=False)
        return ph, np.array(ph.mesh.frequencies), np.array(ph.mesh.weights)

    ph, wd, q = mesh([x["forces"] for x in d["displacements"]])
    sc = ph.supercell; types = np.array([tmap.index(s) for s in sc.symbols])
    fm = []
    for x in d["displacements"]:
        pos = sc.positions.copy(); pos[x["atom"] - 1] += np.array(x["displacement"])
        fm.append(np.asarray(dp.eval(pos[None], sc.cell.reshape(1, 9), types)[1]).reshape(-1, 3))
    _, wm, q2 = mesh(fm)
    assert np.array_equal(q, q2)
    return dict(freq_dft=wd, freq_mlip=wm, weights=q, n_prim=len(ph.primitive))


for src in sys.argv[4:]:
    mid = Path(src).name.split(".")[0]; out = out_dir / f"{mid}.npz"
    if out.exists():
        continue
    try:
        np.savez_compressed(out, **run(src)); print(mid, "ok", flush=True)
    except Exception:
        (out_dir / f"{mid}.FAILED.txt").write_text(traceback.format_exc()); print(mid, "FAILED", flush=True)
