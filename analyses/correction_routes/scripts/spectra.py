"""PROTOCOL.md data step: band-matched mode spectra for DFT, F0 and MatterSim plus structure/composition features (GenSi).

Same phonopy path as mlip_vs_dft_phonons.py (stored displacements, produce_force_constants, mesh 100). Per model,
modes are matched to DFT by sorted order at each q-point and binned by model ln(nu): 480 log bins over 1e-4-60 THz
holding every model mode with nu > 0 (phonopy's thermal cutoff is 0 THz; DEVIATIONS.md D1). The ln-ratio target uses
only matched modes with both frequencies > 0.05 THz (separate weights `_wd`).
Usage: spectra.py <out_dir> <entry.yaml.bz2> [...]   (resume-safe; one npz per entry)
"""
import bz2, sys, time
from pathlib import Path
import numpy as np, torch, yaml, spglib
torch.serialization.add_safe_globals([slice])
from deepmd.infer import DeepPot  # noqa: E402
from mattersim.forcefield import MatterSimCalculator  # noqa: E402
from ase import Atoms  # noqa: E402
from phonopy import Phonopy  # noqa: E402
from phonopy.structure.atoms import PhonopyAtoms  # noqa: E402
from pymatgen.core import Element  # noqa: E402

out_dir = Path(sys.argv[1]); entries = sys.argv[2:]
dp = DeepPot("/GenSIvePFS/users/zirenj/fes_experiment_pbe_d3bj/models/DPA-3.1-3M.pt", head="MP_traj_v024_alldata_mixu"); tmap = dp.get_type_map()
ms = MatterSimCalculator(load_path="/root/.local/mattersim/pretrained_models/mattersim-v1.0.0-5M.pth", device="cuda")
EDGES = np.linspace(np.log(1e-4), np.log(60.0), 481); NB_ = 480; NUMIN = 0.05


def crystal_system(n):
    return int(np.searchsorted([2, 15, 74, 142, 167, 194, 230], n, side="left"))


def run(src, out):
    t0 = time.time(); d = yaml.load(bz2.open(src), Loader=yaml.CSafeLoader)
    uc = d["unit_cell"]; pts = uc["points"]; sym = [p["symbol"] for p in pts]; lat = np.array(uc["lattice"]); frac = np.array([p["coordinates"] for p in pts])
    masses = [p["mass"] for p in pts]
    cell = PhonopyAtoms(symbols=sym, cell=lat, scaled_positions=frac, masses=masses)
    pm = np.array(d.get("primitive_matrix", np.eye(3))); scm = np.array(d["supercell_matrix"])
    probe = Phonopy(cell, supercell_matrix=scm, primitive_matrix=pm); sc = probe.supercell; L = sc.cell
    pos0 = sc.positions; types = np.array([tmap.index(s) for s in sc.symbols], dtype=np.int32)
    forces = {"DFT": [np.array(x["forces"]) for x in d["displacements"]], "F0": [], "MS": []}
    for x in d["displacements"]:
        p = pos0.copy(); p[x["atom"] - 1] += np.array(x["displacement"])
        forces["F0"].append(np.asarray(dp.eval(p[None], L.reshape(1, 9), types)[1]).reshape(-1, 3))
        a = Atoms(list(sc.symbols), positions=p, cell=L, pbc=True); a.calc = ms; forces["MS"].append(a.get_forces())
    nu, Fphon = {}, {}
    for side, fl in forces.items():
        ph = Phonopy(cell, supercell_matrix=scm, primitive_matrix=pm)
        ph.dataset = {"natom": len(sc), "first_atoms": [dict(number=x["atom"] - 1, displacement=np.array(x["displacement"]), forces=f) for x, f in zip(d["displacements"], fl)]}
        ph.produce_force_constants(); ph.run_mesh(100.0, with_eigenvectors=False)
        nu[side] = np.sort(ph.mesh.frequencies, axis=1); wq = ph.mesh.weights / ph.mesh.weights.sum(); n_prim = len(ph.primitive)
        ph.run_thermal_properties(temperatures=np.arange(0, 1501, 50, dtype=float))
        Fphon[side] = np.array(ph.get_thermal_properties_dict()["free_energy"]) * 1000 / 96.48533212 / n_prim   # same-run phonopy F (meV/atom)
    W = np.broadcast_to(wq[:, None], nu["DFT"].shape) / n_prim   # per-atom mode weight
    res = {f"{k}_Fphon": v for k, v in Fphon.items()}
    for side in ("F0", "MS"):
        m, r = nu[side], nu["DFT"]; pos = m > 0; keep = m > NUMIN; tgt = keep & (r > NUMIN)
        lm = np.log(np.where(pos, m, 1.0)); lr = np.log(np.where(tgt, r, 1.0)) - lm
        b = np.clip(np.digitize(lm, EDGES) - 1, 0, NB_ - 1)
        w_b = np.bincount(b[pos], W[pos], NB_); s_m = np.bincount(b[pos], (W * lm)[pos], NB_)
        wd_b = np.bincount(b[tgt], W[tgt], NB_); s_d = np.bincount(b[tgt], (W * lr)[tgt], NB_)
        nz, nzd = w_b > 0, wd_b > 0
        res[f"{side}_w"] = w_b; res[f"{side}_m"] = np.where(nz, s_m / np.where(nz, w_b, 1), 0.0)
        res[f"{side}_wd"] = wd_b; res[f"{side}_d"] = np.where(nzd, s_d / np.where(nzd, wd_b, 1), 0.0)
        res[f"{side}_excluded_w"] = float(W[~keep].sum()); res[f"{side}_nonpositive_w"] = float(W[~pos].sum())
        lv = lm[keep]; wv = W[keep]
        res[f"{side}_spec"] = np.array([np.average(lv, weights=wv), np.sqrt(np.average((lv - np.average(lv, weights=wv)) ** 2, weights=wv)),
                                        W[m < 0.5].sum() / W.sum(), W[m < 1.0].sum() / W.sum(), np.quantile(lv, 0.99), float(m.min()), W[~keep].sum() / W.sum()])
    nums = [sorted(set(sym)).index(s) + 1 for s in sym]
    sg = spglib.get_symmetry_dataset((lat, frac, nums), symprec=1e-3); sgn = int(sg.number if hasattr(sg, "number") else sg["number"])
    Z = np.array([Element(s).Z for s in sym], float); M = np.array(masses, float); X = np.array([Element(s).X for s in sym], float)
    vol = abs(np.linalg.det(lat))
    res["struct"] = np.array([vol / len(sym), M.sum() / vol, n_prim, sgn, crystal_system(sgn)])
    res["comp"] = np.array([Z.mean(), Z.std(), M.mean(), M.std(), np.nanmean(X), np.nanstd(X), len(set(sym))])
    np.savez(out, **res)
    print(f"{out.stem}: {time.time() - t0:.0f} s", flush=True)


for src in entries:
    out = out_dir / (Path(src).name.split(".")[0] + ".npz")
    if out.exists():
        continue
    try:
        run(src, out)
    except Exception as exc:
        (out_dir / (out.stem + ".FAILED.txt")).write_text(repr(exc) + "\n"); print(out.stem, "FAILED", repr(exc), flush=True)
