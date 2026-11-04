"""PROTOCOL_RE1_MECH.md: unobserved-direction thermodynamic error (D3), identifiability and the force-error control,
on the same subsets as RE1 (same seed rule), A1 basis at r_c. Diagnostic only; never an input to any RE1 method.
Usage: re1_mech.py <out_dir> <entry.yaml.bz2> [...]   (resume-safe; runs on GenSi)
"""
import bz2, hashlib, logging, sys, time
from pathlib import Path
import numpy as np, torch, yaml
torch.serialization.add_safe_globals([slice])
from deepmd.infer import DeepPot  # noqa: E402
from mattersim.forcefield import MatterSimCalculator  # noqa: E402
from ase import Atoms  # noqa: E402
from phonopy import Phonopy  # noqa: E402
from phonopy.structure.atoms import PhonopyAtoms  # noqa: E402
from hiphive import ClusterSpace, StructureContainer, ForceConstantPotential  # noqa: E402
logging.getLogger("hiphive").setLevel(logging.ERROR)

out_dir = Path(sys.argv[1]); entries = sys.argv[2:]
dp = DeepPot("/GenSIvePFS/users/zirenj/fes_experiment_pbe_d3bj/models/DPA-3.1-3M.pt", head="MP_traj_v024_alldata_mixu"); tmap = dp.get_type_map()
ms = MatterSimCalculator(load_path="/root/.local/mattersim/pretrained_models/mattersim-v1.0.0-5M.pth", device="cuda")
TG = np.arange(0, 1501, 50, dtype=float); BUDGETS = (0.10, 0.25, 0.50, 1.00); REPEATS = (1, 2); NPMAX = 2000


def unit_map(cs, scA, N, max_elems=1e8):
    """Linear map theta -> supercell force constants from 3N unit displacements (implementation speed-up only)."""
    if 9 * N * N * cs.n_dofs > max_elems: return None
    sct = StructureContainer(cs)
    for a in range(N):
        for x in range(3):
            at = scA.copy(); u = np.zeros((N, 3)); u[a, x] = 1.0; at.new_array("displacements", u); at.new_array("forces", np.zeros((N, 3))); sct.add_structure(at)
    return sct.get_fit_data()[0]   # rows ordered (a, x, i, y)


def run(src, out):
    t0 = time.time(); m_id = Path(src).name.split(".")[0]; d = yaml.load(bz2.open(src), Loader=yaml.CSafeLoader)
    uc = d["unit_cell"]; pts = uc["points"]; sym = [p["symbol"] for p in pts]; lat = np.array(uc["lattice"]); frac = np.array([p["coordinates"] for p in pts])
    cell = PhonopyAtoms(symbols=sym, cell=lat, scaled_positions=frac, masses=[p["mass"] for p in pts])
    pm = np.array(d.get("primitive_matrix", np.eye(3))); scm = np.array(d["supercell_matrix"])
    sc = Phonopy(cell, supercell_matrix=scm, primitive_matrix=pm).supercell; L = sc.cell; pos0 = sc.positions; N = len(sc)
    ssym = list(sc.symbols); types = np.array([tmap.index(s) for s in ssym], dtype=np.int32); disp = d["displacements"]; nd = len(disp)
    dimg = min(np.linalg.norm(np.array([i, j, k]) @ L) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1) if (i, j, k) != (0, 0, 0))
    rc = min(6.0, 0.5 * dimg - 0.05); res = dict(T=TG, natom_sc=N, n_disp=nd, rc=rc)

    def ev(pos):
        f0 = np.asarray(dp.eval(pos[None], L.reshape(1, 9), types)[1]).reshape(-1, 3)
        a = Atoms(ssym, positions=pos, cell=L, pbc=True); a.calc = ms; return {"F0": f0, "MS": a.get_forces()}
    Fd = [np.array(x["forces"]) for x in disp]; Fm = {"F0": [], "MS": []}; R0 = ev(pos0)
    for x in disp:
        p = pos0.copy(); p[x["atom"] - 1] += np.array(x["displacement"]); f = ev(p); Fm["F0"].append(f["F0"] - R0["F0"]); Fm["MS"].append(f["MS"] - R0["MS"])
    for M in ("F0", "MS"): res[f"forceRMS_{M}"] = float(np.sqrt(np.mean(np.concatenate([(a - b).ravel() for a, b in zip(Fm[M], Fd)]) ** 2)))
    partner = {i: j for i, xi in enumerate(disp) for j, xj in enumerate(disp)
               if i != j and xi["atom"] == xj["atom"] and np.allclose(np.array(xi["displacement"]), -np.array(xj["displacement"]))}

    def labels(forces, subset):
        s = set(subset); return [(forces[i] - forces[partner[i]]) / 2 if (i in partner and partner[i] in s) else forces[i] for i in subset]

    def thermo(fc):
        ph = Phonopy(cell, supercell_matrix=scm, primitive_matrix=pm); ph.force_constants = fc.copy()
        ph.run_mesh(100.0, with_eigenvectors=False); fmin = float(ph.mesh.frequencies.min())
        ph.run_thermal_properties(temperatures=TG); n = len(ph.primitive)
        return np.array(ph.get_thermal_properties_dict()["free_energy"]) * 1000 / 96.48533212 / n, fmin

    subsets = {}
    for r in REPEATS:   # identical seed rule to re1.py
        perm = np.random.default_rng(int(hashlib.sha256(f"{m_id}|{r}".encode()).hexdigest()[:8], 16)).permutation(nd)
        for b in BUDGETS:
            k = nd if b == 1.0 else max(1, int(round(b * nd))); subsets[(b, r)] = tuple(sorted(perm[:k].tolist()))
    try:
        cs = ClusterSpace(Atoms(sym, scaled_positions=frac, cell=lat, pbc=True), [rc])
    except ValueError:
        cs = None
    res["npar"] = 0 if cs is None else cs.n_dofs
    if cs is None or cs.n_dofs > NPMAX:
        res["applicable"] = False; np.savez(out, **res); print(f"{m_id}: not applicable (npar {res['npar']})", flush=True); return
    res["applicable"] = True; scA = Atoms(ssym, positions=pos0, cell=L, pbc=True); rows = []
    for x in disp:
        at = scA.copy(); u = np.zeros((N, 3)); u[x["atom"] - 1] = np.array(x["displacement"]); at.new_array("displacements", u); at.new_array("forces", np.zeros((N, 3)))
        sct = StructureContainer(cs); sct.add_structure(at); rows.append(sct.get_fit_data()[0])
    Afull = np.vstack(rows); npar = Afull.shape[1]; AtA = Afull.T @ Afull; reg = 1e-8 * np.trace(AtA) / npar * np.eye(npar)
    full = list(range(nd)); yref = np.concatenate([f.ravel() for f in labels(Fd, full)])
    th_ref = np.linalg.solve(AtA + reg, Afull.T @ yref)
    th0 = {M: np.linalg.solve(AtA + reg, Afull.T @ np.concatenate([f.ravel() for f in labels(Fm[M], full)])) for M in ("F0", "MS")}
    Au = unit_map(cs, scA, N)
    Phi = lambda th: (ForceConstantPotential(cs, th).get_force_constants(scA).get_fc_array(order=2) if Au is None
                      else -(Au @ th).reshape(N, 3, N, 3).transpose(2, 0, 3, 1))
    res["F_ref_A1"], res["fmin_ref_A1"] = thermo(Phi(th_ref))
    for M in ("F0", "MS"): res[f"theta_err_norm_{M}"] = float(np.linalg.norm(th0[M] - th_ref) / max(np.linalg.norm(th_ref), 1e-30))
    done = {}
    for (b, r), s in subsets.items():
        if s not in done:
            A = np.vstack([rows[i] for i in s]); U, sv, Vt = np.linalg.svd(A, full_matrices=False); keep = sv > 1e-8 * sv.max(); V = Vt[keep].T
            P = lambda v: v - V @ (V.T @ v)   # (I - P_B) v
            st = dict(rank_frac=float(keep.sum() / npar), cond=float(sv.max() / sv[keep].min()))
            st["F_D3_zero"], st["fmin_D3_zero"] = thermo(Phi(th_ref + P(-th_ref)))
            for M in ("F0", "MS"):
                st[f"F_D3_{M}"], st[f"fmin_D3_{M}"] = thermo(Phi(th_ref + P(th0[M] - th_ref)))
            done[s] = st
        for kk, vv in done[s].items(): res[f"{kk}_{b}_{r}"] = vv
    np.savez(out, **res)
    print(f"{m_id}: N={N} npar={npar} subsets={len(done)} {time.time() - t0:.0f} s", flush=True)


for src in entries:
    out = out_dir / (Path(src).name.split(".")[0] + ".npz")
    if out.exists(): continue
    try:
        run(src, out)
    except Exception as exc:
        (out_dir / (out.stem + ".FAILED.txt")).write_text(repr(exc) + "\n"); print(out.stem, "FAILED", repr(exc), flush=True)
