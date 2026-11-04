"""PROTOCOL_RE1.md: reference-efficiency development experiment. Per entry: nested random subsets of the archived
reference displacement configurations (budgets 10/25/50/100%, two repeats), and for every (cutoff, lambda) in a fixed
grid the methods
  M1 zero-prior A1 total force constants fitted to reference forces (no MLIP),
  M2 MLIP + A1 residual (prior centred on the MLIP everywhere),
  M3 A1 total force constants with a ridge prior centred on the MLIP's own A1 projection (attribution control),
plus M0 (fixed MLIP). All through the same fixed projection and phonopy path (S1 handled at scoring).
Usage: re1.py <out_dir> <entry.yaml.bz2> [...]   (resume-safe; runs on GenSi)
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
TG = np.arange(0, 1501, 50, dtype=float)
BUDGETS = (0.10, 0.25, 0.50, 1.00); REPEATS = (1, 2); LAMS = (1e-6, 1e-3, 1e-1, 1e1); CUTS = ("rc", "rc4")


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
    ph0 = Phonopy(cell, supercell_matrix=scm, primitive_matrix=pm); sc = ph0.supercell; L = sc.cell; pos0 = sc.positions; N = len(sc)
    ssym = list(sc.symbols); types = np.array([tmap.index(s) for s in ssym], dtype=np.int32); disp = d["displacements"]; nd = len(disp)
    dimg = min(np.linalg.norm(np.array([i, j, k]) @ L) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1) if (i, j, k) != (0, 0, 0))
    rcs = {"rc": min(6.0, 0.5 * dimg - 0.05)}; rcs["rc4"] = min(4.0, rcs["rc"])

    def ev(pos):
        f0 = np.asarray(dp.eval(pos[None], L.reshape(1, 9), types)[1]).reshape(-1, 3)
        a = Atoms(ssym, positions=pos, cell=L, pbc=True); a.calc = ms; return {"F0": f0, "MS": a.get_forces()}
    Fd = [np.array(x["forces"]) for x in disp]; Fm = {"F0": [], "MS": []}; R0 = ev(pos0)
    for x in disp:
        p = pos0.copy(); p[x["atom"] - 1] += np.array(x["displacement"]); f = ev(p); Fm["F0"].append(f["F0"] - R0["F0"]); Fm["MS"].append(f["MS"] - R0["MS"])

    def new_ph():
        return Phonopy(cell, supercell_matrix=scm, primitive_matrix=pm)

    def fc_from(forces):
        ph = new_ph()
        ph.dataset = {"natom": N, "first_atoms": [dict(number=x["atom"] - 1, displacement=np.array(x["displacement"]), forces=np.asarray(f)) for x, f in zip(disp, forces)]}
        ph.produce_force_constants(); return ph.force_constants.copy()

    def project(fc):
        ph = new_ph(); ph.force_constants = fc.copy()
        ph.symmetrize_force_constants(level=1, show_drift=False); ph.symmetrize_force_constants_by_space_group(show_drift=False); return ph.force_constants.copy()

    def thermo(fc):
        """fc must already satisfy the projection (projected MLIP part plus symmetric hiPhive part)."""
        ph = new_ph(); ph.force_constants = fc.copy()
        ph.run_mesh(100.0, with_eigenvectors=False); fmin = float(ph.mesh.frequencies.min())
        ph.run_thermal_properties(temperatures=TG); n = len(ph.primitive)
        return np.array(ph.get_thermal_properties_dict()["free_energy"]) * 1000 / 96.48533212 / n, fmin

    partner = {i: j for i, xi in enumerate(disp) for j, xj in enumerate(disp)
               if i != j and xi["atom"] == xj["atom"] and np.allclose(np.array(xi["displacement"]), -np.array(xj["displacement"]))}

    def labels(forces, subset):
        s = set(subset); return [(forces[i] - forces[partner[i]]) / 2 if (i in partner and partner[i] in s) else forces[i] for i in subset]

    res = dict(T=TG, natom_sc=N, n_disp=nd, dimg=dimg, rc=rcs["rc"], rc4=rcs["rc4"], n_pm_pairs=len(partner) // 2)
    t_th = time.time(); Kraw = {M: fc_from(Fm[M]) for M in ("F0", "MS")}; K = {M: project(Kraw[M]) for M in ("F0", "MS")}   # projection is linear and idempotent
    for M in ("F0", "MS"): res[f"F_M0_{M}"], res[f"fmin_M0_{M}"] = thermo(K[M])
    res["t_thermo"] = (time.time() - t_th) / 2
    prim = Atoms(sym, scaled_positions=frac, cell=lat, pbc=True); scA = Atoms(ssym, positions=pos0, cell=L, pbc=True)
    # nested random subsets per repeat (seeded by mp-id and repeat); 100% is the full set
    subsets = {}
    for r in REPEATS:
        perm = np.random.default_rng(int(hashlib.sha256(f"{m_id}|{r}".encode()).hexdigest()[:8], 16)).permutation(nd)
        for b in BUDGETS:
            k = nd if b == 1.0 else max(1, int(round(b * nd))); subsets[(b, r)] = tuple(sorted(perm[:k].tolist()))
    for (b, r), s in subsets.items(): res[f"nconf_{b}_{r}"] = len(s)
    for cut in CUTS:
        try:
            cs = ClusterSpace(prim, [rcs[cut]])
        except ValueError:
            cs = None
        res[f"npar_{cut}"] = 0 if cs is None else cs.n_dofs
        if cs is None: continue
        rows = []   # sensing rows per configuration (depend only on the displacement)
        for x in disp:
            at = scA.copy(); u = np.zeros((N, 3)); u[x["atom"] - 1] = np.array(x["displacement"]); at.new_array("displacements", u); at.new_array("forces", np.zeros((N, 3)))
            sct = StructureContainer(cs); sct.add_structure(at); rows.append(sct.get_fit_data()[0])
        Afull = np.vstack(rows); npar = Afull.shape[1]; scale = np.trace(Afull.T @ Afull) / npar
        theta0 = {}
        for M in ("F0", "MS"):   # MLIP's own A1 projection: fit to MLIP response forces on all configurations (free)
            yM = np.concatenate([f.ravel() for f in labels(Fm[M], list(range(nd)))]); AtA = Afull.T @ Afull
            theta0[M] = np.linalg.solve(AtA + 1e-8 * scale * np.eye(npar), Afull.T @ yM)
        fcp_cache = {}

        Au = unit_map(cs, scA, N)

        def Phi(theta):   # identical to ForceConstantPotential(cs, theta).get_force_constants(scA) (validated to 1e-15)
            if Au is None:
                return ForceConstantPotential(cs, theta).get_force_constants(scA).get_fc_array(order=2)
            return -(Au @ theta).reshape(N, 3, N, 3).transpose(2, 0, 3, 1)
        done = {}
        for (b, r), s in subsets.items():
            if s in done:   # identical subset (nested budgets collapse for small n_disp): reuse
                for key, val in done[s].items(): res[key.format(b=b, r=r)] = val
                continue
            A = np.vstack([rows[i] for i in s]); AtA = A.T @ A
            yref = np.concatenate([f.ravel() for f in labels(Fd, list(s))]); store = {}
            for lam in LAMS:
                Lm = lam * scale * np.eye(npar)
                th1 = np.linalg.solve(AtA + Lm, A.T @ yref)
                F1, f1 = thermo(Phi(th1)); store[f"F_M1_{cut}_{lam:g}_{{b}}_{{r}}"] = F1; store[f"fmin_M1_{cut}_{lam:g}_{{b}}_{{r}}"] = f1
                for M in ("F0", "MS"):
                    yres = np.concatenate([(fd - fm).ravel() for fd, fm in zip(labels(Fd, list(s)), labels(Fm[M], list(s)))])
                    th2 = np.linalg.solve(AtA + Lm, A.T @ yres); P2 = Phi(th2); F2, f2 = thermo(K[M] + P2)
                    if "proj_equiv_check" not in res and lam == LAMS[0]:   # proj(K + Phi) == proj(K) + Phi
                        Fchk, _ = thermo(project(Kraw[M] + P2)); res["proj_equiv_check"] = float(np.abs(Fchk - F2).max())
                    th3 = np.linalg.solve(AtA + Lm, A.T @ yref + Lm @ theta0[M]); F3, f3 = thermo(Phi(th3))
                    store[f"F_M2_{M}_{cut}_{lam:g}_{{b}}_{{r}}"] = F2; store[f"fmin_M2_{M}_{cut}_{lam:g}_{{b}}_{{r}}"] = f2
                    store[f"F_M3_{M}_{cut}_{lam:g}_{{b}}_{{r}}"] = F3; store[f"fmin_M3_{M}_{cut}_{lam:g}_{{b}}_{{r}}"] = f3
            done[s] = store
            for key, val in store.items(): res[key.format(b=b, r=r)] = val
    np.savez(out, **res)
    print(f"{m_id}: N={N} nd={nd} npar={res['npar_rc']}/{res['npar_rc4']} subsets={len(set(subsets.values()))} t_thermo={res['t_thermo']:.1f}s total {time.time() - t0:.0f} s", flush=True)


for src in entries:
    out = out_dir / (Path(src).name.split(".")[0] + ".npz")
    if out.exists(): continue
    try:
        run(src, out)
    except Exception as exc:
        (out_dir / (out.stem + ".FAILED.txt")).write_text(repr(exc) + "\n"); print(out.stem, "FAILED", repr(exc), flush=True)
