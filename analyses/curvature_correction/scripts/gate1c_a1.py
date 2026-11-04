"""PROTOCOL_1c_A1.md: symmetry-built full-tensor second-order residual (hiPhive cluster space) against R3 and the
truncated oracle, all with the same labels (D1), S1 setting and fixed projection.
Usage: gate1c_a1.py <out_dir> <entry.yaml.bz2> [...]   (resume-safe; runs on GenSi)
"""
import bz2, logging, sys, time
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
TG = np.arange(0, 1501, 50, dtype=float); SH = (1.15, 1.5)


def run(src, out):
    t0 = time.time(); d = yaml.load(bz2.open(src), Loader=yaml.CSafeLoader)
    uc = d["unit_cell"]; pts = uc["points"]; sym = [p["symbol"] for p in pts]; lat = np.array(uc["lattice"]); frac = np.array([p["coordinates"] for p in pts])
    cell = PhonopyAtoms(symbols=sym, cell=lat, scaled_positions=frac, masses=[p["mass"] for p in pts])
    pm = np.array(d.get("primitive_matrix", np.eye(3))); scm = np.array(d["supercell_matrix"])
    ph0 = Phonopy(cell, supercell_matrix=scm, primitive_matrix=pm); sc = ph0.supercell; L = sc.cell; pos0 = sc.positions; N = len(sc)
    ssym = list(sc.symbols); types = np.array([tmap.index(s) for s in ssym], dtype=np.int32); disp = d["displacements"]
    dimg = min(np.linalg.norm(np.array([i, j, k]) @ L) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1) if (i, j, k) != (0, 0, 0))
    rc = min(6.0, 0.5 * dimg - 0.05)

    def ev(pos):
        f0 = np.asarray(dp.eval(pos[None], L.reshape(1, 9), types)[1]).reshape(-1, 3)
        a = Atoms(ssym, positions=pos, cell=L, pbc=True); a.calc = ms; return {"F0": f0, "MS": a.get_forces()}
    F = {"DFT": [np.array(x["forces"]) for x in disp], "F0": [], "MS": []}; R0 = ev(pos0)
    for x in disp:
        p = pos0.copy(); p[x["atom"] - 1] += np.array(x["displacement"]); f = ev(p); F["F0"].append(f["F0"]); F["MS"].append(f["MS"])

    def new_ph():
        return Phonopy(cell, supercell_matrix=scm, primitive_matrix=pm)

    def fc_from(forces):
        ph = new_ph()
        ph.dataset = {"natom": N, "first_atoms": [dict(number=x["atom"] - 1, displacement=np.array(x["displacement"]), forces=np.asarray(f)) for x, f in zip(disp, forces)]}
        ph.produce_force_constants(); return ph.force_constants.copy()

    def thermo(fc, project=True):
        ph = new_ph(); ph.force_constants = fc.copy()
        if project:
            ph.symmetrize_force_constants(level=1, show_drift=False); ph.symmetrize_force_constants_by_space_group(show_drift=False)
        ph.run_mesh(100.0, with_eigenvectors=False); fr = ph.mesh.frequencies; k = np.unravel_index(np.argmin(fr), fr.shape)
        q = ph.mesh.qpoints[k[0]]; q = q - np.round(q); fmin = float(fr[k])
        ph.run_thermal_properties(temperatures=TG); n = len(ph.primitive)
        return np.array(ph.get_thermal_properties_dict()["free_energy"]) * 1000 / 96.48533212 / n, fmin, float(np.linalg.norm(q))

    # pair geometry (minimum image)
    fr_ = sc.scaled_positions; dfr = fr_[None, :, :] - fr_[:, None, :]; dfr -= np.round(dfr); rv = dfr @ L; dist = np.linalg.norm(rv, axis=-1); np.fill_diagonal(dist, np.inf)
    keys = {}; cls = -np.ones((N, N), int); dmin = {}
    for i in range(N):
        for j in range(N):
            if dist[i, j] < 6.0:
                k = tuple(sorted((ssym[i], ssym[j]))); dmin[k] = min(dmin.get(k, np.inf), dist[i, j])
    for i in range(N):
        for j in range(N):
            if dist[i, j] < 6.0:
                k = tuple(sorted((ssym[i], ssym[j]))); r = dist[i, j] / dmin[k]; s = 0 if r < SH[0] else (1 if r < SH[1] else 2); cls[i, j] = keys.setdefault((k, s), len(keys))
    rhat = rv / np.where(np.isfinite(dist), dist, 1.0)[..., None]; RR = np.einsum("ija,ijb->ijab", rhat, rhat); nc = len(keys)
    # hiPhive cluster space on the unit cell, supercell atoms in phonopy order
    prim = Atoms(sym, scaled_positions=frac, cell=lat, pbc=True); scA = Atoms(ssym, positions=pos0, cell=L, pbc=True)
    try:
        cs = ClusterSpace(prim, [rc])
    except ValueError:   # DEVIATIONS_1c D1: no second-order dofs within rc (short supercell) -> A1 inapplicable
        cs = None
    res = dict(T=TG, natom_sc=N, n_disp=len(disp), rc=rc, dimg=dimg, n_classes_R3=nc, npar_A1=0 if cs is None else cs.n_dofs)
    K_dft = fc_from(F["DFT"])
    for M in ("F0", "MS"):
        K = fc_from([f - R0[M] for f in F[M]])
        lab = [np.asarray(fd) - (np.asarray(fm) - R0[M]) for fd, fm in zip(F["DFT"], F[M])]
        part = {i: j for i, xi in enumerate(disp) for j, xj in enumerate(disp)
                if i != j and xi["atom"] == xj["atom"] and np.allclose(np.array(xi["displacement"]), -np.array(xj["displacement"]))}
        lab = [(lab[i] - lab[part[i]]) / 2 if i in part else lab[i] for i in range(len(lab))]; y = np.concatenate(lab).ravel()
        # ---- R3 (as gate 1b) ----
        def design(kind):
            cols = []
            for c in range(nc):
                col = []
                for x in disp:
                    a = x["atom"] - 1; u = np.array(x["displacement"]); v = np.zeros((N, 3)); m = cls[:, a] == c
                    blk = K[:, a] if kind == "scale" else -RR[:, a]
                    v[m] = np.einsum("iab,b->ia", blk[m], u); ma = cls[a, :] == c
                    v[a] = -np.einsum("jab,b->a", (K[a] if kind == "scale" else -RR[a])[ma], u)
                    col.append(v.ravel())
                cols.append(-np.concatenate(col))
            return np.array(cols).T
        A = np.hstack([design("scale"), design("radial")]); npr = A.shape[1]; AtA = A.T @ A
        th = np.linalg.solve(AtA + 1e-8 * np.trace(AtA) / npr * np.eye(npr), A.T @ y)
        off = cls >= 0; ts = np.where(off, th[:nc][np.clip(cls, 0, None)], 0.0); tr = np.where(off, th[nc:][np.clip(cls, 0, None)], 0.0)
        dK3 = ts[..., None, None] * K * off[..., None, None] - tr[..., None, None] * RR * off[..., None, None]
        for i in range(N): dK3[i, i] = -dK3[i].sum(axis=0) + dK3[i, i]
        # ---- A1 (hiPhive full tensor) ----
        scont = StructureContainer(cs) if cs is not None else None
        for x, f in (zip(disp, lab) if cs is not None else []):
            at = scA.copy(); u = np.zeros((N, 3)); u[x["atom"] - 1] = np.array(x["displacement"])
            at.new_array("displacements", u); at.new_array("forces", np.asarray(f)); scont.add_structure(at)
        if cs is not None:
            Ah, yh = scont.get_fit_data(); npa = Ah.shape[1]
        else:
            yh = y; npa = np.inf
        res[f"A1_labels_{M}"] = len(yh)
        if npa > 0.5 * len(yh):
            res[f"A1_applicable_{M}"] = False; res[f"F_{M}_A1_proj"] = np.full(len(TG), np.nan); res[f"fmin_{M}_A1_proj"] = np.nan; res[f"qmin_{M}_A1_proj"] = np.nan
            res[f"F_{M}_A1"] = np.full(len(TG), np.nan); res[f"fmin_{M}_A1"] = np.nan; res[f"qmin_{M}_A1"] = np.nan; res[f"r2_{M}_A1"] = np.nan
        else:
            AtAh = Ah.T @ Ah; thh = np.linalg.solve(AtAh + 1e-8 * np.trace(AtAh) / npa * np.eye(npa), Ah.T @ yh)
            res[f"A1_applicable_{M}"] = True; res[f"r2_{M}_A1"] = float(1 - np.sum((yh - Ah @ thh) ** 2) / np.sum((yh - yh.mean()) ** 2))
            dKA = ForceConstantPotential(cs, thh).get_force_constants(scA).get_fc_array(order=2)
            res[f"F_{M}_A1"], res[f"fmin_{M}_A1"], res[f"qmin_{M}_A1"] = thermo(K + dKA, project=False)
            res[f"F_{M}_A1_proj"], res[f"fmin_{M}_A1_proj"], res[f"qmin_{M}_A1_proj"] = thermo(K + dKA, project=True)
            res[f"A1_label_vs_fc_reldiff_{M}"] = float(np.linalg.norm(Ah @ thh - np.concatenate([-np.einsum("iab,b->ia", dKA[:, x["atom"] - 1], np.array(x["displacement"])) for x in disp]).ravel())
                                                       / max(np.linalg.norm(Ah @ thh), 1e-30))
        # ---- oracle truncated to rc ----
        dKo = (K_dft - K) * (dist < rc)[..., None, None]
        for i in range(N): dKo[i, i] = -dKo[i].sum(axis=0) + dKo[i, i]
        for tag, KK in (("base_proj", K), ("R3_proj", K + dK3), ("oracle6_proj", K + dKo)):
            res[f"F_{M}_{tag}"], res[f"fmin_{M}_{tag}"], res[f"qmin_{M}_{tag}"] = thermo(KK, project=True)
        res[f"npar_{M}_R3"] = npr
    np.savez(out, **res)
    print(f"{out.stem}: N={N} rc={rc:.2f} npar_A1={res['npar_A1']} {time.time() - t0:.0f} s", flush=True)


for src in entries:
    out = out_dir / (Path(src).name.split(".")[0] + ".npz")
    if out.exists(): continue
    try:
        run(src, out)
    except Exception as exc:
        (out_dir / (out.stem + ".FAILED.txt")).write_text(repr(exc) + "\n"); print(out.stem, "FAILED", repr(exc), flush=True)
