"""PROTOCOL_1b.md: physical-consistency check of R3. Same forces, labels (D1) and R3 fit as gate1.py; adds the
assembly self-check, ASR / permutation / space-group errors, location of the minimum frequency, and one fixed
projection (symmetrize_force_constants(level=1) then symmetrize_force_constants_by_space_group()) applied identically
to the baseline (R0-subtracted) and to R3.
Usage: gate1b.py <out_dir> <entry.yaml.bz2> [...]   (resume-safe; runs on GenSi)
"""
import bz2, sys, time
from pathlib import Path
import numpy as np, torch, yaml
torch.serialization.add_safe_globals([slice])
from deepmd.infer import DeepPot  # noqa: E402
from mattersim.forcefield import MatterSimCalculator  # noqa: E402
from ase import Atoms  # noqa: E402
from phonopy import Phonopy  # noqa: E402
from phonopy.structure.atoms import PhonopyAtoms  # noqa: E402

out_dir = Path(sys.argv[1]); entries = sys.argv[2:]
dp = DeepPot("/GenSIvePFS/users/zirenj/fes_experiment_pbe_d3bj/models/DPA-3.1-3M.pt", head="MP_traj_v024_alldata_mixu"); tmap = dp.get_type_map()
ms = MatterSimCalculator(load_path="/root/.local/mattersim/pretrained_models/mattersim-v1.0.0-5M.pth", device="cuda")
TG = np.arange(0, 1501, 50, dtype=float); CUT = 6.0; SH = (1.15, 1.5)


def run(src, out):
    t0 = time.time(); d = yaml.load(bz2.open(src), Loader=yaml.CSafeLoader)
    uc = d["unit_cell"]; pts = uc["points"]; sym = [p["symbol"] for p in pts]
    cell = PhonopyAtoms(symbols=sym, cell=np.array(uc["lattice"]), scaled_positions=np.array([p["coordinates"] for p in pts]), masses=[p["mass"] for p in pts])
    pm = np.array(d.get("primitive_matrix", np.eye(3))); scm = np.array(d["supercell_matrix"])
    ph0 = Phonopy(cell, supercell_matrix=scm, primitive_matrix=pm); sc = ph0.supercell; L = sc.cell; pos0 = sc.positions; N = len(sc)
    ssym = list(sc.symbols); types = np.array([tmap.index(s) for s in ssym], dtype=np.int32); disp = d["displacements"]

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

    def thermo(fc, project):
        ph = new_ph(); ph.force_constants = fc.copy()
        if project:
            ph.symmetrize_force_constants(level=1, show_drift=False); ph.symmetrize_force_constants_by_space_group(show_drift=False)
        ph.run_mesh(100.0, with_eigenvectors=False); fr = ph.mesh.frequencies; k = np.unravel_index(np.argmin(fr), fr.shape)
        q = ph.mesh.qpoints[k[0]]; q = q - np.round(q); fmin = float(fr[k])
        ph.run_thermal_properties(temperatures=TG); n = len(ph.primitive)
        return np.array(ph.get_thermal_properties_dict()["free_energy"]) * 1000 / 96.48533212 / n, fmin, float(np.linalg.norm(q))

    def sym_errors(K):
        s = np.abs(K).max(); asr = float(np.abs(K.sum(axis=1)).max() / s); perm = float(np.abs(K - K.transpose(1, 0, 3, 2)).max() / s)
        ph = new_ph(); ph.force_constants = K.copy(); ph.symmetrize_force_constants_by_space_group(show_drift=False)
        sg = float(np.linalg.norm(K - ph.force_constants) / np.linalg.norm(K)); return asr, perm, sg

    fr = sc.scaled_positions; dfr = fr[None, :, :] - fr[:, None, :]; dfr -= np.round(dfr); rv = dfr @ L; dist = np.linalg.norm(rv, axis=-1)
    np.fill_diagonal(dist, np.inf); keys = {}; cls = -np.ones((N, N), int); dmin = {}
    for i in range(N):
        for j in range(N):
            if dist[i, j] < CUT:
                k = tuple(sorted((ssym[i], ssym[j]))); dmin[k] = min(dmin.get(k, np.inf), dist[i, j])
    for i in range(N):
        for j in range(N):
            if dist[i, j] < CUT:
                k = tuple(sorted((ssym[i], ssym[j]))); r = dist[i, j] / dmin[k]; s = 0 if r < SH[0] else (1 if r < SH[1] else 2)
                cls[i, j] = keys.setdefault((k, s), len(keys))
    rhat = rv / np.where(np.isfinite(dist), dist, 1.0)[..., None]; RR = np.einsum("ija,ijb->ijab", rhat, rhat)
    nc = len(keys); res = dict(T=TG, n_classes=nc, natom_sc=N, n_disp=len(disp))
    K_dft = fc_from(F["DFT"]); res["F_DFT"], res["fmin_DFT"], res["qmin_DFT"] = thermo(K_dft, False)
    res["asr_DFT"], res["perm_DFT"], res["sg_DFT"] = sym_errors(K_dft)
    for M in ("F0", "MS"):
        K = fc_from([f - R0[M] for f in F[M]])
        lab = [np.asarray(fd) - (np.asarray(fm) - R0[M]) for fd, fm in zip(F["DFT"], F[M])]
        part = {i: j for i, xi in enumerate(disp) for j, xj in enumerate(disp)
                if i != j and xi["atom"] == xj["atom"] and np.allclose(np.array(xi["displacement"]), -np.array(xj["displacement"]))}
        lab = [(lab[i] - lab[part[i]]) / 2 if i in part else lab[i] for i in range(len(lab))]; y = np.concatenate(lab).ravel()

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
        A = np.hstack([design("scale"), design("radial")]); npar = A.shape[1]
        AtA = A.T @ A; th = np.linalg.solve(AtA + 1e-8 * np.trace(AtA) / npar * np.eye(npar), A.T @ y)
        th_s, th_r = th[:nc], th[nc:]; off = cls >= 0
        ts = np.where(off, th_s[np.clip(cls, 0, None)], 0.0); tr = np.where(off, th_r[np.clip(cls, 0, None)], 0.0)
        dK = ts[..., None, None] * K * off[..., None, None] - tr[..., None, None] * RR * off[..., None, None]
        for i in range(N): dK[i, i] = -dK[i].sum(axis=0) + dK[i, i]
        # step 1: assembly self-check (design prediction vs explicit -dK u)
        expl = np.concatenate([-np.einsum("iab,b->ia", dK[:, x["atom"] - 1], np.array(x["displacement"])) for x in disp]).ravel()
        res[f"assembly_reldiff_{M}"] = float(np.linalg.norm(A @ th - expl) / max(np.linalg.norm(A @ th), 1e-30))
        Kc = K + dK; res[f"theta_{M}"] = th; res[f"npar_{M}"] = npar
        for tag, KK in (("base", K), ("R3", Kc), ("dK", dK)):
            if tag == "dK":
                s = np.abs(K).max(); res[f"asr_{M}_dK"] = float(np.abs(dK.sum(axis=1)).max() / s); res[f"perm_{M}_dK"] = float(np.abs(dK - dK.transpose(1, 0, 3, 2)).max() / s); continue
            res[f"asr_{M}_{tag}"], res[f"perm_{M}_{tag}"], res[f"sg_{M}_{tag}"] = sym_errors(KK)
        for tag, KK in (("base", K), ("R3", Kc)):
            for proj in (False, True):
                nm = f"{tag}{'_proj' if proj else ''}"; res[f"F_{M}_{nm}"], res[f"fmin_{M}_{nm}"], res[f"qmin_{M}_{nm}"] = thermo(KK, proj)
    np.savez(out, **res)
    print(f"{out.stem}: N={N} {time.time() - t0:.0f} s", flush=True)


for src in entries:
    out = out_dir / (Path(src).name.split(".")[0] + ".npz")
    if out.exists(): continue
    try:
        run(src, out)
    except Exception as exc:
        (out_dir / (out.stem + ".FAILED.txt")).write_text(repr(exc) + "\n"); print(out.stem, "FAILED", repr(exc), flush=True)
