"""PROTOCOL.md (gate 1): local harmonic curvature corrections fitted to all archived reference displacement labels.

Per entry and model (F0, MatterSim): force constants from raw MLIP forces ("none_raw", the v0.1 baseline) and from
response forces F(u) - F(0) ("none_R0"); restricted corrections R1 (global scale), R2 (shell scale), R3 (shell scale
+ central springs) fitted by least squares to dF = F_ref(u) - [F_MLIP(u) - F_MLIP(0)] (F_ref(0) ~ 0 assumed);
an oracle K_DFT check. Each K goes through the same phonopy path (mesh 100, 0-1500 K / 50 K).
Usage: gate1.py <out_dir> <entry.yaml.bz2> [...]   (resume-safe; runs on GenSi)
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

    def fc_from(forces):
        ph = Phonopy(cell, supercell_matrix=scm, primitive_matrix=pm)
        ph.dataset = {"natom": N, "first_atoms": [dict(number=x["atom"] - 1, displacement=np.array(x["displacement"]), forces=np.asarray(f)) for x, f in zip(disp, forces)]}
        ph.produce_force_constants(); fc = ph.force_constants; assert fc.shape == (N, N, 3, 3), fc.shape; return fc

    def thermo(fc):
        ph = Phonopy(cell, supercell_matrix=scm, primitive_matrix=pm); ph.force_constants = fc
        ph.run_mesh(100.0, with_eigenvectors=False); fmin = float(ph.mesh.frequencies.min())
        ph.run_thermal_properties(temperatures=TG); n = len(ph.primitive)
        return np.array(ph.get_thermal_properties_dict()["free_energy"]) * 1000 / 96.48533212 / n, fmin

    # pair classes (minimum image; supercell images are >= 12 A apart, so one image within 6 A)
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
    K_dft = fc_from(F["DFT"]); res["F_DFT"], res["fmin_DFT"] = thermo(K_dft)
    for M in ("F0", "MS"):
        K_raw = fc_from(F[M]); K = fc_from([f - R0[M] for f in F[M]])
        res[f"F_{M}_none_raw"], res[f"fmin_{M}_none_raw"] = thermo(K_raw); res[f"F_{M}_none_R0"], res[f"fmin_{M}_none_R0"] = thermo(K)
        res[f"R0_max_{M}"] = float(np.abs(R0[M]).max())
        lab = [np.asarray(fd) - (np.asarray(fm) - R0[M]) for fd, fm in zip(F["DFT"], F[M])]   # dF labels
        # DEVIATIONS.md D1: antisymmetrise over +-u partners (removes even-order terms and the reference residual)
        part = {}
        for i, xi in enumerate(disp):
            for j, xj in enumerate(disp):
                if i != j and xi["atom"] == xj["atom"] and np.allclose(np.array(xi["displacement"]), -np.array(xj["displacement"])): part[i] = j
        lab = [(lab[i] - lab[part[i]]) / 2 if i in part else lab[i] for i in range(len(lab))]
        res[f"n_pm_paired_{M}"] = len(part)
        y = np.concatenate(lab).ravel()
        # descriptive: how much label variance the full force-constant difference K_DFT - K explains (upper bound for any dK)
        yfc = np.concatenate([-np.einsum("iab,b->ia", (K_dft - K)[:, x["atom"] - 1], np.array(x["displacement"])) for x in disp]).ravel()
        res[f"r2_fullfc_{M}"] = float(1 - np.sum((y - yfc) ** 2) / np.sum((y - y.mean()) ** 2))
        res[f"label_rms_{M}"] = float(np.sqrt(np.mean(y ** 2)))
        # basis columns: v[i] = B[i,a] u for each displacement; prediction dF = -sum_c theta_c v_c
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
        Ascale = design("scale"); Arad = design("radial")
        reps = {"R1": np.sum(Ascale, axis=1, keepdims=True), "R2": Ascale, "R3": np.hstack([Ascale, Arad])}
        for rep, A in reps.items():
            npar = A.shape[1]; res[f"npar_{M}_{rep}"] = npar
            if npar > 0.1 * len(y):
                res[f"F_{M}_{rep}"] = np.full(len(TG), np.nan); res[f"fmin_{M}_{rep}"] = np.nan; res[f"r2_{M}_{rep}"] = np.nan; continue
            AtA = A.T @ A; lam = 1e-8 * np.trace(AtA) / npar; th = np.linalg.solve(AtA + lam * np.eye(npar), A.T @ y)
            res[f"r2_{M}_{rep}"] = float(1 - np.sum((y - A @ th) ** 2) / np.sum((y - y.mean()) ** 2)); res[f"theta_{M}_{rep}"] = th
            dK = np.zeros_like(K)
            if rep == "R1": th_s, th_r = np.full(nc, th[0]), np.zeros(nc)
            elif rep == "R2": th_s, th_r = th, np.zeros(nc)
            else: th_s, th_r = th[:nc], th[nc:]
            off = cls >= 0; ts = np.where(off, th_s[np.clip(cls, 0, None)], 0.0); tr = np.where(off, th_r[np.clip(cls, 0, None)], 0.0)
            dK = ts[..., None, None] * K * off[..., None, None] - tr[..., None, None] * RR * off[..., None, None]
            for i in range(N): dK[i, i] = -dK[i].sum(axis=0) + dK[i, i]
            res[f"F_{M}_{rep}"], res[f"fmin_{M}_{rep}"] = thermo(K + dK)
    np.savez(out, **res)
    print(f"{out.stem}: N={N} disp={len(disp)} classes={nc} {time.time() - t0:.0f} s", flush=True)


for src in entries:
    out = out_dir / (Path(src).name.split(".")[0] + ".npz")
    if out.exists(): continue
    try:
        run(src, out)
    except Exception as exc:
        (out_dir / (out.stem + ".FAILED.txt")).write_text(repr(exc) + "\n"); print(out.stem, "FAILED", repr(exc), flush=True)
