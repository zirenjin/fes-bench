"""PROTOCOL.md part 3: tail audit of the fixed entry list on the identical stored displacements (runs on GenSi).

Sides: DFT (stored forces), F0 (DPA-3.1-3M, MP_traj head), MatterSim-v1.0.0-5M. Same phonopy path as
mlip_vs_dft_phonons.py (produce_force_constants, mesh 100), plus the checks listed in PROTOCOL.md.
Usage: tail_audit.py <raw_dir> <stored_F0_dir> <stored_MS_dir> <out.json> <mp-id> [...]
"""
import bz2, json, sys
from pathlib import Path
import numpy as np, torch, yaml, spglib
torch.serialization.add_safe_globals([slice])
from deepmd.infer import DeepPot  # noqa: E402
from mattersim.forcefield import MatterSimCalculator  # noqa: E402
from ase import Atoms  # noqa: E402
from phonopy import Phonopy  # noqa: E402
from phonopy.structure.atoms import PhonopyAtoms  # noqa: E402

raw, sF0, sMS, out = map(Path, sys.argv[1:5]); ids = sys.argv[5:]
dp = DeepPot("/GenSIvePFS/users/zirenj/fes_experiment_pbe_d3bj/models/DPA-3.1-3M.pt", head="MP_traj_v024_alldata_mixu"); tmap = dp.get_type_map()
ms = MatterSimCalculator(load_path="/root/.local/mattersim/pretrained_models/mattersim-v1.0.0-5M.pth", device="cuda")
KB = 0.08617333262  # meV/K
H_OVER_K = 47.9924307  # h/k_B in K/THz
T = 1000.0; TG = np.arange(0, 1501, 50, dtype=float)


def mode_fth(nu):
    """Per-mode F(T) - F(0) in meV for real frequencies (THz)."""
    return KB * T * np.log1p(-np.exp(-H_OVER_K * nu / T))


res = {}
for m in ids:
    d = yaml.load(bz2.open(raw / f"{m}.yaml.bz2"), Loader=yaml.CSafeLoader)
    uc = d["unit_cell"]; pts = uc["points"]; sym = [p["symbol"] for p in pts]; lat = np.array(uc["lattice"]); frac = np.array([p["coordinates"] for p in pts])
    cell = PhonopyAtoms(symbols=sym, cell=lat, scaled_positions=frac, masses=[p["mass"] for p in pts])
    pm = np.array(d.get("primitive_matrix", np.eye(3))); scm = np.array(d["supercell_matrix"])
    nums = [list(dict.fromkeys(sym)).index(s) + 1 for s in sym]
    sg = spglib.get_spacegroup((lat, frac, nums), symprec=1e-3)
    probe = Phonopy(cell, supercell_matrix=scm, primitive_matrix=pm); sc = probe.supercell; L = sc.cell
    img = min(np.linalg.norm(np.array([i, j, k]) @ L) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1) if (i, j, k) != (0, 0, 0))
    info = dict(space_group_spglib=sg, space_group_yaml=d["space_group"]["type"], natoms_uc=len(sym), vol_per_atom=float(abs(np.linalg.det(lat)) / len(sym)),
                supercell_natoms=len(sc), supercell_lengths=[float(x) for x in np.linalg.norm(L, axis=1)], min_image_distance=float(img),
                n_displacements=len(d["displacements"]), nac_in_yaml=any(k in d for k in ("born", "nac", "dielectric")), phonopy_version_ref=d["phonopy"]["version"])
    # Alexandria's own processed free energy (phonopy kJ/mol per primitive cell) for a processing check
    info["alexandria_free_e_meV_atom"] = {str(t): f * 1000 / 96.48533212 / len(probe.primitive) for t, f in zip(d["temperatures"], d["free_e"])}
    pos0 = sc.positions; types = np.array([tmap.index(s) for s in sc.symbols], dtype=np.int32)
    forces = {"DFT": [np.array(x["forces"]) for x in d["displacements"]], "F0": [], "MS": []}
    for x in d["displacements"]:
        p = pos0.copy(); p[x["atom"] - 1] += np.array(x["displacement"])
        forces["F0"].append(np.asarray(dp.eval(p[None], L.reshape(1, 9), types)[1]).reshape(-1, 3))
        a = Atoms(list(sc.symbols), positions=p, cell=L, pbc=True); a.calc = ms; forces["MS"].append(a.get_forces())
    sides = {}; fc_dft = None
    for side, fl in forces.items():
        ph = Phonopy(cell, supercell_matrix=scm, primitive_matrix=pm)
        ph.dataset = {"natom": len(sc), "first_atoms": [dict(number=x["atom"] - 1, displacement=np.array(x["displacement"]), forces=f) for x, f in zip(d["displacements"], fl)]}
        ph.produce_force_constants(); fc = ph.force_constants.copy(); n_prim = len(ph.primitive)
        if side == "DFT": fc_dft = fc
        asr = float(np.abs(fc.sum(axis=1)).max()) if fc.shape[0] == fc.shape[1] else None

        def thermo(mesh):
            ph.run_mesh(mesh, with_eigenvectors=False); nu = ph.mesh.frequencies; w = ph.mesh.weights / ph.mesh.weights.sum()
            ph.run_thermal_properties(temperatures=TG); Fp = np.array(ph.get_thermal_properties_dict()["free_energy"]) * 1000 / 96.48533212 / n_prim
            return nu, w, float(Fp[TG == T][0] - Fp[0]), Fp
        nu, w, fth, Fp = thermo(100.0)
        W = np.broadcast_to(w[:, None], nu.shape); real = nu > 0.05
        per_mode = np.where(real, mode_fth(np.clip(nu, 1e-6, None)), 0.0)
        man = float((W * per_mode).sum() / n_prim); low = float((W * per_mode * (nu < 0.5)).sum() / n_prim)
        _, _, fth200, _ = thermo(200.0)
        ph.symmetrize_force_constants(); _, _, fth_sym, _ = thermo(100.0)
        sides[side] = dict(fmin=float(nu.min()), frac_imag=float(W[nu < -0.05].sum() / W.sum()), frac_below_0p5THz=float(W[nu < 0.5].sum() / W.sum()),
                           mean_ln_nu_real=float((W * np.log(np.clip(nu, 1e-6, None)) * real).sum() / (W * real).sum()),
                           Fth1000=fth, Fth1000_manual_real_modes=man, Fth1000_from_modes_below_0p5THz=low, Fth1000_mesh200=fth200, Fth1000_symmetrized=fth_sym,
                           asr_violation_max=asr, fc_rel_err_vs_dft=float(np.linalg.norm(fc - fc_dft) / np.linalg.norm(fc_dft)),
                           F300_meV_atom=float(Fp[TG == 300][0]), F600_meV_atom=float(Fp[TG == 600][0]))
    for side, sdir in (("F0", sF0), ("MS", sMS)):
        z = np.load(sdir / f"{m}.npz"); i = int(np.flatnonzero(z["T"] == T)[0])
        sides[side]["stored_Fth1000"] = float(z["F_mlip"][i] - z["F_mlip"][0]); sides["DFT"]["stored_Fth1000"] = float(z["F_dft"][i] - z["F_dft"][0])
    for side in ("F0", "MS"):
        s, r = sides[side], sides["DFT"]
        s["dFth1000"] = s["Fth1000"] - r["Fth1000"]; s["dFth1000_high_T_estimate"] = 3 * KB * T * (s["mean_ln_nu_real"] - r["mean_ln_nu_real"])
        s["dFth1000_from_modes_below_0p5THz"] = s["Fth1000_from_modes_below_0p5THz"] - r["Fth1000_from_modes_below_0p5THz"]
    res[m] = dict(info=info, sides=sides)
    print(m, sg, {k: round(v["Fth1000"], 2) for k, v in sides.items()}, flush=True)
out.write_text(json.dumps(res, indent=1) + "\n")
