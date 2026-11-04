"""PROTOCOL.md (A layer): zero-shot harmonic G(T) for the 9 frozen v2 phases, F0 and MS, geometries G1 (fixed-symmetry
relaxation) and G2 (frozen representative). Response forces F(u) - F(0), phonopy mesh 100, validity f_min >= -0.05 THz.
Usage: a_layer_phonons.py <fes-bench processed dir> <out_dir>   (resume-safe; runs on GenSi)
"""
import sys, time, json
from pathlib import Path
import numpy as np, torch, spglib
torch.serialization.add_safe_globals([slice])
from ase.io import read
from ase.constraints import FixSymmetry
from ase.filters import FrechetCellFilter
from ase.optimize import BFGS
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from deepmd.calculator import DP
from mattersim.forcefield import MatterSimCalculator

proc, out_dir = Path(sys.argv[1]), Path(sys.argv[2]); out_dir.mkdir(parents=True, exist_ok=True)
SC = {"hf/bcc": 4 * np.eye(3), "hf/hcp": np.diag([4, 4, 2]), "ti/bcc": 4 * np.eye(3), "ti/hcp": np.diag([4, 4, 3]),
      "zr/bcc": 3 * np.eye(3), "zr/hcp": np.diag([4, 4, 2]), "sio2/quartz_beta": 2 * np.eye(3),
      "sio2/cristobalite_beta": 2 * np.eye(3), "sio2/tridymite_p63mmc": 2 * np.eye(3)}
TG = np.arange(0, 2601, 1, dtype=float)
CALC = {"F0": lambda: DP(model="/GenSIvePFS/users/zirenj/fes_experiment_pbe_d3bj/models/DPA-3.1-3M.pt", head="MP_traj_v024_alldata_mixu"),
        "MS": lambda: MatterSimCalculator(load_path="/root/.local/mattersim/pretrained_models/mattersim-v1.0.0-5M.pth", device="cuda")}


def sg(a):
    return spglib.get_spacegroup((a.cell[:], a.get_scaled_positions(), a.numbers), symprec=1e-3)


def run(M, calc, key, geom):
    out = out_dir / f"{M}_{geom}_{key.replace('/', '_')}.npz"
    if out.exists(): return
    t0 = time.time(); a = read(proc / key / "structure.extxyz"); a.calc = calc; res = dict(sg_before=sg(a))
    if geom == "G1":
        a.set_constraint(FixSymmetry(a)); opt = BFGS(FrechetCellFilter(a), logfile=None)
        conv = opt.run(fmax=1e-3, steps=1000); a.set_constraint()
        res.update(relax_converged=bool(conv), relax_steps=int(opt.nsteps))
    res.update(sg_after=sg(a), cell=a.cell[:], vol_per_atom=a.get_volume() / len(a), E0_eV_atom=a.get_potential_energy() / len(a),
               fmax_residual=float(np.abs(a.get_forces()).max()))
    ph = Phonopy(PhonopyAtoms(symbols=a.get_chemical_symbols(), cell=a.cell[:], scaled_positions=a.get_scaled_positions()),
                 supercell_matrix=SC[key], primitive_matrix="auto")
    ph.generate_displacements(distance=0.01); s0 = ph.supercell
    from ase import Atoms
    def forces(pos):
        b = Atoms(s0.symbols, positions=pos, cell=s0.cell, pbc=True); b.calc = calc; return b.get_forces()
    F0 = forces(s0.positions)
    ph.forces = [forces(d.positions) - F0 for d in ph.supercells_with_displacements]
    ph.produce_force_constants(); ph.run_mesh(100.0, with_eigenvectors=False)
    res["fmin_THz"] = float(ph.mesh.frequencies.min())
    ph.run_thermal_properties(temperatures=TG); n = len(ph.primitive)
    res.update(T=TG, F_meV_atom=np.array(ph.get_thermal_properties_dict()["free_energy"]) * 1000 / 96.48533212 / n,
               n_prim=n, natom_sc=len(s0), n_disp=len(ph.supercells_with_displacements), seconds=time.time() - t0)
    np.savez(out, **res)
    print(f"{M} {geom} {key}: sg {res['sg_before']} -> {res['sg_after']}, E0 {res['E0_eV_atom']:.5f}, fmin {res['fmin_THz']:.3f} THz, "
          f"V/at {res['vol_per_atom']:.3f}, {res.get('relax_steps', 0)} steps, {res['seconds']:.0f} s", flush=True)


for M in ("F0", "MS"):
    calc = CALC[M]()
    for key in SC:
        for geom in ("G1", "G2"):
            try:
                run(M, calc, key, geom)
            except Exception as exc:
                (out_dir / f"{M}_{geom}_{key.replace('/', '_')}.FAILED.txt").write_text(repr(exc) + "\n"); print(M, geom, key, "FAILED", repr(exc), flush=True)
