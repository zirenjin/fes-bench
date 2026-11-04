"""PROTOCOL.md (relaxation channel), per entry: MLIP relaxation of the DFT unit cell, then harmonic phonons at the MLIP geometry.

Relaxation: FrechetCellFilter + FixSymmetry(symprec 1e-3), FIRE, fmax 1e-3 eV/A, <= 3000 steps, P = 0.
Phonons: stored supercell/primitive matrices, phonopy displacements 0.01 A (default +/- rule), symprec 1e-4, mesh 100.
Usage: relax_phonons.py <ckpt> <head> <out_dir> <entry.yaml.bz2> [more ...]   (model loaded once; resume-safe)
"""
import bz2, sys, time, traceback
from pathlib import Path
import numpy as np, torch, yaml, spglib
LOADER = getattr(yaml, "CSafeLoader", yaml.SafeLoader)
torch.serialization.add_safe_globals([slice])
from ase import Atoms  # noqa: E402
from ase.constraints import FixSymmetry  # noqa: E402
from ase.filters import FrechetCellFilter  # noqa: E402
from ase.optimize import FIRE  # noqa: E402
from deepmd.calculator import DP  # noqa: E402
from phonopy import Phonopy  # noqa: E402
from phonopy.structure.atoms import PhonopyAtoms  # noqa: E402

ckpt, head, out_dir = sys.argv[1:4]; entries = sys.argv[4:]; out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
TGRID = np.arange(0, 1501, 50, dtype=float)
calc = DP(ckpt, head=head); dp = calc.dp; tmap = dp.get_type_map()


def sg(at):
    return spglib.get_spacegroup((at.cell[:], at.get_scaled_positions(), at.numbers), symprec=1e-3)


def forces(cell, symbols, cart):
    return np.asarray(dp.eval(cart[None], np.asarray(cell).reshape(1, 9), np.array([tmap.index(s) for s in symbols]))[1]).reshape(-1, 3)


def run(src):
    d = yaml.load(bz2.open(src), Loader=LOADER)
    uc = d["unit_cell"]; pts = uc["points"]; sym = [p["symbol"] for p in pts]; n = len(pts)
    at = Atoms(sym, scaled_positions=[p["coordinates"] for p in pts], cell=uc["lattice"], pbc=True); at.calc = calc
    sg0, v0, frac0 = sg(at), at.get_volume() / n, at.get_scaled_positions()
    e_dftgeom = at.get_potential_energy() / n * 1000
    sc = d["supercell"]; lat = np.array(sc["lattice"]); r0 = forces(lat, [p["symbol"] for p in sc["points"]], np.array([p["coordinates"] for p in sc["points"]]) @ lat)
    at.set_constraint(FixSymmetry(at, symprec=1e-3))
    opt = FIRE(FrechetCellFilter(at), logfile=None); conv = bool(opt.run(fmax=1e-3, steps=3000))
    at.set_constraint()
    e_rel = at.get_potential_energy() / n * 1000; sg1 = sg(at)
    dfrac = (at.get_scaled_positions() - frac0 + 0.5) % 1 - 0.5; dmax = float(np.linalg.norm(dfrac @ at.cell[:], axis=1).max())
    cell = PhonopyAtoms(symbols=sym, cell=at.cell[:], scaled_positions=at.get_scaled_positions(), masses=[p["mass"] for p in pts])
    ph = Phonopy(cell, supercell_matrix=np.array(d["supercell_matrix"]), primitive_matrix=np.array(d.get("primitive_matrix", np.eye(3))), symprec=1e-4)
    ph.generate_displacements(distance=0.01)
    ph.forces = [forces(s.cell, s.symbols, s.positions) for s in ph.supercells_with_displacements]
    ph.produce_force_constants(); ph.run_mesh(100.0, with_eigenvectors=False); fmin = float(ph.mesh.frequencies.min())
    ph.run_thermal_properties(temperatures=TGRID); n_prim = len(ph.primitive)
    F = np.array(ph.get_thermal_properties_dict()["free_energy"]) * 1000 / 96.48533212 / n_prim
    return dict(T=TGRID, F_vib_relaxed=F, E_dftgeom_mlip=e_dftgeom, E_relaxed_mlip=e_rel, vol_dft=v0, vol_relaxed=at.get_volume() / n,
                dmax_A=dmax, sg_before=sg0, sg_after=sg1, converged=conv, nsteps=opt.nsteps, fmin_relaxed=fmin,
                R0_rms=float(np.sqrt(np.mean(r0 ** 2))), R0_max=float(np.linalg.norm(r0, axis=1).max()), n_disp=len(ph.supercells_with_displacements))


for src in entries:
    mid = Path(src).name.split(".")[0]; out = out_dir / f"{mid}.npz"
    if out.exists():
        continue
    t0 = time.time()
    try:
        r = run(src); np.savez(out, **r)
        print(f"{mid}: conv {r['converged']} steps {r['nsteps']} sg {r['sg_before']}->{r['sg_after']} dV {100 * (r['vol_relaxed'] / r['vol_dft'] - 1):+.2f}% "
              f"dE_relax {r['E_relaxed_mlip'] - r['E_dftgeom_mlip']:.2f} meV/atom R0 {r['R0_rms']:.4f} fmin {r['fmin_relaxed']:.3f} THz, {time.time() - t0:.0f} s", flush=True)
    except Exception:
        (out_dir / f"{mid}.FAILED.txt").write_text(traceback.format_exc()); print(f"{mid}: FAILED", flush=True)
