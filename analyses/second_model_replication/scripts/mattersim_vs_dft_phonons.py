"""Same-geometry DFT vs MLIP harmonic thermodynamics for one Alexandria PBE entry: MatterSim variant of mlip_vs_dft_phonons.py.

Only the force/energy evaluator changes (MatterSimCalculator on the identical displaced supercells and unit cell); phonopy path,
mesh, temperature grid and output keys are unchanged, so the original scorers apply as they are.

Both sides go through the identical phonopy path: the stored displaced supercells, forces -> force constants
(produce_force_constants), mesh density 100, thermal properties on a fixed temperature grid, per primitive cell.
DFT forces are read from the file; MLIP forces and the unit-cell static energy come from DPA-3.1-3M with the
given head. No relaxation: geometries are the DFT ones, so this isolates PES-curvature and static-energy errors
from MLIP geometry errors (a separate experiment).
Usage: mlip_vs_dft_phonons.py <ckpt> <head> <out_dir> <entry.yaml.bz2> [more entries ...]
The model is loaded once per process; entries whose output exists are skipped (resume-safe).
"""
import bz2, resource, sys, time
import numpy as np, torch, yaml
LOADER = getattr(yaml, "CSafeLoader", yaml.SafeLoader)  # libyaml when available; same safe semantics
torch.serialization.add_safe_globals([slice])
from ase import Atoms  # noqa: E402
from mattersim.forcefield import MatterSimCalculator  # noqa: E402
from phonopy import Phonopy  # noqa: E402
from phonopy.structure.atoms import PhonopyAtoms  # noqa: E402
from ase.data import atomic_numbers  # noqa: E402

from pathlib import Path  # noqa: E402
ckpt, head, out_dir = sys.argv[1:4]; entries = sys.argv[4:]
TGRID = np.arange(0, 1501, 50, dtype=float)
calc = MatterSimCalculator(load_path=ckpt, device="cuda")   # head argument unused (kept for the same CLI)


def ev(symbols, cell, cart):
    a = Atoms(list(symbols), positions=cart, cell=np.asarray(cell), pbc=True); a.calc = calc
    return a.get_potential_energy(), a.get_forces()


def run_entry(src, out):
  t0 = time.time()
  d = yaml.load(bz2.open(src), Loader=LOADER)
  uc = d["unit_cell"]; pts = uc["points"]; sym = [p["symbol"] for p in pts]
  cell = PhonopyAtoms(symbols=sym, cell=np.array(uc["lattice"]), scaled_positions=np.array([p["coordinates"] for p in pts]),
                      masses=[p["mass"] for p in pts])
  pm = np.array(d.get("primitive_matrix", np.eye(3)))


  def phonon(forces_list):
      ph = Phonopy(cell, supercell_matrix=np.array(d["supercell_matrix"]), primitive_matrix=pm)
      sc = ph.supercell
      ref = np.array([p["coordinates"] for p in d["supercell"]["points"]])
      assert np.allclose((sc.scaled_positions - ref + 0.5) % 1 - 0.5, 0, atol=1e-6), "supercell ordering differs"
      ph.dataset = {"natom": len(sc), "first_atoms": [dict(number=x["atom"] - 1, displacement=np.array(x["displacement"]), forces=np.array(f))
                                                       for x, f in zip(d["displacements"], forces_list)]}
      ph.produce_force_constants()
      ph.run_mesh(100.0, with_eigenvectors=False)
      fmin = float(ph.mesh.frequencies.min())
      ph.run_thermal_properties(temperatures=TGRID)
      tp = ph.get_thermal_properties_dict(); n_prim = len(ph.primitive)  # dict API works in phonopy 2.x and 4.x
      to_meV_atom = 1000 / 96.48533212 / n_prim
      return dict(F=np.array(tp["free_energy"]) * to_meV_atom, fmin=fmin, n_prim=n_prim, fc=ph.force_constants, sc=ph.supercell)


  dft = phonon([x["forces"] for x in d["displacements"]])
  sc = dft["sc"]; sc_sym = sc.symbols; sc_cart = sc.positions; sc_cell = sc.cell
  mlip_forces = []
  for x in d["displacements"]:
      pos = sc_cart.copy(); pos[x["atom"] - 1] += np.array(x["displacement"])
      e, f = ev(sc_sym, sc_cell, pos)
      mlip_forces.append(np.asarray(f).reshape(-1, 3))
  ml = phonon(mlip_forces)
  e_uc, _ = ev(sym, np.array(cell.cell), np.array(cell.positions))
  n = len(sym)
  np.savez(out, T=TGRID, F_dft=dft["F"], F_mlip=ml["F"], fmin_dft=dft["fmin"], fmin_mlip=ml["fmin"],
           E_dft_eV_atom=d["energy"] / n, E_mlip_eV_atom=float(np.asarray(e_uc).ravel()[0]) / n,
           fc_rel_err=float(np.linalg.norm(ml["fc"] - dft["fc"]) / np.linalg.norm(dft["fc"])),
           natoms_supercell=len(sc_sym), n_displacements=len(d["displacements"]), head=head, seconds=time.time() - t0)
  print(f"{src.split('/')[-1]}: sc {len(sc_sym)} atoms, {len(d['displacements'])} disp, dE0 {(float(np.asarray(e_uc).ravel()[0]) - d['energy']) / n * 1000:+.1f} meV/atom (absolute), "
        f"ZPE dft {dft['F'][0]:.2f} mlip {ml['F'][0]:.2f}, F600 dft {dft['F'][12]:.2f} mlip {ml['F'][12]:.2f} meV/atom, fmin {dft['fmin']:.3f}/{ml['fmin']:.3f} THz, "
        f"{time.time() - t0:.0f} s, peak RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6:.2f} GB", flush=True)


for src in entries:
    out = Path(out_dir) / (Path(src).name.split('.')[0] + '.npz')
    if out.exists():
        continue
    try:
        run_entry(src, out)
    except Exception as exc:  # record and continue; failures are reported, not silently dropped
        (Path(out_dir) / (out.stem + '.FAILED.txt')).write_text(repr(exc) + "\n")
        print(out.stem, 'FAILED', repr(exc), flush=True)
