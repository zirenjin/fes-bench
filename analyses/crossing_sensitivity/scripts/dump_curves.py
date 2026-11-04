"""Collect the per-entry harmonic curves of one MLIP result directory into a single npz (runs on GenSi).
Usage: dump_curves.py <mlip_dir> <out.npz>
"""
import glob, sys
from pathlib import Path
import numpy as np

mdir, out = sys.argv[1:3]
ids, Fd, Fm, Ed, Em, fd, fm = [], [], [], [], [], [], []
for p in sorted(glob.glob(mdir + "/*.npz")):
    z = np.load(p); ids.append(Path(p).stem)
    Fd.append(z["F_dft"]); Fm.append(z["F_mlip"]); Ed.append(float(z["E_dft_eV_atom"]) * 1000); Em.append(float(z["E_mlip_eV_atom"]) * 1000)
    fd.append(float(z["fmin_dft"])); fm.append(float(z["fmin_mlip"])); T = z["T"]
np.savez(out, ids=np.array(ids), T=T, F_dft=np.array(Fd), F_mlip=np.array(Fm), E_dft=np.array(Ed), E_mlip=np.array(Em),
         fmin_dft=np.array(fd), fmin_mlip=np.array(fm))
print(len(ids), "entries")
