"""Diagnostic: on training-pool entries, compare F0 and F1 forces with DFT, split into the equilibrium-offset part
(mean force over displacements, per atom) and the displacement response (force minus that mean)."""
import bz2, sys
from pathlib import Path
import numpy as np, torch, yaml
torch.serialization.add_safe_globals([slice])
from deepmd.infer import DeepPot
LOADER = getattr(yaml, "CSafeLoader", yaml.SafeLoader)
A = "GenSi/project/fes-alexandria-polymorph-audit-20261030"
models = {"F0": DeepPot("GenSi/project/fes_experiment_pbe_d3bj/models/DPA-3.1-3M.pt", head="MP_traj_v024_alldata_mixu"),
          "F1": DeepPot("train_F1/model.ckpt-20000.pt")}
if Path("train_F1b/model.ckpt-20000.pt").exists():
    models["F1b"] = DeepPot("train_F1b/model.ckpt-20000.pt")
for mid in sys.argv[1:]:
    d = yaml.load(bz2.open(f"{A}/data/raw/{mid}.yaml.bz2"), Loader=LOADER)
    sc = d["supercell"]; lat = np.array(sc["lattice"]); cart = np.array([p["coordinates"] for p in sc["points"]]) @ lat
    sym = [p["symbol"] for p in sc["points"]]
    Fd = np.array([x["forces"] for x in d["displacements"]])
    out = [mid, f"n_disp {len(Fd)}", f"|F_dft| rms {np.sqrt((Fd**2).mean()):.4f}"]
    for k, dp in models.items():
        tm = dp.get_type_map(); ty = np.array([tm.index(s) for s in sym])
        Fm = []
        for x in d["displacements"]:
            pos = cart.copy(); pos[x["atom"] - 1] += np.array(x["displacement"])
            Fm.append(np.asarray(dp.eval(pos[None], lat.reshape(1, 9), ty)[1]).reshape(-1, 3))
        Fm = np.array(Fm); e = Fm - Fd
        eq = dp.eval(cart[None], lat.reshape(1, 9), ty)[1].reshape(-1, 3)
        resp_m, resp_d = Fm - eq[None], Fd  # DFT residual at equilibrium assumed ~0
        out.append(f"{k}: rmse {np.sqrt((e**2).mean()):.4f} eq-residual rms {np.sqrt((eq**2).mean()):.4f} response rmse {np.sqrt(((resp_m-resp_d)**2).mean()):.4f}")
    print(" | ".join(out), flush=True)
