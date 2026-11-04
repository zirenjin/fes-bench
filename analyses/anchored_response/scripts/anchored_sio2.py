"""PROTOCOL.md dataset 3: SiO2 ds62 frozen generic-DPA readout (F0), beta-cristobalite held out; arms none / Tc(crist-quartz) / T0ref(851 K).

Rebuilds each arm-B cristobalite model's predicted curves with the unchanged run_abc_generic.py functions (features, basis).
Sources are reference-labelled inputs with b = 0 (their readout fits are used as in the original scorer); only b_crist is set.
Usage (WSL): anchored_sio2.py <run_F0 dir> <out.json>
"""
import importlib.util, json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
from anchor_common import crossings, match_events, tc_anchor, side_correct  # noqa: E402

run, out = Path(sys.argv[1]), Path(sys.argv[2])
spec = importlib.util.spec_from_file_location("rag", run / "run_abc_generic.py"); R = importlib.util.module_from_spec(spec); spec.loader.exec_module(R)
emb = np.load(run / "generic_checkpoint_embeddings.npz"); data = {p: R.read(p) for p in R.PH}; T = data[R.PH[0]][0]; B = R.basis(T)
truth = {p: data[p][1] * 1000 for p in R.PH}; TGT = "cristobalite_beta"; POOL = ["quartz_beta", "tridymite_p63mmc", TGT]
TE = [1000, 1500, 2000]; ti = {t: int(np.argmin(abs(T - t))) for t in TE}
res = dict(window_K=[float(T[0]), float(T[-1])], models={})
for f in sorted((run / "models").glob(f"{TGT}__*__B.npz")):
    _, cfg, scaling, _ = f.stem.split("__")
    src = [p for p in R.PH if p != TGT] if cfg == "high_all_three" else cfg[len("low_"):].split("+")
    z = R.features(emb, src, "B", scaling)[0]; w = np.load(f)["W"]
    pred = {p: B @ (z[p] @ w) * 1000 for p in R.PH}
    dref_q, d0_q = truth[TGT] - truth["quartz_beta"], pred[TGT] - pred["quartz_beta"]
    tc, _ = tc_anchor(T, dref_q, d0_q, T[0])
    bias = {"none": 0.0, "Tc": -float(np.interp(tc, T, d0_q)) if tc is not None else None, "T0ref": float(dref_q[0] - d0_q[0])}
    blk = dict(sources=src, scaling=scaling, provided_Tc=tc)
    for arm, b in bias.items():
        if b is None:
            blk[arm] = "not applicable (no reference crossing)"; continue
        g = dict(pred); g[TGT] = pred[TGT] + b; a = {}
        for other in [p for p in R.PH if p != TGT]:
            dref, d = truth[TGT] - truth[other], g[TGT] - g[other]
            a[other] = dict(rms_e=float(np.sqrt(np.mean((d - dref) ** 2))), dS_rms_ueV=float(np.sqrt(np.mean((np.gradient(d - dref, T) * 1000) ** 2))),
                            sign_correct_at={str(t): bool(np.sign(d[ti[t]]) == np.sign(dref[ti[t]])) for t in TE},
                            events=match_events(crossings(T, dref), crossings(T, d), provided=tc if (arm == "Tc" and other == "quartz_beta") else None))
            if arm == "Tc" and other == "quartz_beta":
                a[other]["side_correct"] = side_correct(T, dref, d, tc, T[0], T[-1])
        stp = np.array(POOL)[np.argmin([g[p] for p in POOL], axis=0)]; str_ = np.array(POOL)[np.argmin([truth[p] for p in POOL], axis=0)]
        pe = [float(T[k + 1]) for k in range(len(T) - 1) if stp[k] != stp[k + 1]]; re_ = [float(T[k + 1]) for k in range(len(T) - 1) if str_[k] != str_[k + 1]]
        a["envelope"] = dict(stable_accuracy=float(np.mean(stp == str_)), events=match_events(re_, pe, provided=tc if arm == "Tc" else None),
                             predicted=[(stp[k], stp[k + 1], float(T[k + 1])) for k in range(len(T) - 1) if stp[k] != stp[k + 1]])
        blk[arm] = a
    res["models"][f"{cfg}__{scaling}"] = blk
out.write_text(json.dumps(res, indent=1) + "\n"); print(json.dumps(res, indent=1)[:3000])
