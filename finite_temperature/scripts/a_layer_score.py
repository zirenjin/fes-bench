"""PROTOCOL.md (A layer) scoring: validity, pair dG RMS, crossings, eps50 budget, SiO2 lowest-phase accuracy, and the
label-using controls under the frozen temp_extrap / phase_lopo / system_loso splits. Truncated (invalid) curves are scored
only as a labelled diagnostic.
Usage: a_layer_score.py <data_dir> <fes-bench processed dir> <out.json>
"""
import csv, json, sys
from itertools import combinations
from pathlib import Path
import numpy as np

data, proc, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
GATE = -0.05
PH = {"ti": ["hcp", "bcc"], "zr": ["hcp", "bcc"], "hf": ["hcp", "bcc"], "sio2": ["quartz_beta", "cristobalite_beta", "tridymite_p63mmc"]}
TE_THR = json.load(open(proc / "splits_v2/temp_extrap.json"))["threshold_T_K"]
SPLIT_SHA = {k: json.load(open(proc / f"splits_v2/{k}.json")).get("sha256") for k in ("temp_extrap", "phase_lopo", "system_loso")}


def read_ref(sysname, tier):
    """Per-phase reference G (meV/atom) on the common grid of the tier."""
    G = {}
    for p in PH[sysname]:
        f = "reference_G_pbe.csv" if sysname == "sio2" and tier == "PBE" else "reference_G.csv"   # SiO2: pure PBE vs PBE-D3(BJ)
        rows = list(csv.DictReader(open(proc / sysname / p / f))); T = np.array([float(r["T_K"]) for r in rows]); g = np.array([float(r["G_eV_per_atom"]) for r in rows])
        G[p] = (T, g * 1000)
    Ts = [G[p][0] for p in G]; Tc = Ts[0]
    for t in Ts[1:]: Tc = np.intersect1d(Tc, t)
    return Tc, {p: np.interp(Tc, *G[p]) for p in G}


def model(M, geom, sysname, p, T):
    z = np.load(data / f"{M}_{geom}_{sysname}_{p}.npz")
    return z["E0_eV_atom"] * 1000 + np.interp(T, z["T"], z["F_meV_atom"]), float(z["fmin_THz"]), z


def roots(T, d):
    """Sign changes between consecutive non-zero points (exact zeros occur in the 0.1 meV-quantized metal tables; DEVIATIONS D1)."""
    nz = np.flatnonzero(d != 0); s = np.sign(d[nz]); k = np.flatnonzero(s[:-1] * s[1:] < 0); i, j = nz[k], nz[k + 1]
    return [(T[a] - d[a] * (T[b] - T[a]) / (d[b] - d[a]), int(np.sign(d[b] - d[a]))) for a, b in zip(i, j)]


def crossing_score(T, dref, dmod):
    rr, rm = roots(T, dref), roots(T, dmod); hits = []
    for tr, sgn in rr:
        w = (T >= tr - 25) & (T <= tr + 25); dS = abs(np.polyfit(T[w], dref[w], 1)[0]); eps = dS * 50
        e_at = float(np.interp(tr, T, dmod - dref)); cand = [tm for tm, s2 in rm if s2 == sgn and abs(tm - tr) <= 50]
        hits.append(dict(T_ref=float(tr), direction=sgn, eps50_meV=float(eps), e_at_Tx_meV=e_at, within_eps50=bool(abs(e_at) <= eps),
                         hit50=bool(cand), Tc_err_K=float(min(cand, key=lambda x: abs(x - tr)) - tr) if cand else None))
    matched = sum(1 for tm, s2 in rm if any(s2 == h["direction"] and abs(tm - h["T_ref"]) <= 50 for h in hits))
    return dict(reference_roots=hits, model_roots=[float(t) for t, _ in rm], false_roots=len(rm) - matched)


def pair_block(T, dref, dmod, mask=None):
    m = np.ones_like(T, bool) if mask is None else mask
    return dict(rms_meV=float(np.sqrt(np.mean((dmod - dref)[m] ** 2))), mean_err_meV=float(np.mean((dmod - dref)[m])), **crossing_score(T[m], dref[m], dmod[m]))


res = dict(protocol_sha256=open(data.parent / "PROTOCOL.sha256").read().split()[0], split_sha256=SPLIT_SHA, phases={}, tiers={})
TIERS = [("ti", "PBE"), ("zr", "PBE"), ("hf", "PBE"), ("sio2", "PBE"), ("sio2", "PBE-D3(BJ)")]
for M in ("F0", "MS"):
    for geom in ("G1", "G2"):
        for s in PH:
            for p in PH[s]:
                z = np.load(data / f"{M}_{geom}_{s}_{p}.npz")
                res["phases"][f"{M}|{geom}|{s}|{p}"] = dict(fmin_THz=float(z["fmin_THz"]), valid=bool(z["fmin_THz"] >= GATE), sg_before=str(z["sg_before"]), sg_after=str(z["sg_after"]),
                                                            vol_per_atom=float(z["vol_per_atom"]), relax_converged=bool(z["relax_converged"]) if "relax_converged" in z else None)

for s, tier in TIERS:
    T, Gr = read_ref(s, tier); key = f"{s}|{tier}"; blk = dict(T_range=[float(T[0]), float(T[-1])], n_T=len(T), pairs={})
    for a, b in combinations(PH[s], 2):
        dref = Gr[a] - Gr[b]; blk["pairs"][f"{a}-{b}"] = pb = dict(reference_roots=[float(t) for t, _ in roots(T, dref)])
        for M in ("F0", "MS"):
            for geom in ("G1", "G2"):
                ga, fa, _ = model(M, geom, s, a, T); gb, fb, _ = model(M, geom, s, b, T); dmod = ga - gb; valid = fa >= GATE and fb >= GATE
                r = dict(valid=bool(valid), fmin=[fa, fb], scored_as="prediction" if valid else "truncated-curve diagnostic only", zero_shot=pair_block(T, dref, dmod))
                # temp_extrap: offsets fitted on training temperatures of the same pair
                tr, te = T < TE_THR[s], T >= TE_THR[s]; e = dref - dmod
                if tr.sum() >= 2 and te.sum() >= 2:
                    c = e[tr].mean(); lin = np.polyfit(T[tr], e[tr], 1)
                    r["temp_extrap"] = dict(threshold_K=TE_THR[s], zero_shot_test=pair_block(T, dref, dmod, te),
                                            TE_c=pair_block(T, dref, dmod + c, te), TE_l=pair_block(T, dref, dmod + np.polyval(lin, T), te))
                pb[f"{M}|{geom}"] = r
        # SiO2 LOPO, reference-partner budget (the predicted-partner policy is identical to zero shot by construction)
    if s == "sio2":
        lopo = {}
        for M in ("F0", "MS"):
            for geom in ("G1", "G2"):
                Gm = {p: model(M, geom, s, p, T)[0] for p in PH[s]}; fm = {p: model(M, geom, s, p, T)[1] for p in PH[s]}
                for h in PH[s]:
                    src = [p for p in PH[s] if p != h]; c = np.mean([Gr[p] - Gm[p] for p in src], axis=0)
                    for p in src:
                        dref = Gr[p] - Gr[h]; dpred = Gr[p] - (Gm[h] + c)
                        lopo[f"{M}|{geom}|holdout={h}|partner={p}"] = dict(valid=bool(min(fm.values()) >= GATE), zero_shot=pair_block(T, dref, Gm[p] - Gm[h]), reference_partner=pair_block(T, dref, dpred))
        blk["phase_lopo_reference_partner"] = lopo
    res["tiers"][key] = blk

# system_loso among metals: transfer of the hcp-bcc pair error
loso = {}
for M in ("F0", "MS"):
    for geom in ("G1", "G2"):
        P = {}
        for s in ("ti", "zr", "hf"):
            T, Gr = read_ref(s, "PBE"); ga, fa, _ = model(M, geom, s, "hcp", T); gb, fb, _ = model(M, geom, s, "bcc", T)
            P[s] = (T, Gr["hcp"] - Gr["bcc"], ga - gb, fa >= GATE and fb >= GATE)
        for h in P:
            src = [s for s in P if s != h]; T, dref, dmod, vh = P[h]
            c = np.mean([np.mean(P[s][1] - P[s][2]) for s in src]); lin = np.polyfit(np.concatenate([P[s][0] for s in src]), np.concatenate([P[s][1] - P[s][2] for s in src]), 1)
            r = dict(valid=bool(vh and all(P[s][3] for s in src)), zero_shot=pair_block(T, dref, dmod), LOSO_c=pair_block(T, dref, dmod + c), LOSO_l=pair_block(T, dref, dmod + np.polyval(lin, T)))
            if h == "hf":
                ov = T <= 1527; r["overlap_T"] = dict(zero_shot=pair_block(T, dref, dmod, ov), LOSO_c=pair_block(T, dref, dmod + c, ov), LOSO_l=pair_block(T, dref, dmod + np.polyval(lin, T), ov))
            loso[f"{M}|{geom}|holdout={h}"] = r
res["system_loso"] = loso

# SiO2 lowest-phase accuracy
for tier in ("PBE", "PBE-D3(BJ)"):
    T, Gr = read_ref("sio2", tier); lo_ref = np.argmin(np.stack([Gr[p] for p in PH["sio2"]]), 0); acc = {}
    for M in ("F0", "MS"):
        for geom in ("G1", "G2"):
            lo = np.argmin(np.stack([model(M, geom, "sio2", p, T)[0] for p in PH["sio2"]]), 0); acc[f"{M}|{geom}"] = float(np.mean(lo == lo_ref))
    res["tiers"][f"sio2|{tier}"]["lowest_phase_accuracy"] = acc
out.write_text(json.dumps(res, indent=1) + "\n"); print("wrote", out)
