"""FES-Bench v0.1: build the fixed manifest (entries, pairs, splits) and the prediction files from stored curves.

No model or DFT computation: everything is read from earlier, frozen outputs.
Usage: python build.py   (paths below are the frozen source outputs)
"""
import csv, itertools, json
from pathlib import Path
import numpy as np

HERE = Path(__file__).parent; OUT = HERE.parent
SRC = dict(
    groups=OUT / "fes-alexandria-polymorph-audit-20261030/metadata/selection/FINAL_GROUPS.csv",
    curves_F0=OUT / "fes-tc-sensitivity-20261102/data/curves_F0.npz",
    curves_MS=OUT / "fes-tc-sensitivity-20261102/data/curves_MS.npz",
    meta=HERE / "manifest/entry_meta_raw.csv",
    excl8=OUT / "fes-mattersim-replication-20261102/groups_excl8.csv",
    excl71=OUT / "fes-mattersim-replication-20261102/groups_excl71.csv",
    relax=OUT / "fes-relaxation-channel-20261101/results/RELAX_RESULTS.json",
    routes=OUT / "fes-correction-routes-20261102/results/OOF_G.npz")
GATE = -0.05; LO, HI = 50.0, 1500.0
(HERE / "manifest").mkdir(exist_ok=True); (HERE / "predictions").mkdir(exist_ok=True)

groups = {r["formula"]: r["mp_ids"].split(";") for r in csv.DictReader(open(SRC["groups"]))}
ids = [m for ms in groups.values() for m in ms]; form = {m: f for f, ms in groups.items() for m in ms}
forms = sorted(groups)


def load(p):
    z = np.load(p); k = {m: n for n, m in enumerate(z["ids"])}; o = [k[m] for m in ids]
    return dict(T=z["T"], E_dft=z["E_dft"][o], F_dft=z["F_dft"][o], E=z["E_mlip"][o], F=z["F_mlip"][o], fmin_dft=z["fmin_dft"][o], fmin=z["fmin_mlip"][o])


F0, MS = load(SRC["curves_F0"]), load(SRC["curves_MS"]); T = F0["T"]
assert np.allclose(F0["F_dft"], MS["F_dft"]) and np.allclose(F0["E_dft"], MS["E_dft"])
meta = {r["mp_id"]: r for r in csv.DictReader(open(SRC["meta"]))}
kept8 = {r["formula"] for r in csv.DictReader(open(SRC["excl8"]))}; kept71 = {r["formula"] for r in csv.DictReader(open(SRC["excl71"]))}
rx = json.load(open(SRC["relax"]))["excluded_ids"]; relax = {m: k for k, v in rx.items() for m in v}

# ---------- entries ----------
rows = []
for n, m in enumerate(ids):
    r = meta[m]; v = [F0["fmin_dft"][n], F0["fmin"][n], MS["fmin"][n]]
    rows.append(dict(mp_id=m, formula=form[m], spacegroup=r["spacegroup"], crystal_system=r["crystal_system"], n_atoms_primitive=r["n_prim"],
                     vol_per_atom_A3=r["vol_per_atom"], natoms_supercell=r["natoms_supercell"], n_displacements=r["n_displacements"],
                     reference="Alexandria PBE phonon benchmark (Loew et al. 2025)", reference_tier="PBE (VASP), fixed DFT geometry, harmonic, phonopy mesh 100",
                     license="CC BY 4.0",
                     fmin_dft_THz=f"{v[0]:.4f}", fmin_F0_THz=f"{v[1]:.4f}", fmin_MatterSim_THz=f"{v[2]:.4f}",
                     harmonic_valid_dft=int(v[0] >= GATE), harmonic_valid_F0=int(v[1] >= GATE), harmonic_valid_MatterSim=int(v[2] >= GATE),
                     harmonic_valid_common=int(min(v) >= GATE),
                     flag_U_scheme_tm_without_O_F=int(form[m] not in kept8), flag_potentially_magnetic=int(form[m] not in kept71),
                     relax_pipeline_F0=relax.get(m, "ok"), seconds_F0=r["seconds_F0"], seconds_MatterSim=r["seconds_MS"]))
with open(HERE / "manifest/ENTRIES.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

# ---------- pairs ----------
Gref = F0["E_dft"][:, None] + F0["F_dft"]; ix = {m: n for n, m in enumerate(ids)}; sel = (T >= LO) & (T <= HI)


def crossings(d):
    s = np.sign(d[sel]); Tt = T[sel]; dd = d[sel]; k = np.flatnonzero(s[:-1] * s[1:] < 0)
    return [float(Tt[i] - dd[i] * (Tt[i + 1] - Tt[i]) / (dd[i + 1] - dd[i])) for i in k]


def on_envelope(f, i, j, tc):
    """Is the pairwise reference crossing at tc also a boundary of the lower envelope of this formula's candidate phases?"""
    g = {m: np.interp(tc, T, Gref[ix[m]]) for m in groups[f]}; gij = 0.5 * (g[i] + g[j])
    return int(all(g[m] >= gij - 1e-6 for m in groups[f] if m not in (i, j)))


prow = []
for f, ms in groups.items():
    for i, j in itertools.combinations(ms, 2):
        d = Gref[ix[i]] - Gref[ix[j]]; c = crossings(d); tc = c[0] if c else None
        S = -np.interp(tc, T, np.gradient(d, T)) * 1000 if tc else None
        valid = int(rows[ix[i]]["harmonic_valid_common"] and rows[ix[j]]["harmonic_valid_common"])
        refv = int(rows[ix[i]]["harmonic_valid_dft"] and rows[ix[j]]["harmonic_valid_dft"])
        prow.append(dict(pair_id=f"{i}|{j}", formula=f, i=i, j=j, **{f"dG_ref_{t}K": f"{d[int(np.flatnonzero(T == t)[0])]:.4f}" for t in (300, 600, 1000, 1500)},
                         n_ref_crossings=len(c), Tc_ref_first_K="" if tc is None else f"{tc:.2f}",
                         dS_ref_at_Tc_ueV_per_atomK="" if S is None else f"{S:.4f}", eps50_meV="" if S is None else f"{abs(S) * 0.05:.5f}",
                         n_candidate_phases=len(ms), crossing_on_lower_envelope="" if tc is None else on_envelope(f, i, j, tc),
                         reference_valid=refv, harmonic_valid_common=valid))
with open(HERE / "manifest/PAIRS.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(prow[0])); w.writeheader(); w.writerows(prow)

# ---------- splits (outer folds of the correction-route protocol, seeds 1001-1003) ----------
srow = {f: dict(formula=f) for f in forms}
for rep in (1, 2, 3):
    rng = np.random.default_rng(1000 + rep); p = list(rng.permutation(forms))
    for k, fs in enumerate(np.array_split(p, 5)):
        for f in fs: srow[f][f"outer_fold_repeat{rep}"] = k
with open(HERE / "manifest/SPLITS.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["formula", "outer_fold_repeat1", "outer_fold_repeat2", "outer_fold_repeat3"]); w.writeheader(); w.writerows(srow.values())

# ---------- predictions ----------
def save(name, E, F, fmin, info, G=None):
    d = dict(ids=np.array(ids), T=T, fmin=fmin, info=json.dumps(info))
    if G is None:
        d.update(E0=E, F=F)
    else:
        d.update(G=G)
    np.savez_compressed(HERE / "predictions" / f"{name}.npz", **d)


save("reference_PBE", F0["E_dft"], F0["F_dft"], F0["fmin_dft"], dict(kind="reference"))
for S, use_ref_E in (("S0", False), ("S1", True)):
    tag = dict(S0="no reference information", S1="static dE0 taken from the reference (common control)")[S]
    E0F = F0["E_dft"] if use_ref_E else F0["E"]; E0M = MS["E_dft"] if use_ref_E else MS["E"]
    save(f"F0_DPA-3.1-3M_{S}", E0F, F0["F"], F0["fmin"], dict(kind="model", model="DPA-3.1-3M (MP_traj_v024_alldata_mixu head)", setting=S, setting_text=tag))
    save(f"MatterSim-v1.0.0-5M_{S}", E0M, MS["F"], MS["fmin"], dict(kind="model", model="MatterSim-v1.0.0-5M", setting=S, setting_text=tag))
    save(f"AVG_F0+MatterSim_{S}", (E0F + E0M) / 2, (F0["F"] + MS["F"]) / 2, np.minimum(F0["fmin"], MS["fmin"]),
         dict(kind="baseline", model="equal-weight average of the two models' G(T) (frozen default combination baseline)", setting=S, setting_text=tag))
# interventions (outer out-of-fold G, all three repeats) from the correction-route protocol; repeat 1 is the primary
z = np.load(SRC["routes"]); zi = {m: n for n, m in enumerate(z["ids"])}; o = [zi[m] for m in ids]
for M, fm in (("F0", F0["fmin"]), ("MS", MS["fmin"])):
    for S in ("S0", "S1"):
        for r in ("A", "B", "C"):
            for rep in (1, 2, 3):
                save(f"route{r}_{M}_{S}" + ("" if rep == 1 else f"_rep{rep}"), None, None, fm,
                     dict(kind="intervention", route=r, model=M, setting=S, repeat=rep, source=f"fes-correction-routes-20261102 OOF repeat {rep}"),
                     G=z[f"{rep}_{M}_{S}_{r}"][o])
# consistency: the D1 out-of-fold G equals the average built here (D1 has no fitted parameter)
for S in ("S0", "S1"):
    a = np.load(HERE / "predictions" / f"AVG_F0+MatterSim_{S}.npz"); assert np.allclose(a["E0"][:, None] + a["F"], z[f"1_D_{S}_D1"][o], atol=1e-9)
cr = [p for p in prow if p["n_ref_crossings"] > 0]
print("entries", len(rows), "pairs", len(prow), "crossing pairs", len(cr), "reference-valid pairs", sum(p["reference_valid"] for p in prow),
      "reference-valid crossing pairs", sum(p["reference_valid"] for p in cr), "crossings on lower envelope", sum(p["crossing_on_lower_envelope"] for p in cr),
      "crossings in 2-phase formulas", sum(p["n_candidate_phases"] == 2 for p in cr),
      "common-valid entries", sum(r["harmonic_valid_common"] for r in rows), "common-valid pairs", sum(p["harmonic_valid_common"] for p in prow))
