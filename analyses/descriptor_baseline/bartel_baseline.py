"""PROTOCOL.md: Bartel (2018) descriptor baseline on the database polymorph pairs.
Writes predictions/bartel_S{0,1}.npz in the scorer format and results/BARTEL.json.
Usage: python analyses/descriptor_baseline/bartel_baseline.py   (from the repository root)
"""
import csv, json, sys
from itertools import combinations
from math import gcd
from functools import reduce
from pathlib import Path
import numpy as np
from ase.formula import Formula
from ase.data import atomic_masses, atomic_numbers

HERE = Path(__file__).parent; DB = HERE.parents[1] / "database_pairs"; sys.path.insert(0, str(DB))
from score import Bench  # noqa: E402


def reduced_mass(formula):
    """Eq. 6 of Bartel et al. (2018), as implemented in pymatgen GibbsComputedStructureEntry._reduced_mass."""
    c = Formula(formula).count(); g = reduce(gcd, c.values()); c = {k: v // g for k, v in c.items()}
    if len(c) < 2: return None
    s = sum((a_i + a_j) * (m_i * m_j) / (m_i + m_j) for (e_i, a_i), (e_j, a_j) in combinations(c.items(), 2)
            for m_i, m_j in [(atomic_masses[atomic_numbers[e_i]], atomic_masses[atomic_numbers[e_j]])])
    return s / ((len(c) - 1) * sum(c.values()))


B = Bench(); E = B.E; T = B.T
V = np.array([float(r["vol_per_atom_A3"]) for r in E]); m = np.array([reduced_mass(r["formula"]) or np.nan for r in E])
slope = -2.48e-4 * np.log(V) - 8.94e-5 * m / V                                     # eV/atom/K; ln T and constant cancel within a formula
F = slope[:, None] * T[None, :] * 1000                                             # meV/atom
fmin = np.where(np.isfinite(m), 1.0, -np.inf)                                     # elemental formulas: reduced mass undefined -> invalid
F = np.where(np.isfinite(F), F, 0.0)
ref = np.load(DB / "predictions/reference_PBE.npz"); f0 = np.load(DB / "predictions/F0_DPA-3.1-3M_S0.npz")
(HERE / "predictions").mkdir(exist_ok=True); (HERE / "results").mkdir(exist_ok=True)
for S, E0, txt in (("S1", ref["E0"], "reference static energy"), ("S0", f0["E0"], "DPA-3.1 static energy")):
    np.savez(HERE / f"predictions/bartel_{S}.npz", ids=np.array(B.ids), T=T, fmin=fmin, E0=E0, F=F,
             info=json.dumps(dict(kind="baseline", model="Bartel 2018 descriptor (Eq. 4, 6)", setting=S, setting_text=txt)))

KEYS = ("crossing_success50", "crossing_hit50", "near_Tc_within_energy_budget50", "near_Tc_rms_e", "pair_dFth_rms_1000K", "pair_dS_rms_1000K",
        "wrong_ordering_mean", "ordering_success_mean", "failure_model_invalid")
res = dict(protocol_sha256=open(HERE / "PROTOCOL.sha256").read().split()[0], n_entries_elemental_invalid=int((~np.isfinite(m)).sum()),
           n_pairs_elemental=int(np.sum(~np.isfinite(m[B.pi]) & B.ref_valid)))
ge300 = B.ref_valid & B.cross & (B.tc >= 300)
for S in ("S1", "S0"):
    for sub in ("model_valid", "all"):
        mm = B.score(HERE / f"predictions/bartel_{S}.npz", sub)["metrics"]
        res[f"{S}|{sub}"] = {k: dict(value=mm[k]["value"], ci95=mm[k]["ci95"]) for k in KEYS if k in mm}
    mm = B.score(HERE / f"predictions/bartel_{S}.npz", ge300)["metrics"]
    res[f"{S}|crossings_ge300K"] = {k: dict(value=mm[k]["value"], ci95=mm[k]["ci95"]) for k in ("near_Tc_rms_e", "crossing_hit50", "near_Tc_within_energy_budget50") if k in mm}
    res[f"{S}|crossings_ge300K"]["n_pairs"] = int(ge300.sum())
(HERE / "results/BARTEL.json").write_text(json.dumps(res, indent=1) + "\n"); print(json.dumps(res, indent=1))
