"""FES-Bench v0.1 unified scorer: one prediction file -> all main metrics, each labelled with its weighting.

A prediction file (npz) holds ids, T (K), fmin (THz) and either E0 (meV/atom) + F (meV/atom, harmonic F incl. ZPE)
or G = E0 + F. Component metrics (single-phase, cancellation, static) need E0 + F; all others use G only.
Uncertainty: formula-cluster bootstrap, 2000 draws, fixed seed (the same draws for every prediction, so
differences between predictions are paired). Intervals are conditional on the stored predictions: they do not
cover retraining, re-tuning or development-stage selection, and they are not adjusted for multiple comparisons.

Validity (review M2). Pairs touching a reference entry below the harmonic gate have no trustworthy target and are
never scored. Among the remaining reference-valid pairs, a pair is model-valid if both entries pass the gate for this
prediction. Accuracy metrics use model-valid pairs ("model_valid"), or the pairs valid for DFT and both models
("common_valid", for paired comparisons); "all" keeps truncated-harmonic curves and is a diagnostic only.
Failure and success are reported over all reference-valid pairs: success = model-valid and correct.
Crossings (review M3) are pairwise free-energy crossings; `crossing_on_lower_envelope` marks the ones that are
also boundaries of the candidate-phase lower envelope.
Usage: python score.py <prediction.npz> [...]   -> scores/<name>.json
"""
import csv, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).parent
GATE, LO, HI, NB, SEED = -0.05, 50.0, 1500.0, 2000, 20261102
TS = (300, 600, 1000, 1500)


class Bench:
    def __init__(self):
        self.E = list(csv.DictReader(open(HERE / "manifest/ENTRIES.csv"))); self.P = list(csv.DictReader(open(HERE / "manifest/PAIRS.csv")))
        self.ids = [r["mp_id"] for r in self.E]; self.ix = {m: n for n, m in enumerate(self.ids)}
        self.forms = sorted({r["formula"] for r in self.E}); fi = {f: k for k, f in enumerate(self.forms)}
        self.pi = np.array([self.ix[p["i"]] for p in self.P]); self.pj = np.array([self.ix[p["j"]] for p in self.P])
        self.pf = np.array([fi[p["formula"]] for p in self.P]); self.ef = np.array([fi[r["formula"]] for r in self.E])
        ref = np.load(HERE / "predictions/reference_PBE.npz"); assert list(ref["ids"]) == self.ids
        self.T = ref["T"]; self.ref = dict(E0=ref["E0"], F=ref["F"], G=ref["E0"][:, None] + ref["F"], fmin=ref["fmin"])
        self.it = {t: int(np.flatnonzero(self.T == t)[0]) for t in TS}; self.sel = (self.T >= LO) & (self.T <= HI)
        self.dref = self.ref["G"][self.pi] - self.ref["G"][self.pj]
        self.tc = np.array([float(p["Tc_ref_first_K"]) if p["Tc_ref_first_K"] else np.nan for p in self.P]); self.cross = np.isfinite(self.tc)
        self.eps50 = np.array([float(p["eps50_meV"]) if p["eps50_meV"] else np.nan for p in self.P])
        self.common_valid = np.array([int(p["harmonic_valid_common"]) for p in self.P], bool)
        self.ref_valid = np.array([int(p["reference_valid"]) for p in self.P], bool)
        self.envelope = np.array([p["crossing_on_lower_envelope"] == "1" for p in self.P], bool)
        self.boot = np.random.default_rng(SEED).integers(0, len(self.forms), (NB, len(self.forms)))

    # ---------- aggregation helpers ----------
    def per_formula(self, v, group):
        ok = np.isfinite(v); n = np.bincount(group[ok], None, len(self.forms)); s = np.bincount(group[ok], v[ok], len(self.forms)); return s, n

    def fw_mean(self, v, group, mask=None):
        """Formula-weighted mean: mean over formulas of the within-formula mean (formulas without values are skipped)."""
        v = np.where(mask if mask is not None else True, v, np.nan); s, n = self.per_formula(v, group)
        pf = np.where(n > 0, s / np.where(n > 0, n, 1), np.nan)
        return float(np.nanmean(pf)), np.nanmean(pf[self.boot], axis=1)

    def pooled_rms(self, v, group, mask=None):
        """Pooled RMS over items (entries or pairs); the bootstrap resamples formulas."""
        v = np.where(mask if mask is not None else True, v, np.nan); s, n = self.per_formula(v ** 2, group)
        return float(np.sqrt(s.sum() / n.sum())), np.sqrt(s[self.boot].sum(1) / n[self.boot].sum(1))

    @staticmethod
    def pack(pt, draws, weighting, unit=""):
        lo, hi = np.nanpercentile(draws, [2.5, 97.5]); return dict(value=pt, ci95=[float(lo), float(hi)], weighting=weighting, unit=unit, _draws=draws)

    def crossings(self, d):
        s = np.sign(d[self.sel]); Tt = self.T[self.sel]; dd = d[self.sel]; k = np.flatnonzero(s[:-1] * s[1:] < 0)
        return [float(Tt[i] - dd[i] * (Tt[i + 1] - Tt[i]) / (dd[i + 1] - dd[i])) for i in k]

    # ---------- scoring ----------
    def score(self, path, subset="all"):
        z = np.load(path); assert list(z["ids"]) == self.ids; info = json.loads(str(z["info"]))
        G = z["G"] if "G" in z else z["E0"][:, None] + z["F"]; comp = "E0" in z
        fmin = z["fmin"]; mv = self.ref_valid & (fmin[self.pi] >= GATE) & (fmin[self.pj] >= GATE)
        pm = subset if isinstance(subset, np.ndarray) else dict(all=np.ones(len(self.P), bool), model_valid=mv, common_valid=self.common_valid)[subset]
        em = np.zeros(len(self.E), bool); em[self.pi[pm]] = True; em[self.pj[pm]] = True
        d = G[self.pi] - G[self.pj]; e = d - self.dref; it = self.it; m = {}
        # Q1 precision near the transition
        etc = np.array([np.interp(t, self.T, x) if np.isfinite(t) else np.nan for t, x in zip(self.tc, d)])
        cm = self.cross & pm
        m["near_Tc_rms_e"] = self.pack(*self.fw_rms(etc, cm), "formula-weighted RMS over reference-crossing pairs", "meV/atom")
        m["near_Tc_within_energy_budget50"] = self.pack(*self.fw_mean((np.abs(etc) <= self.eps50).astype(float), self.pf, cm),
                                                        "formula-weighted fraction with |e(Tx_ref)| <= |dS_ref| x 50 K (energy budget, not a temperature hit)", "")
        m["near_Tc_within_energy_budget100"] = self.pack(*self.fw_mean((np.abs(etc) <= 2 * self.eps50).astype(float), self.pf, cm), "formula-weighted fraction", "")
        hit50 = np.array([self.hit(d[k], self.tc[k], 50.0) if self.cross[k] else np.nan for k in range(len(self.P))])
        m["crossing_hit50"] = self.pack(*self.fw_mean(hit50, self.pf, cm), "formula-weighted fraction with a model crossing within +-50 K of the reference crossing", "")
        m["crossing_hit50_envelope_only"] = self.pack(*self.fw_mean(hit50, self.pf, cm & self.envelope), "as crossing_hit50, lower-envelope crossings only", "")
        m["near_Tc_median_e_over_eps50"] = dict(value=float(np.nanmedian(np.abs(etc[cm]) / self.eps50[cm])), weighting="pair median", unit="")
        hits = miss = false = 0; dT = []
        for k in np.flatnonzero(pm):
            rc, pc = self.crossings(self.dref[k]), self.crossings(d[k]); used = set()
            for c in rc:
                j = next((j for j, p in enumerate(pc) if j not in used and abs(p - c) <= 100), None)
                if j is None: miss += 1
                else: hits += 1; used.add(j)
            false += len(pc) - len(used)
            if rc and pc: dT.append(min(abs(p - rc[0]) for p in pc))
        m["events"] = dict(hits=hits, missed=miss, false=false, n_reference=hits + miss, median_abs_dTc_K=float(np.median(dT)) if dT else None,
                           weighting="pair counts, 50-1500 K, +-100 K matching")
        # Q2 capability profile
        for t in (300, 1000):
            m[f"pair_dG_rms_{t}K"] = self.pack(*self.pooled_rms(e[:, it[t]], self.pf, pm), "pooled over pairs", "meV/atom")
        m["pair_dG_rms_1000K_fw"] = self.pack(*self.fw_rms(e[:, it[1000]], pm), "formula-weighted RMS", "meV/atom")
        dS = -np.gradient(e, self.T, axis=1) * 1000
        m["pair_dS_rms_1000K"] = self.pack(*self.pooled_rms(dS[:, it[1000]], self.pf, pm), "pooled over pairs", "ueV/(atom K)")
        wr = {}
        for t in (600, 1000, 1500):
            nd = pm & (np.abs(self.dref[:, it[t]]) >= 1.0); w = (np.sign(d[:, it[t]]) != np.sign(self.dref[:, it[t]])).astype(float)
            wr[t] = self.fw_mean(w, self.pf, nd); m[f"wrong_ordering_{t}K"] = self.pack(*wr[t], "formula-weighted, |dG_ref| >= 1 meV/atom", "")
        m["wrong_ordering_mean"] = self.pack(float(np.mean([wr[t][0] for t in wr])), np.mean([wr[t][1] for t in wr], axis=0), "mean of the three formula-weighted rates", "")
        if comp:
            E0, F = z["E0"], z["F"]; R = self.ref
            dFth = {t: (F[:, it[t]] - F[:, 0]) - (R["F"][:, it[t]] - R["F"][:, 0]) for t in TS}
            m["single_dFth_rms_1000K"] = self.pack(*self.pooled_rms(dFth[1000], self.ef, em), "pooled over entries", "meV/atom")
            m["single_dFth_rms_1000K_fw"] = self.pack(*self.fw_rms_entries(dFth[1000], em), "formula-weighted RMS over entries", "meV/atom")
            pth = dFth[1000][self.pi] - dFth[1000][self.pj]
            m["pair_dFth_rms_1000K"] = self.pack(*self.pooled_rms(pth, self.pf, pm), "pooled over pairs (ZPE excluded)", "meV/atom")
            m["pair_dFth_rms_1000K_fw"] = self.pack(*self.fw_rms(pth, pm), "formula-weighted RMS (ZPE excluded)", "meV/atom")
            num = pth ** 2; den = dFth[1000][self.pi] ** 2 + dFth[1000][self.pj] ** 2
            sn, nn = self.per_formula(np.where(pm, num, np.nan), self.pf); sd, _ = self.per_formula(np.where(pm, den, np.nan), self.pf)
            cf = np.where((nn > 0) & (sd > 0), 1 - sn / np.where(sd > 0, sd, 1), np.nan)
            m["cancellation_c_1000K"] = self.pack(float(np.nanmean(cf)), np.nanmean(cf[self.boot], axis=1), "formula-weighted mean of per-formula c", "")
            de0 = (E0 - R["E0"]); pe0 = de0[self.pi] - de0[self.pj]
            m["pair_dE0_rms"] = self.pack(*self.pooled_rms(pe0, self.pf, pm), "pooled over pairs", "meV/atom")
        # tail and coverage
        a1000 = np.abs(e[pm, it[1000]])
        m["tail_p95_abs_e_1000K"] = dict(value=float(np.quantile(a1000, 0.95)), weighting="pair quantile", unit="meV/atom")
        m["tail_n_pairs_max_abs_e_gt20"] = dict(value=int(np.sum(np.abs(e[pm][:, self.sel]).max(1) > 20)), weighting="pair count (50-1500 K)", unit="pairs")
        m["decided_fraction_lofo_quantile_1000K"] = self.lofo(d[:, it[1000]], self.dref[:, it[1000]], pm)
        # failure and success over all reference-valid pairs (independent of the accuracy subset)
        R = self.ref_valid; fail = (~mv).astype(float)
        m["failure_model_invalid"] = self.pack(*self.fw_mean(fail, self.pf, R), "formula-weighted fraction of reference-valid pairs that are model-invalid", "")
        succ = {}
        for t in (600, 1000, 1500):
            nd = R & (np.abs(self.dref[:, it[t]]) >= 1.0); ok = (mv & (np.sign(d[:, it[t]]) == np.sign(self.dref[:, it[t]]))).astype(float)
            succ[t] = self.fw_mean(ok, self.pf, nd)
        m["ordering_success_mean"] = self.pack(float(np.mean([succ[t][0] for t in succ])), np.mean([succ[t][1] for t in succ], axis=0),
                                               "formula-weighted, over all reference-valid non-degenerate pairs: model-valid AND correctly ordered", "")
        tsucc = np.where(mv, np.nan_to_num(hit50), 0.0)
        m["crossing_success50"] = self.pack(*self.fw_mean(tsucc, self.pf, R & self.cross), "formula-weighted, over all reference-valid crossing pairs: model-valid AND crossing within +-50 K", "")
        m["validity"] = dict(entries_invalid_model=int(np.sum(fmin < GATE)), entries_invalid_reference=int(np.sum(self.ref["fmin"] < GATE)),
                             entries=len(self.ids), pairs_in_subset=int(pm.sum()), pairs_reference_valid=int(R.sum()),
                             pairs_model_invalid_among_reference_valid=int((R & ~mv).sum()), crossing_pairs_in_subset=int(cm.sum()),
                             pairs_common_valid=int(self.common_valid.sum()),
                             note="invalid = dense-mesh minimum frequency < -0.05 THz at the DFT geometry; phonopy integrates positive modes only")
        return dict(name=Path(path).stem, info=info, subset=subset if isinstance(subset, str) else "custom", metrics=m)

    def fw_rms(self, v, mask):
        s, n = self.per_formula(np.where(mask, v ** 2, np.nan), self.pf); pf = np.where(n > 0, s / np.where(n > 0, n, 1), np.nan)
        return float(np.sqrt(np.nanmean(pf))), np.sqrt(np.nanmean(pf[self.boot], axis=1))

    def fw_rms_entries(self, v, mask):
        s, n = self.per_formula(np.where(mask, v ** 2, np.nan), self.ef); pf = np.where(n > 0, s / np.where(n > 0, n, 1), np.nan)
        return float(np.sqrt(np.nanmean(pf))), np.sqrt(np.nanmean(pf[self.boot], axis=1))

    def hit(self, d, tc, tol):
        pc = self.crossings(d); return float(bool(pc) and min(abs(p - tc) for p in pc) <= tol)

    def lofo(self, dpred, dref, pm, alphas=(0.05,)):
        """Leave-one-formula-out empirical-quantile abstention (review M5): threshold = plain (1-alpha) quantile of |error|
        pooled over the other formulas' pairs. No finite-sample correction and no formula-level exchangeability argument,
        so alpha is a tuning level, not a guaranteed coverage or selective-error bound; 'coverage' = decided fraction."""
        nd = pm & (np.abs(dref) >= 1.0); err = np.abs(dpred - dref); out = {}
        for a in alphas:
            dec = np.zeros(len(dpred), bool)
            for f in np.unique(self.pf[nd]):
                tr, te = nd & (self.pf != f), nd & (self.pf == f); dec[te] = np.abs(dpred[te]) > np.quantile(err[tr], 1 - a)
            wrong = np.sign(dpred) != np.sign(dref)
            out[str(a)] = dict(coverage=float(dec[nd].mean()), selective_risk=float((dec & wrong)[nd].sum() / max(dec[nd].sum(), 1)),
                               weighting="pair-weighted over non-degenerate pairs")
        return out


def strip(o):
    if isinstance(o, dict): return {k: strip(v) for k, v in o.items() if not k.startswith("_")}
    if isinstance(o, (np.floating, np.integer)): return o.item()
    return o


if __name__ == "__main__":
    B = Bench(); (HERE / "scores").mkdir(exist_ok=True)
    for p in sys.argv[1:]:
        for sub in ("model_valid", "common_valid", "all"):
            r = B.score(p, sub); (HERE / "scores" / f"{r['name']}__{sub}.json").write_text(json.dumps(strip(r), indent=1) + "\n")
            mm = r["metrics"]; print(f"{r['name']:28s} {sub:12s} nearTc {mm['near_Tc_rms_e']['value']:.2f} eps50 {mm['near_Tc_within_energy_budget50']['value']:.2f} hit50 {mm['crossing_hit50']['value']:.2f} succ50 {mm['crossing_success50']['value']:.2f} fail {mm['failure_model_invalid']['value']:.2f} osucc {100 * mm['ordering_success_mean']['value']:.1f} "
                                     f"wrong {100 * mm['wrong_ordering_mean']['value']:.1f} pair1000 {mm['pair_dG_rms_1000K']['value']:.2f} p95 {mm['tail_p95_abs_e_1000K']['value']:.1f} "
                                     f"dec {mm['decided_fraction_lofo_quantile_1000K']['0.05']['coverage']:.2f}" + (f" single {mm['single_dFth_rms_1000K']['value']:.2f} pairth {mm['pair_dFth_rms_1000K']['value']:.2f} c {mm['cancellation_c_1000K']['value']:.2f}" if 'single_dFth_rms_1000K' in mm else ""))
