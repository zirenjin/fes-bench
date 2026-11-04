"""Shared pair-level helpers for PROTOCOL.md (anchored temperature-response test). Units: meV/atom, K."""
import numpy as np


def crossings(T, d):
    """All zero crossings of d(T) (linear interpolation), in grid order."""
    s = np.sign(d); k = np.flatnonzero(s[:-1] * s[1:] < 0)
    return [float(T[i] - d[i] * (T[i + 1] - T[i]) / (d[i + 1] - d[i])) for i in k]


def match_events(ref, pred, provided=None, tol=100.0):
    """Greedy +-tol matching of predicted to reference crossings; a provided crossing is removed from both lists first."""
    ref = [r for r in ref if provided is None or abs(r - provided) > 1e-6]
    pred = [p for p in pred if provided is None or abs(p - provided) > tol]
    hits, used = 0, set()
    for r in ref:
        j = next((j for j, p in enumerate(pred) if j not in used and abs(p - r) <= tol), None)
        if j is not None:
            hits += 1; used.add(j)
    return dict(reference=len(ref), predicted=len(pred), hits=hits, missed=len(ref) - hits, false=len(pred) - hits)


def tc_anchor(T, d_ref, d_pred, lo):
    """Tc arm: first reference crossing at T >= lo; returns (Tc, corrected prediction) or (None, None) if not applicable."""
    sel = T >= lo; c = crossings(T[sel], d_ref[sel])
    if not c:
        return None, None
    tc = c[0]
    return tc, d_pred - np.interp(tc, T, d_pred)


def side_correct(T, d_ref, d_corr, tc, lo, hi, gap=100.0):
    m = (T >= lo) & (T <= hi) & (np.abs(T - tc) >= gap)
    return float(np.mean(np.sign(d_corr[m]) == np.sign(d_ref[m]))) if m.any() else None
