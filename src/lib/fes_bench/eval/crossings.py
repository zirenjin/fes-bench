"""Deterministic one-to-one matching of one-dimensional crossing events."""
from __future__ import annotations

import math
from typing import Any


def match_root_events(reference: list[float], predicted: list[float],
                      max_distance_K: float | None = None) -> dict[str, Any]:
    """Maximize matches, then minimize total |Tc error| in monotone order.

    With no distance limit this preserves the legacy geometric notion of a
    hit anywhere within the evaluated temperature interval. A caller may use
    a prespecified tolerance; this helper never estimates it from test data.
    """
    if max_distance_K is not None and (not math.isfinite(max_distance_K) or max_distance_K < 0):
        raise ValueError('max_distance_K must be finite and nonnegative')
    if not all(math.isfinite(v) for v in [*reference, *predicted]):
        raise ValueError('crossing temperatures must be finite')
    rs=sorted(enumerate(reference),key=lambda a:a[1]);ps=sorted(enumerate(predicted),key=lambda a:a[1])
    states={(0,0):(0,0.,[])}
    def put(key,value):
        old=states.get(key)
        if old is None or (-value[0],value[1])<(-old[0],old[1]):states[key]=value
    for i in range(len(rs)+1):
        for j in range(len(ps)+1):
            n,c,pairs=states[i,j]
            if i<len(rs):put((i+1,j),(n,c,pairs))
            if j<len(ps):put((i,j+1),(n,c,pairs))
            if i<len(rs) and j<len(ps):
                distance=abs(rs[i][1]-ps[j][1])
                if max_distance_K is None or distance<=max_distance_K:
                    put((i+1,j+1),(n+1,c+distance,pairs+[(rs[i][0],ps[j][0])]))
    pairs=states[len(rs),len(ps)][2]
    return {'matches':[{'reference_index':i,'predicted_index':j,'Tc_error_K':float(predicted[j]-reference[i])} for i,j in pairs],
            'n_reference_events':len(reference),'n_predicted_events':len(predicted),
            'n_matched_events':len(pairs),'missed_crossings':len(reference)-len(pairs),
            'false_crossings':len(predicted)-len(pairs)}
