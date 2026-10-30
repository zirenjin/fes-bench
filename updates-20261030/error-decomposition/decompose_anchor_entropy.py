#!/usr/bin/env python3
"""Post-fit phase-pair anchor and integrated-entropy error decomposition."""
import csv
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
T0 = 1542.7087047758675
WINDOW = 21
ARMS = {
    'B_DPA': 'B_DPA',
    'C_DPA_plus_volume_equal_group': 'C_DPA_plus_volume_equal_group',
}
# Orient each gap as first phase minus second.
PAIRS = [
    ('cristobalite_beta', 'quartz_beta'),
    ('tridymite_p63mmc', 'quartz_beta'),
    ('cristobalite_beta', 'tridymite_p63mmc'),
]

def read_csv(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f))

def basis(t):
    x = np.asarray(t, dtype=float) / 1000.0
    return np.column_stack([np.ones_like(x), x, x*np.log(x), x*x])

def local_slope_error(t, eg, n=WINDOW):
    out = np.empty_like(t, dtype=float)
    for i, ti in enumerate(t):
        ix = np.argsort(np.abs(t-ti))[:n]
        c, *_ = np.linalg.lstsq(basis(t[ix]), eg[ix], rcond=None)
        x = ti/1000.0
        dGdT_meV_K = (c[1] + (np.log(x)+1)*c[2] + 2*x*c[3]) / 1000.0
        out[i] = -1000.0*dGdT_meV_K  # -d(error G)/dT in µeV/(atom K)
    return out

def cumulative_entropy_component(t, ds_err):
    # Insert T0 and integrate forward using trapezoids; convert µeV to meV.
    tx = np.unique(np.r_[t, T0])
    sx = np.interp(tx, t, ds_err)
    increments = 0.5*(sx[1:]+sx[:-1])*np.diff(tx)*1e-3
    cum = np.r_[0.0, np.cumsum(increments)]
    c0 = np.interp(T0, tx, cum)
    return tx, -(cum-c0)

def interpolate(t, y, x):
    return float(np.interp(x, t, y))

phase_rows = read_csv(ROOT/'TARGET_PHASE_CURVES.csv')
by_arm_phase = {}
for r in phase_rows:
    if r['source_arm'] != 'three_source':
        continue
    by_arm_phase.setdefault(r['representation'], {}).setdefault(r['phase'], []).append(r)
for reps in by_arm_phase.values():
    for ph in reps:
        reps[ph] = sorted(reps[ph], key=lambda r: float(r['T_K']))

summary, curve_out = [], []
for rep in ARMS:
    phases = by_arm_phase[rep]
    t = np.array([float(r['T_K']) for r in phases['quartz_beta']])
    for first, second in PAIRS:
        a, b = phases[first], phases[second]
        assert np.array_equal(t, np.array([float(r['T_K']) for r in a]))
        ref = np.array([float(x['reference_G_meV_atom']) for x in a]) - np.array([float(x['reference_G_meV_atom']) for x in b])
        pred = np.array([float(x['predicted_G_meV_atom']) for x in a]) - np.array([float(x['predicted_G_meV_atom']) for x in b])
        eg = pred-ref
        ds_err = local_slope_error(t, eg)
        anchor = interpolate(t, eg, T0)
        tx, entropy_component_grid = cumulative_entropy_component(t, ds_err)
        entropy_component = np.interp(t, tx, entropy_component_grid)
        reconstructed = anchor + entropy_component
        closure = eg-reconstructed
        for i in range(len(t)):
            curve_out.append([rep, first, second, float(t[i]), float(eg[i]), float(ds_err[i]), anchor, float(entropy_component[i]), float(reconstructed[i]), float(closure[i])])
        # Summaries quantify how B->C changes anchor and integrated response separately.
        endpoints = [T0, 851.0, 2499.0]
        records = {}
        for temp in endpoints:
            records[temp] = (interpolate(t, eg, temp), interpolate(t, entropy_component, temp), interpolate(t, ds_err, temp))
        summary.append({
            'representation': rep, 'gap_first_minus_second': f'{first}-{second}',
            'T0_K': T0, 'anchor_error_at_T0_meV_atom': anchor,
            'deltaG_error_MAE_meV_atom': float(np.mean(np.abs(eg))),
            'deltaG_error_bias_meV_atom': float(np.mean(eg)),
            'deltaS_error_MAE_micro_eV_atom_K': float(np.mean(np.abs(ds_err))),
            'deltaS_error_bias_micro_eV_atom_K': float(np.mean(ds_err)),
            'at_851K_deltaG_error_meV_atom': records[851.0][0],
            'at_851K_integrated_entropy_contribution_meV_atom': records[851.0][1],
            'at_851K_deltaS_error_micro_eV_atom_K': records[851.0][2],
            'at_T0_deltaG_error_meV_atom': records[T0][0],
            'at_T0_integrated_entropy_contribution_meV_atom': records[T0][1],
            'at_T0_deltaS_error_micro_eV_atom_K': records[T0][2],
            'at_2499K_deltaG_error_meV_atom': records[2499.0][0],
            'at_2499K_integrated_entropy_contribution_meV_atom': records[2499.0][1],
            'at_2499K_deltaS_error_micro_eV_atom_K': records[2499.0][2],
            'decomposition_max_abs_closure_meV_atom': float(np.max(np.abs(closure))),
            'decomposition_RMSE_closure_meV_atom': float(np.sqrt(np.mean(closure**2))),
        })

# Write stable, self-describing tables.
fields = list(summary[0])
with (ROOT/'ANCHOR_ENTROPY_PAIR_SUMMARY.csv').open('w', newline='') as f:
    w=csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(summary)
with (ROOT/'ANCHOR_ENTROPY_ERROR_CURVES.csv').open('w', newline='') as f:
    w=csv.writer(f); w.writerow(['representation','gap_first_minus_second','T_K','deltaG_error_meV_atom','deltaS_error_local21_micro_eV_atom_K','anchor_error_at_T0_meV_atom','integrated_entropy_contribution_from_T0_meV_atom','reconstructed_deltaG_error_meV_atom','closure_residual_meV_atom']); w.writerows(curve_out)

# Paired C minus B decomposition table.
lookup={(r['representation'],r['gap_first_minus_second']):r for r in summary}
paired=[]
for first,second in PAIRS:
    pair=f'{first}-{second}'
    b=lookup[('B_DPA',pair)]; c=lookup[('C_DPA_plus_volume_equal_group',pair)]
    paired.append({
      'gap_first_minus_second':pair,
      'C_minus_B_anchor_error_at_T0_meV_atom':float(c['anchor_error_at_T0_meV_atom'])-float(b['anchor_error_at_T0_meV_atom']),
      'C_minus_B_deltaG_error_MAE_meV_atom':float(c['deltaG_error_MAE_meV_atom'])-float(b['deltaG_error_MAE_meV_atom']),
      'C_minus_B_deltaG_error_bias_meV_atom':float(c['deltaG_error_bias_meV_atom'])-float(b['deltaG_error_bias_meV_atom']),
      'C_minus_B_deltaS_error_MAE_micro_eV_atom_K':float(c['deltaS_error_MAE_micro_eV_atom_K'])-float(b['deltaS_error_MAE_micro_eV_atom_K']),
      'C_minus_B_deltaS_error_bias_micro_eV_atom_K':float(c['deltaS_error_bias_micro_eV_atom_K'])-float(b['deltaS_error_bias_micro_eV_atom_K']),
      'C_minus_B_integrated_entropy_contribution_at_851K_meV_atom':float(c['at_851K_integrated_entropy_contribution_meV_atom'])-float(b['at_851K_integrated_entropy_contribution_meV_atom']),
      'C_minus_B_integrated_entropy_contribution_at_2499K_meV_atom':float(c['at_2499K_integrated_entropy_contribution_meV_atom'])-float(b['at_2499K_integrated_entropy_contribution_meV_atom']),
      'C_minus_B_total_deltaG_error_at_851K_meV_atom':float(c['at_851K_deltaG_error_meV_atom'])-float(b['at_851K_deltaG_error_meV_atom']),
      'C_minus_B_total_deltaG_error_at_2499K_meV_atom':float(c['at_2499K_deltaG_error_meV_atom'])-float(b['at_2499K_deltaG_error_meV_atom']),
    })
with (ROOT/'ANCHOR_ENTROPY_PAIRED_B_TO_C.csv').open('w', newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(paired[0]));w.writeheader();w.writerows(paired)
for r in summary:
    print(r['representation'],r['gap_first_minus_second'],'anchor',round(r['anchor_error_at_T0_meV_atom'],6),'G MAE',round(r['deltaG_error_MAE_meV_atom'],6),'S MAE',round(r['deltaS_error_MAE_micro_eV_atom_K'],6),'closure max',r['decomposition_max_abs_closure_meV_atom'])
print('PAIRED C-B')
for r in paired: print(r)
