#!/usr/bin/env python3
"""Record the existing two-source rank limit and the C2221 admission gate."""
from pathlib import Path
import csv,json,hashlib,os
import numpy as np
RUN=Path(__file__).resolve().parents[1]
P=Path(os.environ.get('FES_PORTAL_ARCHIVE','fes-paper-evidence-portal-20261026'))
R=Path(os.environ.get('FES_TRANSFER_ARCHIVE','decisive_transfer_C_lopo_loso_20261020'))
BASE=Path(os.environ.get('FES_BENCH_ROOT','.'))
RANK_REL=1e-4

def csvrows(p):
    with p.open(newline='') as f:return list(csv.DictReader(f))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def svd_report(a):
    s=np.linalg.svd(np.asarray(a,float),compute_uv=False)
    r=int(np.sum(s>=s[0]*RANK_REL)) if len(s) and s[0]>0 else 0
    return {'shape':list(np.asarray(a).shape),'singular_values':s.tolist(),'effective_rank_threshold_relative':RANK_REL,'effective_rank':r,'numpy_matrix_rank_default':int(np.linalg.matrix_rank(a))}
def thermo_basis(t):
    x=np.asarray(t,float)/1000
    return np.column_stack([np.ones_like(x),x,x*np.log(x),x*x])
def extxyz_features(path):
    lines=path.read_text().splitlines();n=int(lines[0]);import re
    m=re.search(r'Lattice="([^"]+)"',lines[1]);cell=np.array([float(x) for x in m.group(1).split()]).reshape(3,3)
    species=[row.split()[0] for row in lines[2:2+n]]
    masses={'Hf':178.49,'Zr':91.224,'Ti':47.867,'O':15.999,'Si':28.085}
    comp={el:species.count(el)/n for el in ['Hf','Zr','Ti','O','Si']}
    return np.array([abs(np.linalg.det(cell))/n,sum(masses[e] for e in species)/n,*[comp[e] for e in ['Hf','Zr','Ti','O','Si']]])
def main():
    structures=R/'inputs/data/processed/sio2'; drows=csvrows(R/'inputs/representations/DPA_DESCRIPTOR_MEAN.csv')
    desc={(r['system'],r['phase']):np.array([float(r[f'd{i:03d}']) for i in range(128)]) for r in drows}
    phases=['quartz_beta','tridymite_p63mmc','cristobalite_beta']; pairorder=[('quartz_beta','tridymite_p63mmc'),('cristobalite_beta','quartz_beta'),('cristobalite_beta','tridymite_p63mmc')]
    rank={}
    for method in ['C-simple','C-DPA']:
        z=np.load(R/'results'/f'sio2_lopo_cristobalite_beta__{method.replace("-","_")}_readout.npz')
        if method=='C-simple': raw={p:extxyz_features(structures/p/'structure.extxyz') for p in phases}
        else: raw={p:desc[('sio2',p)] for p in phases}
        mu=np.asarray(z['feature_mean']); sd=np.asarray(z['feature_scale']); norm={p:(raw[p]-mu)/sd for p in phases}
        feat=np.vstack([norm[p] for p in phases])
        X=np.vstack([norm[a]-norm[b] for a,b in pairorder])
        W=np.asarray(z['W']); ysrc=np.asarray(z['source_pair_coefficients'])
        pred=np.vstack([(norm[a]-norm[b])@W for a,b in pairorder])
        # Diagnostic-only actual target coefficients. The complete outer labels are not read by fitting or selecting W.
        coeff={}
        for phase in phases:
            rr=csvrows(structures/phase/'reference_G.csv')
            t=np.array([float(x['T_K']) for x in rr]);g=np.array([float(x['G_eV_per_atom']) for x in rr])
            coeff[phase]=np.linalg.lstsq(thermo_basis(t),g,rcond=None)[0]
        actual=np.vstack([coeff[a]-coeff[b] for a,b in pairorder])
        rank[method]={
          'feature_name':'standardized 7D volume/mass/composition' if method=='C-simple' else 'standardized frozen DPA descriptor mean (128D)',
          'independent_source_phase_contrasts':1,
          'feature_contrast_design_all_three_phases':svd_report(X),
          'source_supervised_feature_design':svd_report(np.asarray(z['source_pair_features'])),
          'source_label_coefficient_contrast':{'matrix':ysrc.tolist(),'svd':svd_report(ysrc)},
          'W':{'svd':svd_report(W),'shape':list(W.shape),'exact_linear_ridge_rank_ceiling_from_one_source_contrast':1},
          'predicted_all_pair_coefficient_contrasts':{'pair_order':[list(x) for x in pairorder],'matrix':pred.tolist(),'svd':svd_report(pred)},
          'reference_all_pair_coefficient_contrasts_diagnostic_only':{'temperature_basis':'[1,x,x ln x,x²], x=T/1000 K','fit_support':'full 1649-point archived reference curves; outer labels consumed only for this diagnostic after the archived readout was frozen','matrix':actual.tolist(),'svd':svd_report(actual)},
          'source_phase_feature_names':'C-simple: [volume/atom, mean atomic mass, fractions Hf,Zr,Ti,O,Si]; C-DPA: descriptor columns d000..d127'}
    # Quantify the source table pressure normalization for every exact common-grid pair.
    root=BASE/'fes-bench/data/processed/sio2_source_v2/phases';v={}
    for phase in ['quartz_beta','cristobalite_beta','tridymite_p63mmc','tridymite_c2221']:
        rr=csvrows(root/phase/'thermodynamic_properties.csv');v[phase]={float(x['T_K']):float(x['V_A3_per_atom']) for x in rr}
    names=list(v);common=sorted(set.intersection(*[set(x) for x in v.values()]));pressure=1e-4;ev=1.602176634e-19;pressure_rows=[]
    for i,a in enumerate(names):
      for b in names[i+1:]:
        dmev=[pressure*1e9*(v[a][t]-v[b][t])*1e-30/ev*1000 for t in common]
        pressure_rows.append({'phase_a':a,'phase_b':b,'n_common_states':len(common),'pressure_GPa':pressure,'max_abs_delta_G_for_GPa_rounding_meV_atom':max(abs(x) for x in dmev),'min_delta_meV_atom':min(dmev),'max_delta_meV_atom':max(dmev)})
    admission={
      'experiment':'S1 independent source contrast admission',
      'source_G_table_status':'ADMITTED for same-source thermodynamic G and phase-pair baseline evaluation at exact common state nodes; uncertainty remains unknown, not zero.',
      'structure_migration_status':'NOT ADMITTED as a C2221 structure-to-G training pair: source-row geometry is not mapped to a unique structure.',
      'evidence':{
        'raw_table_reproduction':'1649 common G rows reproduced exactly from publisher-verified source archive; table phase identity/historical source level support thermodynamic-only comparisons.',
        'FCP':'A 12-atom O8Si4 C2221 PBE-D3(BJ) FCP reference exists; no evidence ties every thermodynamic row to this exact reference geometry.',
        'COD_candidate':'COD 8104538 conventional cell is not identical to the FCP-derived cell: aligned axis mismatch up to 4.53%; mapped-site RMS 0.384 Å; max displacement 0.558 Å. Do not attach the original G labels to COD geometry.',
        'pressure':'Published raw P is approximately 0.0001 GPa (0.1 MPa); processed P coordinate is rounded to 0. Quantified pair PΔV shifts are in S1_PRESSURE_ROUNDING_IMPACT.csv.',
        'uncertainty':'Per-row uncertainty/covariance are not supplied; processed size_error=0 is a placeholder.',
        'vibration':'Corrected C2221 harmonic free energy is a fixed-cell harmonic candidate, not complete QHA; source output had a factor-24 per-cell/per-atom error and includes unstable Gamma modes.'},
      'minimum_missing_for_source_migration':'Original source geometry associated with the table/FCP, or an explicitly separate newly computed reference layer on an explicitly selected structure and compatible protocol. The COD candidate cannot inherit the existing G table.',
      'decision':'S1-B/D additional-source contrast experiment blocked; proceed on the original three-phase S1-A/C baseline and S2/S3.'}
    result={'protocol':{'gauge':'shared quartz_beta at each T','basis':'[1,x,x ln x,x²], x=T/1000 K','feature_standardization':'source structures only, archived means/scales','effective_rank_threshold_relative':RANK_REL,'exact_ridge_rank_argument':'With one independent source feature contrast and positive ridge, the exact solution W has rank at most one; preserve observed numerical singular values.'},'methods':rank,'c2221_admission':admission,'source_pressure_pair_effects':pressure_rows,'provenance':{'C_DPA_npz_sha256':sha(R/'results/sio2_lopo_cristobalite_beta__C_DPA_readout.npz'),'C_simple_npz_sha256':sha(R/'results/sio2_lopo_cristobalite_beta__C_simple_readout.npz'),'descriptor_sha256':sha(R/'inputs/representations/DPA_DESCRIPTOR_MEAN.csv'),'phase_lopo_archive_sha256':sha(R/'fold_manifest.json'),'C2221_review_sha256':sha(P/'assets/sources/9c5acf9b435f02_C2221_SOURCE_STRUCTURE_REVIEW_20261023.md'),'C2221_contract_sha256':sha(P/'assets/sources/4e82bfe5bec2b0_C2221_SOURCE_CONTRACT_AUDIT_20261024.json'),'C2221_reaudit_sha256':sha(P/'assets/sources/eacb903402801e_C2221_STRUCTURE_AND_QH_REAUDIT.md'),'corrected_vib_sha256':sha(P/'assets/sources/b23523aadd3475_C2221_CORRECTED_VIBRATIONAL_FREE_ENERGY.csv')}}
    (RUN/'data/S1_RANK_SVD_AND_C2221_ADMISSION.json').write_text(json.dumps(result,indent=2)+'\n')
    with (RUN/'data/S1_PRESSURE_ROUNDING_IMPACT.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(pressure_rows[0]));w.writeheader();w.writerows(pressure_rows)
    print(json.dumps({'methods':{m:{'source_feature_sv':v['source_supervised_feature_design']['singular_values'],'W_sv':v['W']['svd']['singular_values'],'predicted_pair_sv':v['predicted_all_pair_coefficient_contrasts']['svd']['singular_values'],'reference_pair_sv':v['reference_all_pair_coefficient_contrasts_diagnostic_only']['svd']['singular_values']} for m,v in rank.items()},'pressure_pair_effects':pressure_rows,'decision':admission['decision']},indent=2))
if __name__=='__main__':main()
