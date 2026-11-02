"""Frozen paired SiO2 A/B/C readouts; all target roles are development cases."""
import csv
import hashlib
import itertools
import json
from pathlib import Path
import platform
import numpy as np

OUT = Path(__file__).resolve().parent
OLD = Path('project/outputs/fes-sio2-supervision-physics-goal-20261028')
PD = Path('project/fes-bench/data/processed/sio2_source_v2/phases')
FOLD = Path('project/outputs/fes-bench-cross-system-20261020/results/decisive_transfer_C_lopo_loso_20261020/fold_manifest.json')
PH = ['quartz_beta', 'tridymite_p63mmc', 'tridymite_c2221', 'cristobalite_beta']
V = dict(zip(PH, [13.459893903879466, 15.404241970876686, 16., 15.35390250733693]))
GRID = [1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,.01,.1,1,10,100,1000,10000]

def table(path, rows):
    if not rows:
        return
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

def basis(t):
    x = np.asarray(t)/1000
    return np.c_[np.ones(len(x)), x, x*np.log(x), x*x]

def derivative(t):
    x = np.asarray(t)/1000
    return np.c_[np.zeros(len(x)), np.ones(len(x)), 1+np.log(x), 2*x]/1000

def solve(X, Y, H, alpha):
    # Exact reduced-rank solution of ||X W-Y||_H^2/n + alpha ||W||_F^2.
    # Same H for all source pairs because their temperature support is identical.
    lam, Q = np.linalg.eigh(H)
    U, s, VT = np.linalg.svd(X, full_matrices=False)
    weights = s[:, None]*lam[None, :]/(s[:, None]**2*lam[None, :]+len(X)*alpha)
    return VT.T @ (weights * (U.T @ Y @ Q)) @ Q.T

def events(t, pool, curves):
    lab = np.array(pool)[np.argmin(np.column_stack([curves[p] for p in pool]), axis=1)]
    out = []
    for i in np.flatnonzero(lab[:-1] != lab[1:]):
        a,b = str(lab[i]),str(lab[i+1]); d0=curves[a][i]-curves[b][i]; d1=curves[a][i+1]-curves[b][i+1]
        f = abs(d0)/(abs(d0)+abs(d1)) if d0*d1<0 else .5
        out.append({'from':a,'to':b,'T_K':float(t[i]+f*(t[i+1]-t[i]))})
    return lab, out

def match(ref, pred):
    n,m=len(ref),len(pred); dp=[[(0,0.) for _ in range(m+1)] for _ in range(n+1)]
    def better(a,b): return a[0]>b[0] or (a[0]==b[0] and a[1]<b[1])
    for i in range(n+1):
        for j in range(m+1):
            v=dp[i][j]
            for ni,nj in [(i+1,j),(i,j+1)]:
                if ni<=n and nj<=m and better(v,dp[ni][nj]): dp[ni][nj]=v
            if i<n and j<m and ref[i]['from']==pred[j]['from'] and ref[i]['to']==pred[j]['to']:
                d=abs(ref[i]['T_K']-pred[j]['T_K']); nv=(v[0]+1,v[1]+d)
                if d<=100 and better(nv,dp[i+1][j+1]): dp[i+1][j+1]=nv
    return dp[n][m]

def features(emb, src, arm, scaling):
    d=np.vstack([emb[p+'_mean_descriptor'] for p in PH]); v=np.array([V[p] for p in PH])[:,None]
    ids=[PH.index(p) for p in src]; blocks=[]; means=[]; scales=[]
    for raw in ([v] if arm=='A' else [d] if arm=='B' else [d,v]):
        mu=raw[ids].mean(0); sd=raw[ids].std(0); sd[(sd<1e-12)|~np.isfinite(sd)]=1.
        group=np.sqrt(raw.shape[1]) if scaling=='equal_block' else 1.
        blocks.append((raw-mu)/sd/group); means.append(mu); scales.append(sd*group)
    z=np.column_stack(blocks)
    return dict(zip(PH,z)), np.concatenate(means), np.concatenate(scales)

def read(p):
    rows=list(csv.DictReader(open(PD/p/'reference_G.csv')))
    return np.array([float(r['T_K']) for r in rows]),np.array([float(r['G_eV_per_atom']) for r in rows])

def main():
    emb=np.load(OLD/'inputs/structures/task_checkpoint_embeddings.npz')
    fold=json.load(open(FOLD))['folds']['sio2_lopo_cristobalite_beta']
    fit=sorted({int(x['T_index']) for x in fold['fit'] if x['phase']=='quartz_beta'})
    val=sorted({int(x['T_index']) for x in fold['source_validation'] if x['phase']=='quartz_beta'})
    models=[]; source_rows=[]; selection=[]; hashes=[]
    inputs=[OUT/'PROTOCOL.json',Path(__file__),FOLD,OLD/'inputs/structures/task_checkpoint_embeddings.npz']+[PD/p/'reference_G.csv' for p in PH]
    for p in inputs:
        hashes.append({'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
    table(OUT/'INPUT_HASHES.csv',hashes)
    (OUT/'models').mkdir(exist_ok=True)
    for target in PH:
        sources=[p for p in PH if p!=target]
        configs=[('low_'+'+'.join(s),list(s)) for s in itertools.combinations(sources,2)]+[('high_all_three',sources)]
        for name,src in configs:
            # Only source labels enter coefficient fitting, feature scaling and selection.
            data={p:read(p) for p in src}; t=data[src[0]][0]
            assert all(np.array_equal(t,data[p][0]) for p in src)
            b=basis(t); H=b[fit].T@b[fit]/len(fit)
            coef={p:np.linalg.lstsq(b[fit],data[p][1][fit],rcond=None)[0] for p in src}
            pairs=list(itertools.combinations(src,2)); Y=np.array([coef[a]-coef[c] for a,c in pairs])
            for scaling in ['coordinate','equal_block']:
                for arm in ['A','B','C']:
                    z,mu,sd=features(emb,src,arm,scaling)
                    X=np.array([z[a]-z[c] for a,c in pairs]); scale=max(np.trace(X.T@X)/X.shape[1],1e-15)
                    options=[]
                    for mult in GRID:
                        w=solve(X,Y,H,mult*scale)
                        errors=np.concatenate([b[val]@((z[a]-z[c])@w)-(data[a][1][val]-data[c][1][val]) for a,c in pairs])
                        options.append((float(np.mean(abs(errors))),mult,w))
                    score,mult,w=min(options,key=lambda x:(x[0],x[1]))
                    key=f'{target}__{name}__{scaling}__{arm}'
                    np.savez(OUT/'models'/f'{key}.npz',W=w,feature_mean=mu,feature_scale=sd,source_pair_features=X,source_pair_coefficients=Y,H=H,alpha=mult*scale)
                    common={'target':target,'configuration':name,'sources':';'.join(src),'scaling':scaling,'arm':arm}
                    selection.append(dict(common,source_validation_MAE_eV_atom=score,alpha_multiplier=mult,alpha=mult*scale,source_rank=int(np.linalg.matrix_rank(X)),n_source_fit_labels=len(src)*len(fit)))
                    for a,c in pairs:
                        dc=(z[a]-z[c])@w-(coef[a]-coef[c]); errs=b[fit]@((z[a]-z[c])@w)-(data[a][1][fit]-data[c][1][fit])
                        source_rows.append(dict(common,pair=a+'-'+c,source_pair_fit_RMSE_meV_atom=float(np.sqrt(np.mean(errs**2))*1000),**{f'coefficient_error_{i}_eV_atom':float(x) for i,x in enumerate(dc)}))
                    models.append((common,src,z,w,coef))
    table(OUT/'FROZEN_SELECTION.csv',selection); table(OUT/'SOURCE_FIT.csv',source_rows)
    selection_hash=hashlib.sha256((OUT/'FROZEN_SELECTION.csv').read_bytes()).hexdigest()
    (OUT/'SELECTION_FREEZE.json').write_text(json.dumps({'sha256':selection_hash,'target_scoring_started':False},indent=2)+'\n')
    # All model selections are serialized above. Labels below are evaluation-only.
    data={p:read(p) for p in PH}; t=data[PH[0]][0]; b=basis(t); db=derivative(t); pairrows=[]; envrows=[]; mechanism=[]; eventrows=[]
    for common,src,z,w,coef in models:
        target=common['target']; cp={p:z[p]@w for p in PH}; pred={p:b@cp[p] for p in PH}; truth={p:data[p][1] for p in PH}
        # Affine source projection with the sum-to-one constraint enforced explicitly.
        S=np.vstack([z[p] for p in src]); base=S[-1]; D=(S[:-1]-base).T
        weights=np.linalg.lstsq(D,z[target]-base,rcond=None)[0]; weights=np.r_[weights,1-weights.sum()]
        mechanism.append(dict(common,source_weights=json.dumps(weights.tolist()),target_affine_residual=float(np.linalg.norm(z[target]-weights@S)),prediction_combination_residual=float(np.max(abs(cp[target]-weights@np.vstack([cp[p] for p in src]))))))
        for anchor in src:
            dc=cp[target]-cp[anchor]; gap=b@dc; ref=truth[target]-truth[anchor]; err=gap-ref
            local=np.argsort(abs(t-1000))[:21]; rc=np.linalg.lstsq(b[local],ref[local],rcond=None)[0]
            entropy=-(db@dc); refentropy=-np.gradient(ref,t,edge_order=2)
            pairrows.append(dict(common,anchor=anchor,gap_MAE_meV_atom=float(np.mean(abs(err))*1000),gap_bias_meV_atom=float(np.mean(err)*1000),offset_error_1000K_meV_atom=float(np.interp(1000,t,err)*1000),slope_error_1000K_microeV_atom_K=float(abs(derivative([1000])[0]@(dc-rc))*1e6),relative_entropy_MAE_microeV_atom_K=float(np.mean(abs(entropy-refentropy))*1e6),relative_entropy_bias_microeV_atom_K=float(np.mean(entropy-refentropy)*1e6)))
        subsets=[src] if len(src)==2 else list(itertools.combinations(src,2))
        for subset in subsets+[tuple(p for p in PH if p!=target)]:
            pool=list(subset)+[target]
            if len(src)==2 and len(pool)==3 and list(subset)==src and any(r['target']==target and r['configuration']==common['configuration'] and r['scaling']==common['scaling'] and r['arm']==common['arm'] and r['candidate_pool']==';'.join(pool) for r in envrows): continue
            rl,re=events(t,pool,truth); pl,pe=events(t,pool,pred); hits,te=match(re,pe)
            envrows.append(dict(common,candidate_pool=';'.join(pool),stable_accuracy=float(np.mean(rl==pl)),reference_events=len(re),predicted_events=len(pe),hits=hits,missed=len(re)-hits,false=len(pe)-hits,Tc_MAE_K=te/hits if hits else ''))
            eventrows.append(dict(common,candidate_pool=';'.join(pool),reference=json.dumps(re),predicted=json.dumps(pe)))
    table(OUT/'PAIR_METRICS.csv',pairrows); table(OUT/'ENVELOPE_METRICS.csv',envrows); table(OUT/'EVENTS.csv',eventrows); table(OUT/'SOURCE_RESPONSE_WEIGHTS.csv',mechanism)
    increments=[]
    metrics=['gap_MAE_meV_atom','relative_entropy_MAE_microeV_atom_K','slope_error_1000K_microeV_atom_K']
    for row in pairrows:
        if row['arm']!='C':continue
        a=next(x for x in pairrows if all(x[k]==row[k] for k in ['target','configuration','scaling','anchor']) and x['arm']=='A')
        d=next(x for x in pairrows if all(x[k]==row[k] for k in ['target','configuration','scaling','anchor']) and x['arm']=='B')
        common={k:row[k] for k in ['target','configuration','sources','scaling','anchor']}
        increments.append(dict(common,**{f'DPA_given_V__{k}':a[k]-row[k] for k in metrics},**{f'V_given_DPA__{k}':d[k]-row[k] for k in metrics}))
    table(OUT/'PAIRED_INCREMENTS.csv',increments)
    qa={'n_models':len(models),'selection_sha256':selection_hash,'python':platform.python_version(),'numpy':np.__version__,'source_only_selection':True,'max_affine_prediction_closure':max(x['prediction_combination_residual'] for x in mechanism)}
    # Check archived primary DPA reconstruction. This is a numerical regression,
    # not a new fit or an opportunity to select models on target performance.
    oldrows=list(csv.DictReader(open(OLD/'CDPA_TASKCKPT_CRISTOBALITE_PAIR_SCORES.csv')))
    diffs=[]
    for row in pairrows:
        if row['target']=='cristobalite_beta' and row['arm']=='B' and row['scaling']=='coordinate':
            old=next(x for x in oldrows if x['configuration']==row['configuration'] and x['anchor']==row['anchor'])
            diffs.append(abs(row['gap_MAE_meV_atom']-float(old['gap_MAE'])*1000))
    qa['max_archived_primary_DPA_gap_difference_meV_atom']=max(diffs)
    assert max(diffs)<1e-5, qa
    (OUT/'QA.json').write_text(json.dumps(qa,indent=2)+'\n')
    print(json.dumps(qa,indent=2))

if __name__=='__main__': main()
