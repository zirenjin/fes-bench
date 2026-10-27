#!/usr/bin/env python3
"""Show source-only and full envelopes on one common phase gauge and T grid."""
import json,csv,numpy as np
import plot_evidence as p
R=p.ROOT;p.SELECT.clear();plt=p.plt
with (R/'figures/source_data/sio2_envelope_audit.csv').open() as f:rows=list(csv.DictReader(f))
for zoom in [False,True]:
 fig,axes=plt.subplots(1,3,figsize=(p.W,4.15),sharey=True);fig.subplots_adjust(left=.11,right=.98,bottom=.21,top=.75,wspace=.22)
 for ax,title,method in zip(axes,['Reference','C-simple','C-DPA'],['reference','C-simple','C-DPA']):
  rr=[r for r in rows if r['method']==('C-simple' if method=='reference' else method)];t=np.array([float(r['T_K']) for r in rr]);prefix='reference' if method=='reference' else 'predicted';tr=np.array([float(r[prefix+'_tridymite_minus_quartz_eV_atom'])*1000 for r in rr]);cr=np.array([float(r[prefix+'_cristobalite_minus_quartz_eV_atom'])*1000 for r in rr]);source=np.minimum(0,tr);full=np.minimum(source,cr)
  ax.axhline(0,color=p.C['QH'],lw=.8,label='Quartz gauge');ax.plot(t,tr,color=p.C['T3'],lw=1.3,label='Tridymite − quartz');ax.plot(t,cr,color=p.C['T2'],lw=1.3,label='Cristobalite − quartz');ax.plot(t,source,color='#8A939B',lw=2.3,ls='--',alpha=.8,label='Source-only envelope');ax.plot(t,full,color=p.C['ink'],lw=1.0,label='Full envelope');ax.set(xlabel='Temperature (K)',xticks=[1000,1500,2000,2500]);p.clean(ax,False);p.panel(ax,'abc'[list(axes).index(ax)],title)
 axes[0].set_ylabel('Free energy relative to quartz (meV/atom)')
 if zoom:
  vals=[]
  for rr in rows:
   if rr['method'] not in ['C-simple','C-DPA'] or not 1450<=float(rr['T_K'])<=1620:continue
   vals.extend(float(rr[k])*1000 for k in ['reference_tridymite_minus_quartz_eV_atom','reference_cristobalite_minus_quartz_eV_atom','predicted_tridymite_minus_quartz_eV_atom','predicted_cristobalite_minus_quartz_eV_atom'])
  lo,hi=min(vals),max(vals);pad=(hi-lo)*.12
  for ax in axes:ax.set(xlim=(1450,1620),ylim=(lo-pad,hi+pad),xticks=[1450,1500,1550,1600])
 fig.legend(*axes[0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.53,.99),ncol=2,frameon=False,fontsize=7)
 p.export('sio2_envelope_audit_zoom' if zoom else 'sio2_envelope_audit','SiO₂ LOPO: crossing-region detail' if zoom else 'SiO₂ LOPO: target participation and coupled pair crossings',fig,rows,('View: 1450–1620 K. ' if zoom else 'View: full temperature support. ')+'All 1,649 outer-test states per readout. Curves use quartz as a common additive gauge. Dashed grey: the two source phases only. Black: the envelope after adding cristobalite. No phase-specific alignment or new fit is used.','transfer-audit')
path=R/'figures/explainer_figures.json';x=json.loads(path.read_text());x['figures']=[f for f in x['figures'] if not f['key'].startswith('sio2_envelope_audit')]+p.FIGS;path.write_text(json.dumps(x,ensure_ascii=False,indent=2));(R/'figures/explainer_figures.js').write_text('window.FES_EXPLAINER_FIGURES='+json.dumps(x,ensure_ascii=False)+';\n')
