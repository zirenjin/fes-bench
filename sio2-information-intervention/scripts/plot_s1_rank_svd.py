#!/usr/bin/env python3
"""Plot archived S1 singular spectra without refitting or using labels for selection."""
from pathlib import Path
import argparse, json, sys, os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
SKILL_DIR=os.environ.get('NATURE_FIGURE_SKILL')
if SKILL_DIR: sys.path.insert(0,str(Path(SKILL_DIR)/'scripts'))
try:
 from audit_panel_alignment import require_matplotlib_panel_alignment
except ImportError:
 require_matplotlib_panel_alignment=None

COLORS={'C-simple':'#285F8F','C-DPA':'#9B4E68'}
THRESH=1e-4
POSITIVITY_EPSILON=1e-30

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--run-dir',required=True);a=ap.parse_args()
 root=Path(a.run_dir);data=json.loads((root/'data/S1_RANK_SVD_AND_C2221_ADMISSION.json').read_text())['methods']
 mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans','Liberation Sans'],'font.size':7,'pdf.fonttype':42,'svg.fonttype':'none','axes.linewidth':.7,'axes.spines.right':False,'axes.spines.top':False,'legend.frameon':False})
 panels=[('W','W','Fitted readout matrix W','1e−30'),('predicted_all_pair_coefficient_contrasts','predicted_all_pair_coefficient_contrasts','Predicted pair coefficient contrasts','1e−22'),('reference_all_pair_coefficient_contrasts_diagnostic_only','reference_all_pair_coefficient_contrasts_diagnostic_only','Reference pair coefficients · diagnostic only','1e−22')]
 fig,axs=plt.subplots(1,3,figsize=(7.25,2.75))
 spectra={}
 for ax,(key,_,title,ymin) in zip(axs,panels):
  for method in ['C-simple','C-DPA']:
   spec=data[method][key]['svd'];vals=np.asarray(spec['singular_values'],float);floor=float(ymin.replace('−','-'))
   # Log plotting requires a strictly positive floor; exact zeros stay exact in the JSON.
   plotvals=np.where(vals>0,vals,max(floor,POSITIVITY_EPSILON))
   assert np.all(plotvals>0), 'log-axis inputs must be strictly positive'
   x=np.arange(1,len(vals)+1)
   ax.plot(x,plotvals,'o-',lw=1.1,ms=3.5,color=COLORS[method],label=method)
   for xi,vi,raw in zip(x,plotvals,vals):
    if raw == 0: ax.annotate('0', (xi,vi), xytext=(0,10), textcoords='offset points',ha='center',fontsize=5.5,color=COLORS[method])
   cutoff=(vals[0]*THRESH) if len(vals) and vals[0]>0 else floor
   ax.axhline(cutoff,color=COLORS[method],ls=':',lw=.7,alpha=.8)
   spectra.setdefault(method,{})[key]={'singular_values':vals.tolist(),'effective_rank':spec['effective_rank'],'cutoff':cutoff,'numpy_matrix_rank_default':spec['numpy_matrix_rank_default']}
  ax.tick_params(axis='both',labelsize=7.5)
  ax.set_yscale('log');ax.set_ylim(float(ymin.replace('−','-')),max(float(np.max([data[m][key]['svd']['singular_values'][0] for m in ['C-simple','C-DPA'] if data[m][key]['svd']['singular_values']]))*4,1e-2))
  maxdim=max(len(data[m][key]['svd']['singular_values']) for m in ['C-simple','C-DPA'])
  ax.set_xticks(np.arange(1,maxdim+1));ax.set_xlabel('Singular value index');ax.set_title(title,loc='left',fontsize=7.3,fontweight='bold');ax.grid(axis='y',which='major',color='#e3e3e3',lw=.45)
  ax.text(.02,.02,'Dotted lines: relative cutoff = 10⁻⁴ × σ₁',transform=ax.transAxes,fontsize=5.5,color='#555')
 axs[0].set_ylabel('Singular value (log scale)')
 axs[0].legend(loc='upper right',fontsize=6)
 fig.suptitle('S1 | One supervised source contrast constrains saved C readouts',x=.06,ha='left',fontsize=9,fontweight='bold')
 fig.text(.06,.008,'C-simple / C-DPA · cristobalite_beta LOPO. Reference coefficients are post hoc diagnostic fits to held-out curves, not training targets. Exact zeros are shown at the axis floor.',fontsize=6.2)
 fig.tight_layout(rect=[0,.1,1,.88],w_pad=1.2)
 out=root/'figures'/'S1_RANK_SVD';out.parent.mkdir(parents=True,exist_ok=True)
 fig.canvas.draw()
 if require_matplotlib_panel_alignment:
  require_matplotlib_panel_alignment(fig,json_out=str(out)+'.alignment.json',overlay_svg=str(out)+'.alignment.svg',tolerance_pt=1.5,gutter_tolerance_pt=1.5,strict=True)
 fig.savefig(str(out)+'.pdf',bbox_inches='tight');fig.savefig(str(out)+'.svg',bbox_inches='tight');fig.savefig(str(out)+'.png',dpi=600,bbox_inches='tight');fig.savefig(str(out)+'.tiff',dpi=600,bbox_inches='tight');plt.close(fig)
 report={'figure':str(out),'threshold_relative':THRESH,'spectra':spectra,'scope':'existing saved readouts only; no retraining or target selection','rank_summary':{m:{k:v['effective_rank'] for k,v in spectra[m].items()} for m in spectra}}
 (root/'data/S1_RANK_SVD_FIGURE_DATA.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps(report,indent=2))
if __name__=='__main__':main()
