#!/usr/bin/env python3
"""Render frozen S2 information-ladder figures using the saved Python backend."""
from pathlib import Path
import argparse, csv, sys, os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans']
plt.rcParams['svg.fonttype'] = 'none'
# Keep editable vector text: svg.fonttype = "none"; pdf.fonttype = 42.
plt.rcParams.update({'pdf.fonttype':42,'font.size':7,'axes.spines.right':False,'axes.spines.top':False,
                     'axes.linewidth':0.7,'legend.frameon':False,'xtick.major.width':0.7,'ytick.major.width':0.7})
COLORS={'zero_shot':'#4C78A8','oracle_gap':'#F2A541','oracle_gap_slope':'#54A24B','few_shot_1_label':'#B279A2','few_shot_2_label':'#E45756'}
LABELS={'zero_shot':'Zero-shot (0)','oracle_gap':'Oracle gap (1 label)','oracle_gap_slope':'Oracle gap+slope (11 labels)','few_shot_1_label':'Few-shot (1 label)','few_shot_2_label':'Few-shot (2 labels)'}
METHODS=['C-simple','C-DPA']; ARMS=list(COLORS)
SHORT_LABELS={'zero_shot':'Zero-shot','oracle_gap':'Gap\noracle','oracle_gap_slope':'Gap + slope\noracle','few_shot_1_label':'Few-shot\n1 label','few_shot_2_label':'Few-shot\n2 labels'}

def readcsv(path):
    with path.open(newline='') as f:return list(csv.DictReader(f))
def save(fig,path):
    fig.canvas.draw()
    fig.savefig(str(path)+'.svg',bbox_inches='tight')
    fig.savefig(str(path)+'.pdf',bbox_inches='tight')
    fig.savefig(str(path)+'.png',dpi=600,bbox_inches='tight')
    fig.savefig(str(path)+'.tiff',dpi=600,bbox_inches='tight')
    plt.close(fig)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-dir',type=Path,required=True);out=ap.parse_args().run_dir; figdir=out/'figures'; datadir=out/'data'
    summary=readcsv(datadir/'S2_INFORMATION_LADDER_SUMMARY.csv')
    # Figure 1: paired information ladder. Same samples and information ordering in both panels.
    fig,axs=plt.subplots(1,2,figsize=(7.1,2.75),sharex=True)
    x=np.arange(len(ARMS)); offsets={'C-simple':-0.12,'C-DPA':0.12}; mcols={'C-simple':'#285F8F','C-DPA':'#9B4E68'}
    for ax,key,ylabel,title in [(axs[0],'target_quartz_gap_MAE_meV_atom','|ΔG error| (meV atom⁻¹)','Relative gap'),(axs[1],'target_quartz_slope_MAE_meV_atom_K','|relative slope error| (meV atom⁻¹ K⁻¹)','Relative slope')]:
        for method in METHODS:
            rr=[next(r for r in summary if r['method']==method and r['arm']==arm) for arm in ARMS]
            vals=np.array([float(r[key]) for r in rr]); xx=x+offsets[method]
            ax.plot(xx,vals,'o-',color=mcols[method],lw=1.3,ms=4,label=method,zorder=3)
        ax.set_xticks(x, [SHORT_LABELS[a] for a in ARMS],rotation=0,rotation_mode='anchor')
        ax.set_ylabel(ylabel);ax.set_title(title,loc='left',fontweight='bold',fontsize=8);ax.grid(axis='y',color='#dddddd',lw=.5)
    axs[0].legend(loc='upper right',fontsize=6.5)
    fig.suptitle('S2 | Gap and slope information injected into held-out cristobalite',x=.06,ha='left',fontsize=9,fontweight='bold')
    fig.text(.06,.005,'C-simple / C-DPA; LOPO cristobalite_beta; common score support excludes 1263 and 2087 K. Oracle slope uses 11 target temperatures.',fontsize=6.5)
    fig.tight_layout(rect=[0,.10,1,.9],w_pad=1.5);save(fig,figdir/'S2_INFORMATION_LADDER')
    # Figure 2: all phase-pair curves. Visual curves use the full support; formal scores exclude both few-shot coordinates.
    pts=readcsv(datadir/'S2_POINT_PREDICTIONS.csv')
    T=np.array(sorted({float(r['T_K']) for r in pts if r['method']=='C-simple'}))
    pair_names=[('cristobalite_beta','quartz_beta','Cristobalite − quartz'),('cristobalite_beta','tridymite_p63mmc','Cristobalite − tridymite'),('quartz_beta','tridymite_p63mmc','Quartz − tridymite')]
    for method in METHODS:
        fig,axs=plt.subplots(1,3,figsize=(7.25,3.0),sharex=True)
        for ax,(a,b,title) in zip(axs,pair_names):
            ref=[]
            for arm in ARMS:
                rows=[r for r in pts if r['method']==method and r['arm']==arm and r['phase']=='cristobalite_beta']
                byT={float(r['T_K']):r for r in rows}
                # Quartz is the zero gauge. Tridymite is frozen across all arms.
                if a=='cristobalite_beta': va=np.array([float(byT[tt]['target_minus_quartz_corrected_prediction_eV_atom']) for tt in T])
                else: va=np.zeros(len(T))
                if b=='quartz_beta': vb=np.zeros(len(T))
                elif b=='tridymite_p63mmc':
                    base=__import__('csv').DictReader(open(Path(os.environ.get('FES_PORTAL_ARCHIVE','fes-paper-evidence-portal-20261026'))/'figures/source_data/sio2_envelope_audit.csv'))
                    base_map={(r['method'],float(r['T_K'])):float(r['predicted_tridymite_minus_quartz_eV_atom']) for r in base}
                    vb=np.array([base_map[(method,tt)] for tt in T])
                else: vb=np.zeros(len(T))
                ax.plot(T,(va-vb)*1000,color=COLORS[arm],lw=1.1,label=LABELS[arm])
            base=readcsv(Path(os.environ.get('FES_PORTAL_ARCHIVE','fes-paper-evidence-portal-20261026'))/'figures/source_data/sio2_envelope_audit.csv')
            base_map={(r['method'],float(r['T_K'])):r for r in base}
            if a=='cristobalite_beta': va=np.array([float(base_map[(method,tt)]['reference_cristobalite_minus_quartz_eV_atom']) for tt in T])
            else: va=np.zeros(len(T))
            if b=='quartz_beta': vb=np.zeros(len(T))
            elif b=='tridymite_p63mmc': vb=np.array([float(base_map[(method,tt)]['reference_tridymite_minus_quartz_eV_atom']) for tt in T])
            else: vb=np.zeros(len(T))
            ax.plot(T,(va-vb)*1000,color='#222222',lw=1.6,ls='--',label='Reference')
            ax.axhline(0,color='#777777',lw=.6)
            ax.set_title(title,loc='left',fontsize=7.5,fontweight='bold');ax.set_xlabel('Temperature (K)');ax.grid(axis='y',color='#e3e3e3',lw=.4)
            ax.set_xlim(T.min(),T.max())
        axs[0].set_ylabel('ΔG (meV atom⁻¹)')
        handles,labels=axs[0].get_legend_handles_labels();fig.legend(handles,labels,ncol=3,loc='lower center',bbox_to_anchor=(.5,.075),fontsize=6.2,frameon=False)
        fig.suptitle(f'S2 | Pair curves after phase-level information injection — {method}',x=.06,ha='left',fontsize=9,fontweight='bold')
        fig.text(.06,.035,'Full curves shown for interpretation; few-shot coordinates are labeled inputs. Formal scores use common support with both coordinates removed.',fontsize=6.2)
        fig.tight_layout(rect=[0,.23,1,.88],w_pad=1.1);save(fig,figdir/f'S2_PAIR_CURVES_{method.replace("-","_")}')
    # Figure 3: event position and envelope decision metrics, preserving the distinct metric denominators.
    env=readcsv(datadir/'S2_ENVELOPE_METRICS.csv')
    fig,axs=plt.subplots(1,2,figsize=(7.1,2.8),sharex=True)
    x=np.arange(len(ARMS))
    for method,off,col in [('C-simple',-.12,'#285F8F'),('C-DPA',.12,'#9B4E68')]:
        rr=[next(r for r in env if r['method']==method and r['arm']==arm) for arm in ARMS]
        refTc=float(rr[0]['reference_roots_K'])
        tc=[float(r['predicted_roots_K']) if r['predicted_roots_K'] else np.nan for r in rr]
        axs[0].plot(x+off,tc,'o-',color=col,lw=1.2,ms=3.5,label=method)
        axs[1].plot(x+off,[float(r['stable_envelope_accuracy'])*100 for r in rr],'o-',color=col,lw=1.2,ms=3.5,label=method+' accuracy')
        axs[1].plot(x+off,[float(r['target_stable_recall'])*100 for r in rr],'s--',color=col,lw=.9,ms=3,label=method+' target recall')
    axs[0].axhline(refTc,color='#222222',ls='--',lw=1,label=f'Reference {refTc:.2f} K')
    axs[0].set_ylabel('Cristobalite boundary (K)');axs[1].set_ylabel('Correct decisions / target recall (%)');axs[1].set_ylim(90,101)
    axs[0].set_title('Boundary position',loc='left',fontsize=8,fontweight='bold');axs[1].set_title('Envelope and target-phase metrics',loc='left',fontsize=8,fontweight='bold')
    for ax in axs:
        ax.set_xticks(x,[SHORT_LABELS[a] for a in ARMS],rotation=0,rotation_mode='anchor');ax.grid(axis='y',color='#ddd',lw=.5)
    axs[0].legend(fontsize=6.2,loc='best');axs[1].legend(fontsize=6.2,loc='lower right')
    fig.suptitle('S2 | Boundary recovery and stable-envelope decisions',x=.06,ha='left',fontsize=9,fontweight='bold')
    fig.text(.06,.005,'LOPO cristobalite_beta; V2 tie rule; common support excludes 1263 and 2087 K; event match tolerance 100 K. Root error is also reported explicitly in CSV.',fontsize=6.2)
    fig.tight_layout(rect=[0,.12,1,.9],w_pad=1.4);save(fig,figdir/'S2_ENVELOPE_AND_BOUNDARIES')
    print('Rendered S2 PDF/SVG/PNG figures and alignment JSON.')
if __name__=='__main__':main()
