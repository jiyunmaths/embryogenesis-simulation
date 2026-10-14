"""Plot fixed-time qualification errors from unchanged completed reports."""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path('output/studies/repeated-development-gpu-qualification-v2')
OUT=Path('docs/figures')
BRANCHES=('off','tension_adhesion','polarity','full')
LABELS=('Coupling off','Tension / adhesion','Polarity','Full')
COLORS=('#2676b8','#dd8b26','#34966f','#8758a5')


def main():
    summary=json.loads((ROOT/'summary.json').read_text())
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
        'axes.spines.right':False,'svg.fonttype':'none'})
    fig,axes=plt.subplots(1,3,figsize=(14,4.9),gridspec_kw={'width_ratios':[1.2,1.2,1.]})
    for i,(branch,name,color) in enumerate(zip(BRANCHES,LABELS,COLORS)):
        hs=[json.loads((ROOT/(branch+'_'+level)/'gpu-history.json').read_text()) for level in ('coarse','fine')]
        times=[];errors=[]
        for a,b in zip(*hs):
            if len(a['ids'])!=len(b['ids']):
                times.append(a['time']);errors.append(np.nan);continue
            x,y=np.asarray(a['spectrum']),np.asarray(b['spectrum'])
            times.append(a['time']);errors.append(100*np.linalg.norm(x-y)/(np.linalg.norm(y) or 1.))
        axes[0].plot(times,errors,color=color,label=name,lw=1.7)
        field=summary['time_rows'][i]['fields']
        axes[1].plot([f['time'] for f in field],[100*f['errors']['occupancy'] for f in field],
            color=color,marker=('o','s','^','D')[i],ms=4,lw=1.5)
        pairs=summary['space_rows'][2*i:2*i+2]
        axes[2].bar(np.arange(2)+(i-1.5)*.19,[100*r['field_max'] for r in pairs],width=.18,color=color)
    axes[0].set_title('Timestep: transport spectrum',loc='left',fontweight='bold')
    axes[1].set_title('Timestep: occupancy field',loc='left',fontweight='bold')
    axes[2].set_title('Spatial: maximum field error',loc='left',fontweight='bold')
    for ax in axes:
        ax.axhline(1,color='#444444',ls='--',lw=1)
        ax.set_ylabel('Relative discrepancy (%)')
        ax.set_ylim(bottom=0);ax.grid(axis='y',alpha=.16);ax.set_axisbelow(True)
    for ax in axes[:2]:
        ax.set_xlabel('Development time (dimensionless)');ax.set_xlim(0,18)
    axes[2].set_xticks([0,1],['72 → 88','88 → 104'])
    axes[2].set_xlabel('Grid refinement at dt = 0.001875')
    axes[0].legend(loc='upper left',fontsize=8.8,frameon=False)
    fig.suptitle('Repeated developmental qualification: backend agrees, refinement limits fail',
        x=.06,ha='left',fontsize=15,fontweight='bold')
    fig.text(.06,.035,'Dashed lines: preset 1% limits for these measurements. '
        'Same physical times; no temporal or spatial alignment. One synthetic seed, zero new independent histories.',fontsize=9,color='#555555')
    fig.subplots_adjust(top=.80,bottom=.23,left=.06,right=.98,wspace=.32)
    OUT.mkdir(parents=True,exist_ok=True)
    stem=OUT/'repeated-development-qualification-summary'
    fig.savefig(stem.with_suffix('.png'),dpi=180,facecolor='white')
    fig.savefig(stem.with_suffix('.svg'),facecolor='white')
    plt.close(fig)
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    sources=[ROOT/'summary.json',ROOT/'protocol.json',Path(__file__)]
    sources.extend(ROOT/(branch+'_'+level)/'gpu-history.json' for branch in BRANCHES for level in ('coarse','fine'))
    stem.with_suffix('.json').write_text(json.dumps({'source_sha256':{str(p.resolve()):sha(p) for p in sources},
        'artifact_sha256':{str(p.resolve()):sha(p) for p in (stem.with_suffix('.png'),stem.with_suffix('.svg'))},
        'scope':'Descriptive plot of unchanged production comparisons; independent review is separate.'},indent=2)+'\n')


if __name__=='__main__':main()
