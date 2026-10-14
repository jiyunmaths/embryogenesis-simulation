"""Plot independently reconstructed shared-context responses; no integration."""
from pathlib import Path
import csv
import hashlib
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

FOLDER = Path(__file__).resolve().parent

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    source=FOLDER/'independent_review.json'
    data=json.loads(source.read_text())
    assert data['passed'] and data['inputs_unchanged']
    rows=data['scientific_comparisons']
    pairs=list(dict.fromkeys((r['seed'],r['cell']) for r in rows))
    labels=[f'H{s}\ncell {c}' for s,c in pairs]
    colors={'S':'#326a8a','R':'#d99137'}
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(12.5,9.4),layout='constrained')
    x=np.arange(len(pairs));width=.35
    for ax,dose,title in zip(axes[0],('fractional','amount'),('A   Fractional depletion (10%)','B   Common removed amount')):
        selected=[next(r for r in rows if (r['seed'],r['cell'],r['dose'])==(*pair,dose)) for pair in pairs]
        for shift,P,name in ((-.5,'S','Shared sham physics'),(.5,'R','Shared redistributed physics')):
            ax.bar(x+shift*width,[r['chemical']['fine'][P]['rms'] for r in selected],width,color=colors[P],label=name)
        ax.axhline(.01,color='#777',ls='--',lw=1,label='Primary screen = 0.01')
        ax.set_yscale('log');ax.set_xticks(x,labels)
        ax.set_ylabel('Chemical-state response difference (RMS)')
        ax.set_title(title,loc='left',fontweight='bold')
        ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
        ax.legend(frameon=False,fontsize=8.5,loc='upper left')
        values=[r['chemical']['fine'][P]['rms'] for r in selected for P in ('S','R')]
        ax.set_ylim(min(values)*.55,max(.5,max(values)*5))
    ax=axes[1,0]
    for offset,dose in ((-.13,'fractional'),(.13,'amount')):
        selected=[next(r for r in rows if (r['seed'],r['cell'],r['dose'])==(*pair,dose)) for pair in pairs]
        for C,marker,color in (('S','o','#377857'),('R','s','#377857')):
            ax.scatter(x+offset,[r['physical']['fine'][C] for r in selected],marker=marker,s=37,
                       facecolors=color if dose=='fractional' else 'none',edgecolors=color,
                       label=f'Physical contrast, {C} chemistry, '+('fraction' if dose=='fractional' else 'amount'))
        ax.scatter(x+offset,[r['interaction']['fine'] for r in selected],marker='^',s=43,
                   facecolors='#925da1' if dose=='fractional' else 'none',edgecolors='#925da1',
                   label='Signed interaction, '+('fraction' if dose=='fractional' else 'amount'))
    ax.set_yscale('log');ax.set_xticks(x,labels)
    ax.set_title('C   Physical contrasts and signed interaction',loc='left',fontweight='bold')
    ax.set_ylabel('Response contrast RMS (descriptive)')
    ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
    ax.legend(frameon=False,fontsize=8,loc='upper left',ncol=2)
    ax.text(.98,.02,'S = sham; R = redistributed',transform=ax.transAxes,
            ha='right',va='bottom',fontsize=8)
    vals=[r['physical']['fine'][C] for r in rows for C in ('S','R')]+[r['interaction']['fine'] for r in rows]
    ax.set_ylim(min(v for v in vals if v>0)*.45,max(vals)*40)
    ax=axes[1,1]
    for shift,P,name in ((-.08,'S','Shared sham physics'),(.08,'R','Shared redistributed physics')):
        selected=[next(r for r in data['baseline_descriptions'] if (r['seed'],r['cell'],r['physical'])==(*pair,P)) for pair in pairs]
        ax.scatter(x+shift,[r['fine']['final_activator_ratio'] for r in selected],s=46,color=colors[P],label=name)
    ax.axhline(1.,color='#777',ls='--',lw=1)
    ax.set_yscale('log');ax.set_xticks(x,labels)
    ax.set_title('D   Unperturbed target states at t=570',loc='left',fontweight='bold')
    ax.set_ylabel('Final activator: redistributed / sham chemistry')
    ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
    ax.legend(frameon=False,fontsize=9,loc='upper left')
    vals=[r['fine']['final_activator_ratio'] for r in data['baseline_descriptions']]
    ax.set_ylim(min(vals)*.65,max(vals)*2)
    fig.suptitle('Chemical-state responses in reciprocal shared physical contexts\n'
                 'Actual t=510 concentrations; fine dt=0.001875; six recipients nested in three histories',
                 fontsize=14,fontweight='bold')
    for suffix in ('png','svg'):
        fig.savefig(FOLDER/f'shared-context-response.{suffix}',dpi=180,
                    metadata={'Description':'From independent saved-array reconstruction. Chemistry includes concentration levels and distribution; reassignment preserves concentrations, not amounts. No statistical replication by paths.'})
    plt.close(fig)
    exported=[]
    for row in rows:
        exported.append(dict(history=row['seed'],cell=row['cell'],dose=row['dose'],
            chemical_on_sham_physics=row['chemical']['fine']['S']['rms'],
            chemical_on_redistributed_physics=row['chemical']['fine']['R']['rms'],
            chemical_on_sham_above_screen=row['chemical']['fine']['S']['detected'],
            chemical_on_redistributed_above_screen=row['chemical']['fine']['R']['detected'],
            physical_with_sham_chemistry=row['physical']['fine']['S'],
            physical_with_redistributed_chemistry=row['physical']['fine']['R'],
            signed_interaction_rms=row['interaction']['fine']))
    with (FOLDER/'recipient-summary.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(exported[0]));writer.writeheader();writer.writerows(exported)
    report=dict(input_sha256=sha(source),script_sha256=sha(Path(__file__)),
        outputs_sha256={name:sha(FOLDER/name) for name in ('shared-context-response.png','shared-context-response.svg','recipient-summary.csv')},
        matplotlib=matplotlib.__version__,scope='Fine full60 own-control-normalized signed two-species response contrasts, plus final unperturbed target activator ratios. Chemical screens shown only in primary panels; physical/interaction metrics descriptive. Three existing histories, no error bars or new replication.')
    (FOLDER/'figure-provenance.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
