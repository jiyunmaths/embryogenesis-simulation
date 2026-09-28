"""Plot finite reservoir formation, release, and reservoir adaptation."""
import argparse
import json
from pathlib import Path
import numpy as np


def summarize(output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(13,4))
    for scale,color in zip([.25,1,4],['tab:blue','tab:orange','tab:green']):
        data=np.load(output/f'formation_s{scale:g}_R1_seed0.npz')
        t=data['times'];path=data['trajectory']
        axes[0].plot(t,np.std(np.log(path[:,0,:-1]),axis=1),label=f's = {scale:g}',color=color)
        filename=f'release_unit_{scale:g}_R1.npz'
        if scale==4:filename='verified_'+filename
        data=np.load(output/filename);t=data['times'];path=data['trajectory']
        axes[1].plot(t,np.std(np.log(path[:,0,:-1]),axis=1),color=color,label=f's = {scale:g}')
        axes[2].plot(t,path[:,0,-1],color=color,label=f'A, s = {scale:g}')
        axes[2].plot(t,path[:,1,-1],color=color,ls='--',label=f'B, s = {scale:g}')
    axes[0].set(title='Formation from small perturbations',xlabel='Time',ylabel='SD of log activator',xlim=(0,160))
    axes[1].set(title='Release of established differences',xlabel='Time',ylabel='SD of log activator',xlim=(0,1200))
    axes[2].set(title='Reservoir responds to cells',xlabel='Time',ylabel='Reservoir concentrations',xlim=(0,1200))
    for ax in axes[1:]:
        ax.set_xscale("symlog",linthresh=10)
        ax.set_xlabel("Time (logarithmic beyond 10)")
    for ax in axes:ax.legend(fontsize=8)
    fig.suptitle('Reservoir volume equals total cell volume; strong release includes small perturbation')
    fig.tight_layout();fig.savefig(output/'finite-reservoir.png',dpi=160);plt.close(fig)
    rows=json.loads((output/'results.json').read_text())
    check=json.loads((output/'verification.json').read_text())
    summary=dict(baseline_numerical_failures=[r['name'] for r in rows if not r['numerical_pass']],
        baseline_nonstationary=[r['name'] for r in rows if not r['stationary']],
        independent_failures=[r['name'] for r in check if not r['numerical_pass']],
        verified_stable_persistent=sum(r['numerical_pass'] and r['persistent'] and r['stable_stationary'] for r in check))
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=Path('outputs/attribute-finite-reservoir'))
    summarize(p.parse_args().output)
