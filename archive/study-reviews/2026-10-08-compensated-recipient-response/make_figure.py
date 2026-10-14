"""Plot reviewed delayed responses and the resolved roundoff control."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path('output/studies/compensated-recipient-amount-delayed-response')
ARCH=Path(__file__).parent

def main():
    r=json.loads((ARCH/'independent_review.json').read_text());assert r['passed']
    immediate=json.loads(Path('docs/phase_carry_recipient_amount_response_assessment.json').read_text())
    fig,axes=plt.subplots(1,2,figsize=(11,5.2),layout='constrained')
    ax=axes[0];cells=[(7,22),(7,23),(8,19),(8,20),(9,29),(9,30)]
    for offset,rows,label,color in [(-.12,immediate['scientific_comparisons'],'Immediate reset response','#d97706'),(.12,r['scientific_comparisons'],'After a 60-unit wait','#0369a1')]:
        centers=[];lower=[];upper=[]
        for seed,cell in cells:
            values=[x['full_effects']['fine'] for x in rows if (x['seed'],x['cell'])==(seed,cell) and x.get('dose','fractional')=='fractional']
            lo,hi=min(values),max(values);center=(lo+hi)/2
            centers.append(center);lower.append(center-lo);upper.append(hi-center)
        ax.errorbar(np.arange(6)+offset,centers,yerr=[lower,upper],fmt='o',capsize=4,color=color,label=label)
    ax.axhline(.01,color='#475569',ls='--',lw=1,label='Prespecified response screen')
    ax.set(yscale='log',xticks=np.arange(6),xticklabels=[f'H{h}\nCell {c}' for h,c in cells],ylabel='Normalized response contrast, E',title='Responses before and after waiting')
    ax.legend(fontsize=8,loc='upper center',bbox_to_anchor=(.5,-.16));ax.grid(axis='y',alpha=.2)
    key=r['diagnosed_control']['new_order']['trial'];ax=axes[1]
    for root,label,color,ls in [(Path('output/studies/phase-carry-recipient-amount-delayed-response'),'Original fine calculation','#94a3b8','--'),(ROOT,'Compensated fine calculation','#0369a1','-')]:
        with np.load(root/'production/fine/paths.npz') as z:
            i=list(z['trial_keys']).index(key);u=z['trajectory'][:,i];times=z['times'];drift=np.max(abs(np.log(u/u[0])),axis=(1,2))
        ax.plot(times,drift/1e-13,label=label,color=color,ls=ls)
    with np.load(ROOT/'references'/key/'paths.npz') as z:
        u=z['dop853'];ax.plot(z['times'],np.max(abs(np.log(u/u[0])),axis=(1,2))/1e-13,label='Independent DOP853 reference',color='#7c3aed',ls=':',lw=2)
    ax.set(xlabel='Time after delayed challenge',ylabel='Maximum log relaxation / 1e-13',title='Tiny control relaxation is resolved')
    ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.suptitle('Initial-amount control: compensated delayed repeat',fontsize=14)
    fig.savefig(ARCH/'compensated-delayed-response.png',dpi=170)
    fig.savefig(ARCH/'compensated-delayed-response.svg')
    above=[str(c) for h,c in cells if all(x['above_screen'] for x in r['scientific_comparisons'] if (x['seed'],x['cell'])==(h,c))]
    caption=('Left: reviewed fractional-dose response contrasts. Ranges span both original matching references and both shared physical contexts; numerical repetitions are not independent histories. Immediate and delayed assays use their own continuing unpulsed controls. Cells '+('/'.join(above))+' remain above the strict E>0.01 screen after waiting; below-screen responses do not establish equivalence. Right: the previously stuck near-equilibrium control now follows independent-reference relaxation. Its convergence order remains unresolved at the unchanged reference-noise floor. Three existing mature histories, no new developmental histories. Original failed study retained.')
    (ARCH/'figure-caption.md').write_text(caption+'\n')

if __name__=='__main__':main()
