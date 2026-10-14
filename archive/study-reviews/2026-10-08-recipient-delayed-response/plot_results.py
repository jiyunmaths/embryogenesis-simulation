"""Plot accepted immediate and delayed measures without further integration."""
import csv
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

ROOT=Path(__file__).resolve().parent


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


report=ROOT/'independent_review.json';waves_path=ROOT/'waveforms.npz'
r=json.loads(report.read_text());assert r['passed'] and r['inputs_unchanged']
assert sha(waves_path)==r['waveform_sha256']
w=np.load(waves_path)
cells=[(7,22),(7,23),(8,19),(8,20),(9,29),(9,30)]
combinations=[('S','S'),('S','R'),('R','S'),('R','R')]
effect={(x['seed'],x['cell'],x['matching_reference'],x['physics'],x['dose']):x for x in r['scientific_comparisons']}
immediate={(x['seed'],x['cell'],x['matching_reference'],'S' if x['physics']=='sham_physics' else 'R'):x for x in r['immediate_vs_delayed']}
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
fig,axs=plt.subplots(2,2,figsize=(12,9.5))
fig.subplots_adjust(left=.08,right=.90,top=.87,bottom=.15,hspace=.52,wspace=.36)
fig.suptitle('Cell responses after a fixed wait following initial state matching',fontsize=16,y=.965)
fig.text(.5,.92,'Immediate pulse: t=510. Delayed pulse: t=570, after 60 units of free chemical evolution; no second reset.',ha='center',fontsize=10)
matrices=[np.array([[immediate[s,c,L,P]['immediate_rms'] for L,P in combinations] for s,c in cells]),
          np.array([[effect[s,c,L,P,'fractional']['full_effects']['fine'] for L,P in combinations] for s,c in cells]),
          np.array([[effect[s,c,L,P,'amount']['full_effects']['fine'] for L,P in combinations] for s,c in cells])]
titles=['A  Immediate: fractional depletion','B  Delayed: fractional depletion','C  Delayed: equal removed amount']
for ax,values,title in zip((axs[0,0],axs[0,1],axs[1,0]),matrices,titles):
    im=ax.imshow(values,cmap='viridis',norm=LogNorm(vmin=.0001,vmax=.15),aspect='auto')
    ax.set_xticks(range(4),[f'Match {L}\nPhysics {P}' for L,P in combinations],fontsize=9)
    ax.set_yticks(range(6),[f'H{s}, cell {c}' for s,c in cells])
    ax.set_title(title,loc='left',fontsize=11)
    for i in range(6):
        for j in range(4):
            value=values[i,j]
            ax.text(j,i,f'{value:.4f}'+(' *' if value>.01 else ''),ha='center',va='center',
                    color='black' if value>.025 else 'white',fontsize=9,fontweight='bold' if value>.01 else 'normal')
    ax.spines[['bottom','left']].set_visible(False)
cax=fig.add_axes([.92,.57,.018,.28])
colorbar=fig.colorbar(im,cax=cax)
colorbar.set_label('Normalized two-species response RMS',fontsize=9)

ax=axs[1,1];times=w['times']-w['times'][0]
parent=ROOT.parent/'2026-10-08-recipient-state-response/waveforms.npz'
old=np.load(parent)
for color,(seed,cell) in zip(['#2367a1','#df852e'],[(8,20),(9,30)]):
    key=f'{seed}_{cell}_S_fractional_chemical_sham_physics'
    ax.plot(times,old[key][:,0],color=color,ls='--',label=f'Cell {cell}, immediate')
    ax.plot(times,w[key][:,0],color=color,ls='-',label=f'Cell {cell}, delayed')
ax.axhline(0,color='#777777',lw=.7);ax.set_xlim(0,12)
ax.set_xlabel('Elapsed time after each pulse')
ax.set_ylabel('Signed activator response difference\n(redistributed − sham; normalized)')
ax.set_title('D  Response examples: sham match and physics',loc='left',fontsize=10.5)
ax.legend(fontsize=8,frameon=False);ax.grid(alpha=.12)
fig.text(.08,.077,'* Above the prespecified primary screen (RMS > 0.01). S/R: original sham/redistributed matching reference or shared physics.',fontsize=9)
fig.text(.08,.052,'Three existing mature histories; targets are nested within histories. All backgrounds pass the predefined settling screens.',fontsize=9)
fig.text(.08,.027,'Target concentrations and surrounding species totals can differ at the delayed pulse. Frozen mechanics; no new developmental histories.',fontsize=9)
for suffix in ('png','svg'):fig.savefig(ROOT/f'recipient-delayed-response.{suffix}',dpi=180)
plt.close(fig)
csv_path=ROOT/'recipient-summary.csv'
fields=['history','cell','original_matching_reference','shared_physics','dose','delayed_rms','above_primary_screen','immediate_fractional_rms','fine_target_recovery_S','fine_target_recovery_R','fine_network_recovery_S','fine_network_recovery_R']
with csv_path.open('w',newline='') as stream:
    writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader()
    for x in r['scientific_comparisons']:
        s,c,L,P=x['seed'],x['cell'],x['matching_reference'],x['physics'];m=x['fine_pulse_metrics']
        writer.writerow(dict(history=s,cell=c,original_matching_reference=L,shared_physics=P,dose=x['dose'],
            delayed_rms=x['full_effects']['fine'],above_primary_screen=x['above_screen'],
            immediate_fractional_rms=immediate[s,c,L,P]['immediate_rms'] if x['dose']=='fractional' else '',
            fine_target_recovery_S=m['S']['target_recovery_time'],fine_target_recovery_R=m['R']['target_recovery_time'],
            fine_network_recovery_S=m['S']['network_recovery_time'],fine_network_recovery_R=m['R']['network_recovery_time']))
files=[Path(__file__),report,waves_path,parent,csv_path,ROOT/'recipient-delayed-response.png',ROOT/'recipient-delayed-response.svg']
(ROOT/'figure-provenance.json').write_text(json.dumps(dict(passed=True,no_new_integrations=True,no_new_acceptance_threshold=True,
    scope='All primary reference/context/dose effects shown; panel D uses two named examples, not additional scientific tests.',
    files_sha256={str(p):sha(p) for p in files}),indent=2)+'\n')
print('Saved PNG/SVG, 48-row comparison CSV and figure provenance.')
