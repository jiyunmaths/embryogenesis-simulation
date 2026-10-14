"""Show all accepted amount-matching comparisons without new integration."""
import csv,hashlib,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
ROOT=Path(__file__).resolve().parent
sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()
new_files=[ROOT.parent/'2026-10-08-recipient-amount-response/independent_review.json',ROOT/'independent_review.json']
old_files=[Path('docs/phase_carry_recipient_state_response_assessment.json'),Path('docs/phase_carry_recipient_delayed_response_assessment.json')]
reports=[json.loads(f.read_text()) for f in new_files];old=[json.loads(f.read_text()) for f in old_files]
assert reports[0]['passed'] and all(r['inputs_unchanged'] for r in reports) and all(r['full_study_accepted'] for r in old)
assert len([x for x in reports[1]['orders'] if not x['passed']])==1
assert all(x['passed'] for k in ('production_accuracy','raw_timestep_comparisons','response_comparisons','duplicate_continuity') for x in reports[1][k])
rows=[]
for phase,new,prior in zip(('immediate','delayed'),reports,old):
 for r in new['scientific_comparisons']:
  dose=r.get('dose','fractional');q=next(q for q in prior['scientific_comparisons'] if all(q[k]==r[k] for k in ('seed','cell','matching_reference','physics')) and q.get('dose','fractional')==dose)
  rows.append(dict(phase=phase,history=r['seed'],cell=r['cell'],matching_reference=r['matching_reference'],physics=r['physics'],dose=dose,unmatched_rms=q['full_effects']['fine'],amount_matched_rms=r['full_effects']['fine'],unmatched_above=q['above_screen'],matched_above=r['above_screen']))
cells=[(7,22),(7,23),(8,19),(8,20),(9,29),(9,30)];combos=[('S','S'),('S','R'),('R','S'),('R','R')]
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
fig,axs=plt.subplots(2,3,figsize=(14.5,9.2));fig.subplots_adjust(left=.08,right=.91,top=.85,bottom=.19,hspace=.49,wspace=.34)
fig.suptitle('Does chemical allocation matter at a common initial amount?',fontsize=17,y=.97)
fig.text(.5,.93,'Local chemistry and both species totals matched once at t=510. Amounts then evolve freely on shared frozen physics.',ha='center',fontsize=10)
columns=[('immediate','fractional','Immediate: 10% depletion'),('delayed','fractional','After 60 units: 10% depletion'),('delayed','amount','After 60 units: common amount')]
for j,(phase,dose,title) in enumerate(columns):
 for i,(value,label) in enumerate([('unmatched_rms','Original unmatched amounts'),('amount_matched_rms','Initially matched amounts')]):
  mat=np.array([[next(r[value] for r in rows if (r['phase'],r['dose'],r['history'],r['cell'],r['matching_reference'],r['physics'])==(phase,dose,s,c,L,P)) for L,P in combos] for s,c in cells])
  ax=axs[i,j];im=ax.imshow(mat,norm=LogNorm(vmin=1e-5,vmax=.2),cmap='viridis',aspect='auto')
  ax.set_xticks(range(4),[f'Match {L}\nPhysics {P}' for L,P in combos],fontsize=9);ax.set_yticks(range(6),[f'H{s}, cell {c}' for s,c in cells],fontsize=9)
  ax.set_title(title if i==0 else 'Same challenge, initial amounts matched',loc='left',fontsize=10.5)
  if j==0:ax.set_ylabel(label)
  for y in range(6):
   for x in range(4):
    v=mat[y,x];ax.text(x,y,f'{v:.4f}'+(' *' if v>.01 else ''),ha='center',va='center',fontsize=8.7,color='black' if v>.03 else 'white',fontweight='bold' if v>.01 else 'normal')
  ax.spines[['bottom','left']].set_visible(False)
cax=fig.add_axes([.94,.34,.015,.36]);cb=fig.colorbar(im,cax=cax);cb.set_label('Two-species normalized response RMS',fontsize=10)
fig.text(.08,.115,'* Above the prespecified screen E>0.01. All matching references, physical contexts and doses shown.',fontsize=10)
fig.text(.08,.085,'Three existing histories; six recipients are nested measurements. All scientific decisions must agree across five solvers/timesteps.',fontsize=10)
fig.text(.08,.055,'Delayed challenges use actual unpulsed endpoints without rematching. This isolates an initial budget intervention, not continuously fixed dosage.',fontsize=10)
fig.text(.08,.025,'Delayed phase remains provisional: one near-stationary control fails its order gate at ~1.6e-13. No autonomous identity claim.',fontsize=10)
for ext in ('png','svg'):fig.savefig(ROOT/f'initial-amount-control.{ext}',dpi=180)
plt.close(fig)
file=ROOT/'comparison-summary.csv'
with file.open('w',newline='') as stream:
 writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
# Show endogenous amount changes directly rather than implying continued budget control.
fig,axs=plt.subplots(1,3,figsize=(14,4.6),sharey=True);fig.subplots_adjust(left=.08,right=.97,top=.78,bottom=.22,wspace=.2)
fig.suptitle('Initially equal amounts can diverge through reactions',fontsize=16,y=.95)
for ax,seed in zip(axs,(7,8,9)):
 for cell,color in zip([c for s,c in cells if s==seed],('#276da4','#db8629')):
  r=next(r for r in reports[0]['amount_evolution'] if (r['seed'],r['cell'],r['matching_reference'],r['physics'])==(seed,cell,'S','S'))
  ratios=np.array(r['global_species_amount_ratios']);t=np.array(r['times'])
  for species,style in ((0,'-'),(1,'--')):ax.plot(t,np.log(ratios[:,species]),color=color,ls=style,label=f'Cell {cell}: '+('activator' if species==0 else 'inhibitor'))
 ax.axhline(0,color='#777',lw=.8);ax.set_title(f'History {seed}',loc='left');ax.set_xlabel('Elapsed time after initial matching');ax.legend(frameon=False,fontsize=8);ax.grid(alpha=.13)
axs[0].set_ylabel('Log whole-network amount ratio\n(redistributed / sham)')
fig.text(.08,.07,'Named examples: sham matching reference and physics. All 24 paired amount trajectories are saved and independently reconstructed.',fontsize=10)
for ext in ('png','svg'):fig.savefig(ROOT/f'amount-evolution.{ext}',dpi=180)
plt.close(fig)
files=[Path(__file__),*new_files,*old_files,file,*[ROOT/(stem+'.'+ext) for stem in ('initial-amount-control','amount-evolution') for ext in ('png','svg')]]
(ROOT/'figure-provenance.json').write_text(json.dumps(dict(passed=True,no_new_integrations=True,no_new_acceptance_threshold=True,comparisons=len(rows),delayed_full_study_accepted=False,scope='All72 scientific comparisons plotted; named initial-amount curves are descriptive. Delayed phase retains one failed order gate.',files_sha256={str(f):sha(f) for f in files}),indent=2)+'\n')
print('All72 response comparisons plotted; named amount examples use the24 archived paired curves.')
