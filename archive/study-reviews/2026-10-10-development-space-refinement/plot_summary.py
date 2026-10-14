"""Standalone scientific figure from independently reviewed resolution data."""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path('output/studies/development-space-refinement')
HERE = Path(__file__).parent
source = HERE/'independent_review.json'
r = json.loads(source.read_text()); inputs = {str(source):hashlib.sha256(source.read_bytes()).hexdigest()}
spaces = r['space_rows']; temporal = r['time_row']
regular = []
for name in temporal['cases']:
    path = ROOT/name/'network.jsonl'; digest = hashlib.sha256(path.read_bytes()).hexdigest()
    assert r['consumed_sha256'][str(path.resolve())] == digest
    inputs[str(path)] = digest
    records = [json.loads(l) for l in path.open()]
    regular.append({round(x['time'],9):x for x in records if abs(x['time']/.015-round(x['time']/.015))<1e-7})
times = sorted(regular[0]); errors = []
for time in times:
    a,b = regular[0][time],regular[1][time]
    if len(a['ids']) != len(b['ids']): errors.append(None); continue
    x,y = np.array(a['spectrum']),np.array(b['spectrum'])
    errors.append(float(np.linalg.norm(x-y)/(np.linalg.norm(y) or 1.)))
data = dict(spatial_cases=[x['cases'] for x in spaces],
    spatial_field_max=[x['field_max'] for x in spaces],
    spatial_spectrum_max=[x['errors']['relative_spectrum_max'] for x in spaces],
    spatial_dense_spectrum_max=[x['dense']['peak']['relative_spectrum'] for x in spaces],
    temporal_field_checkpoints=temporal['field_checks'],dense_times=times,dense_temporal_spectrum=errors,
    dense_unmatched_times=temporal['dense']['unmatched_times'],inputs_sha256=inputs,
    scope='One synthetic seed7, feedback off, grid120 time check and adjacent spatial screen. No fresh developmental histories.')
target = Path('docs/figures/development-space-refinement-summary')
target.with_suffix('.json').write_text(json.dumps(data,indent=2)+'\n')
fig,axs = plt.subplots(2,2,figsize=(12.4,7.6))
x = np.arange(3); labels = ['72 → 88','88 → 104','104 → 120']
values = 100*np.array(data['spatial_field_max'])
axs[0,0].bar(x,values,color=['#657e9a','#657e9a','#287b72'],width=.6)
axs[0,0].axhline(1,color='#b73e32',ls='--',label='1% field limit')
for a,v in zip(x,values): axs[0,0].text(a,v+.07,f'{v:.3f}%',ha='center',fontsize=10)
axs[0,0].set_xticks(x,labels);axs[0,0].set_ylim(0,3.75)
axs[0,0].set_title('Spatial field discrepancies decrease');axs[0,0].set_ylabel('Maximum relative L2 discrepancy (%)')
axs[0,0].legend(fontsize=9)
values = 100*np.array(data['spatial_spectrum_max'])
axs[0,1].bar(x,values,color=['#657e9a','#657e9a','#287b72'],width=.6,label='Original observations')
axs[0,1].plot(x,100*np.array(data['spatial_dense_spectrum_max']),'D',color='#222',ms=5,label='Dense diagnostic maxima')
axs[0,1].axhline(5,color='#b73e32',ls='--',label='5% spatial spectrum limit')
for a,v in zip(x,values): axs[0,1].text(a,v+.27,f'{v:.3f}%',ha='center',fontsize=10)
axs[0,1].set_xticks(x,labels);axs[0,1].set_ylim(0,13.8)
axs[0,1].set_title('Finest spatial spectrum pair passes');axs[0,1].set_ylabel('Relative spectrum discrepancy (%)')
axs[0,1].legend(fontsize=8,loc='upper right')
cp = temporal['field_checks']
axs[1,0].plot([f['time'] for f in cp],[100*f['errors']['occupancy'] for f in cp],'o-',color='#287b72',label='Occupancy')
axs[1,0].axhline(1,color='#b73e32',ls='--',label='1% temporal field limit')
axs[1,0].set_ylim(0,1.08);axs[1,0].set_xlim(0,18);axs[1,0].set_xticks(range(0,19,3))
axs[1,0].set_title('Grid 120: native field checkpoints');axs[1,0].set_xlabel('Model time')
axs[1,0].set_ylabel('Relative L2 discrepancy (%)');axs[1,0].legend(fontsize=9)
values = [np.nan if e is None else 100*e for e in errors]
axs[1,1].plot(times,values,color='#287b72',lw=1.3)
axs[1,1].axhline(1,color='#b73e32',ls='--',label='1% temporal spectrum limit')
missing = temporal['dense']['unmatched_times']
axs[1,1].scatter(missing,np.zeros(len(missing)),marker='x',color='#555',s=28,label='Unequal cell counts')
axs[1,1].set_ylim(-.04,1.08);axs[1,1].set_xlim(0,18)
axs[1,1].set_title('Grid 120: dense temporal spectra');axs[1,1].set_xlabel('Model time')
axs[1,1].set_ylabel('Relative spectrum discrepancy (%)');axs[1,1].legend(fontsize=8)
fig.suptitle('Feedback-off repeated development: spatial and temporal screens\nSame synthetic seed 7 • t=0–18 • four new paths and one reused path',fontsize=13)
fig.tight_layout()
fig.savefig(target.with_suffix('.png'),dpi=170)
fig.savefig(target.with_suffix('.svg'))
plt.close(fig)
print(target.with_suffix('.png'))
