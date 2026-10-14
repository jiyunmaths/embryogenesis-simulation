"""Report-only diagnostics from independently reviewed saved data."""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path('output/studies/development-time-diagnosis')
HERE = Path(__file__).parent
review = json.loads((HERE/'independent_review.json').read_text())
names = [p['case'] for p in review['paths']]
traces = {}; lineages = {}; inputs = {str(HERE/'independent_review.json'):
    hashlib.sha256((HERE/'independent_review.json').read_bytes()).hexdigest()}
for name in names:
    trace_path = ROOT/name/'network.jsonl'; lineage_path = ROOT/name/'gpu-lineage.json'
    for path in (trace_path, lineage_path):
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        assert review['consumed_sha256'][str(path.resolve())] == digest
        inputs[str(path)] = digest
    traces[name] = [json.loads(l) for l in trace_path.open()]
    lineages[name] = json.loads(lineage_path.read_text())
rows = []
for name in names:
    first = lineages[name][0]
    trace = next(r for r in traces[name] if abs(r['time']-3.) < 1e-9)
    active = next(a for a in trace['active_divisions'] if a['parent'] == 0)
    rows.append(dict(case=name, first_start=first['division_start'], first_abscission=first['division'],
        first_duration=first['division']-first['division_start'], first_axis=first['division_axis'],
        time3_cells=len(trace['ids']), time3_edges=trace['edges'], time3_neck=active['neck'],
        time3_overlap=active['overlap'], time3_progress=active['progress']))
assert all(r['time3_cells'] == 1 and r['time3_edges'] == 0 for r in rows)
parent_rows = []
for i,(a,b) in enumerate(zip(names,names[1:])):
    for x,y in zip(lineages[a],lineages[b]):
        if 'division_start' not in x: continue
        assert x['id'] == y['id'] and x['parent'] == y['parent']
        parent_rows.append(dict(pair=[a,b],parent=x['id'],
            start_difference=x['division_start']-y['division_start'],
            abscission_difference=x['division']-y['division'],
            furrow_duration_difference=(x['division']-x['division_start'])-(y['division']-y['division_start'])))
data = dict(first_furrow=rows, parent_clocks=parent_rows,
    adjacent_occupancy=[r['field_max'] for r in review['pairs']],
    adjacent_spectrum=[r['errors']['relative_spectrum_max'] for r in review['pairs']],
    dense_spectrum=[r['dense_spectrum_peak'] for r in review['pairs']],
    inputs_sha256=inputs,
    interpretation='Errors exist before first abscission and first cell-cell graph; cutoff changes cannot be the sole origin. Pressure correction/preprojection errors halve with dt, but no isolated projection intervention was made. Error decay is descriptive, not a proved order or continuum limit.')
(HERE/'cleavage_diagnosis.json').write_text(json.dumps(data,indent=2)+'\n')
fig,axs = plt.subplots(1,3,figsize=(14.4,4.5))
x = np.arange(3); labels = ['.00375\n→ .001875','.001875\n→ .0009375','.0009375\n→ .00046875']
axs[0].bar(x-.18,100*np.array(data['adjacent_occupancy']),.36,label='Occupancy field')
axs[0].bar(x+.18,100*np.array(data['adjacent_spectrum']),.36,label='Transport spectrum')
axs[0].axhline(1,color='#bc352a',ls='--',label='Unchanged 1% limits')
axs[0].set_xticks(x,labels);axs[0].set_ylabel('Maximum discrepancy (%)')
axs[0].set_title('Adjacent timestep screen');axs[0].legend(fontsize=8)
a,b = names[-2:]
regular = [{round(r['time'],9):r for r in traces[k]
    if abs(r['time']/.015-round(r['time']/.015))<1e-7} for k in (a,b)]
t = sorted(regular[0]); er = []; missing = []
for time in t:
    x,y = regular[0][time],regular[1][time]
    if len(x['ids']) != len(y['ids']): er.append(np.nan);missing.append(time);continue
    sx,sy = np.array(x['spectrum']),np.array(y['spectrum'])
    er.append(100*np.linalg.norm(sx-sy)/(np.linalg.norm(sy) or 1.))
axs[1].plot(t,er,color='#2c7589',lw=1.4)
axs[1].axhline(1,color='#bc352a',ls='--')
axs[1].scatter(missing,np.zeros(len(missing)),marker='x',s=24,color='#555',label='Unequal cell counts')
axs[1].set_ylim(-.04,1.08);axs[1].set_xlim(0,18)
axs[1].set_xlabel('Model time');axs[1].set_ylabel('Relative spectrum discrepancy (%)')
axs[1].set_title('Finest pair: dense diagnostic');axs[1].legend(fontsize=8)
for name,label in zip(names,['.00375','.001875','.0009375','.00046875']):
    ts=[];ns=[]
    for row in traces[name]:
        for active in row['active_divisions']:
            if active['parent'] == 0: ts.append(row['time']);ns.append(active['neck'])
    first = lineages[name][0]
    ts.append(first['division']); ns.append(first['neck_at_abscission'])
    axs[2].plot(ts,ns,label='dt='+label,lw=1.3)
axs[2].axhline(.15,color='#bc352a',ls='--',label='Neck limit')
axs[2].set_xlim(3.10,3.21);axs[2].set_ylim(.08,1.05)
axs[2].set_xlabel('Model time');axs[2].set_ylabel('First-mother neck diagnostic')
axs[2].set_title('First furrow, before cell-cell graph');axs[2].legend(fontsize=8)
fig.suptitle('Feedback off • same synthetic seed 7 • grid 72 • fifteen natural divisions',fontsize=12)
fig.tight_layout()
target = Path('docs/figures/development-time-diagnosis')
fig.savefig(target.with_suffix('.png'),dpi=180)
fig.savefig(target.with_suffix('.svg'))
plt.close(fig)
print(json.dumps({'first_furrow':rows,'parent_clock_max':[r['parent_abscission_time_max'] for r in review['pairs']],
    'angles_max':[r['axis_angle_degrees_max'] for r in review['pairs']], 'figure':str(target.with_suffix('.png'))},indent=2))
