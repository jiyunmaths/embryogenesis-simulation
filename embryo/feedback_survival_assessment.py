"""Immutable matched-time assessment of an ongoing feedback-switch study."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from .resolution import write_json


def assess(source,output=None,until=None):
    source=Path(source);protocol=json.loads((source/'protocol.json').read_text())
    histories={arm:json.loads((source/arm/'history.json').read_text()) for arm in ['switch_on','keep_off']}
    dt=protocol['dt']
    by_step={arm:{round(row['metrics']['time']/dt):row for row in history} for arm,history in histories.items()}
    common=set.intersection(*(set(rows) for rows in by_step.values()))
    if until is not None:common={step for step in common if step*dt<=until+1e-9}
    if not common:raise ValueError('No common observation time')
    end=max(common);time=end*dt
    output=Path(output) if output is not None else source/f'interim-{time:g}'
    output.mkdir(parents=True,exist_ok=False)
    prefix={arm:[row for row in history if round(row['metrics']['time']/dt)<=end] for arm,history in histories.items()}
    rows={}
    for arm,history in prefix.items():
        write_json(output/(arm+'-history.json'),history)
        last=by_step[arm][end]
        stable=[row['metrics']['time'] for row in history if row['spectrum']['maximum_spatial_growth']<0]
        rows[arm]=dict(time=last['metrics']['time'],log_activator_sd=last['attribute_std'][0],
            initial_state_correlation=last['similarity']['initial_log_activator_correlation'],
            chemical_distance_from_initial=last['similarity']['initial_chemical_log_rms'],
            homogeneous_growth_rate=last['spectrum']['maximum_spatial_growth'],
            first_sampled_negative_growth=min(stable) if stable else None,
            aggregate_axis_ratio=last['metrics']['axis_ratio'],
            max_sampled_cell_volume_error=max(row['metrics']['max_cell_volume_error'] for row in history),
            max_sampled_boundary=max(row['metrics']['boundary_occupancy'] for row in history))
    late_start=protocol['start']+protocol['duration']-protocol['late_duration']
    finished=time>=protocol['start']+protocol['duration']-1e-9
    summary=dict(matched_time=time,final_time=protocol['start']+protocol['duration'],
        late_window_start=late_start,late_window_reached=time>=late_start,full_horizon_reached=finished,
        arms=rows,protocol_sha256=hashlib.sha256((source/'protocol.json').read_bytes()).hexdigest(),
        scope='Matched history prefixes; sampled diagnostics only. Interim state does not satisfy the final survival assessment.',
        interpretation='An established pattern can remain heterogeneous after the instantaneous homogeneous-state growth rate becomes negative; long-horizon persistence is still pending.' if not finished else 'Full-horizon prefixes available; consult the production comparison and every-step audits.')
    write_json(output/'assessment.json',summary)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(12,3.7))
    for arm,history in prefix.items():
        t=[r['metrics']['time'] for r in history]
        axes[0].plot(t,[r['attribute_std'][0] for r in history],'-o',label=arm,ms=3)
        axes[1].plot(t,[r['spectrum']['maximum_spatial_growth'] for r in history],'-o',label=arm,ms=3)
        axes[2].plot(t,[r['similarity']['initial_log_activator_correlation'] for r in history],'-o',label=arm,ms=3)
    axes[1].axhline(0,color='grey',ls='--')
    correlations=[r['similarity']['initial_log_activator_correlation'] for h in prefix.values() for r in h]
    correlations=[x for x in correlations if x is not None]
    if correlations:
        lower=min(correlations);axes[2].set_ylim(max(-1.01,lower-max(1e-5,.1*(1-lower))),1.000005)
    for ax,label in zip(axes,['SD of log activator','Homogeneous-state growth rate','Correlation with original cell states']):
        ax.set(xlabel='Model time',ylabel=label);ax.legend(fontsize=8)
    fig.suptitle(f'Interim matched comparison through t={time:g}; final assessment t={late_start:g}–{summary["final_time"]:g}')
    fig.tight_layout();fig.savefig(output/'assessment.png',dpi=160);plt.close(fig)
    lines=[f'# Interim feedback-switch assessment through t={time:g}','','**Final survival outcome remains pending.**' if not finished else 'Full-horizon data available.',
        '',f'The declared late window is t={late_start:g}–{summary["final_time"]:g}.','',
        '| Arm | Log-activator SD | Initial-state correlation | Homogeneous growth rate |',
        '|---|---:|---:|---:|']
    lines += [f'| {arm} | {row["log_activator_sd"]:.6f} | {row["initial_state_correlation"]:.8f} | {row["homogeneous_growth_rate"]:+.6f} |' for arm,row in rows.items()]
    lines += ['',summary['interpretation'],'','![Matched interim trajectories](assessment.png)']
    (output/'ASSESSMENT.md').write_text('\n'.join(lines)+'\n')
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,default=Path('outputs/feedback-survival'))
    p.add_argument('--output',type=Path);p.add_argument('--until',type=float)
    a=p.parse_args();print(json.dumps(assess(a.source,a.output,a.until),indent=2))
