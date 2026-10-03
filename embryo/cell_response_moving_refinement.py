"""Same-state timestep halving for the moving single-cell response pilot."""
import argparse
import copy
import json
from pathlib import Path
import shutil
import numpy as np
from . import cell_response_moving as moving
from .feedback_long import digest
from .feedback_survival_validation import retime
from .resolution import write_json

CRITERIA = dict(max_normalized_response_error=.01, max_relative_auc_error=.02,
                max_recovery_time_error=.30, max_raw_chemical_log_error=.01)


def verify(root):
    p=json.loads((root/'protocol.json').read_text())
    for f,h in {**p['source_sha256'],**p['input_sha256']}.items():
        if digest(f)!=h: raise ValueError('Changed input: '+f)
    return p


def prepare(root, baseline):
    root,baseline=Path(root).resolve(),Path(baseline).resolve()
    if root.exists(): raise FileExistsError(root)
    old=verify(baseline)
    if json.loads((baseline/'status.json').read_text())['state']!='completed':
        raise ValueError('Requires a completed pilot')
    root.mkdir(parents=True)
    sim=retime(old['checkpoint'],root/'refined-source.npz',old['dt']/2)
    p=copy.deepcopy(old)
    for filename in ('initial_states.npz','frozen_reference.npz'):
        shutil.copy2(baseline/filename,root/filename)
    p.update(checkpoint=str(root/'refined-source.npz'),dt=sim.config.dt,
             baseline=str(baseline),refinement_criteria=CRITERIA,
             time_refinement='Same-state timestep halving; assessment in refinement.json after completion.')
    p['source_config']['dt']=sim.config.dt
    from dataclasses import asdict
    p['source_config']=asdict(sim.config)
    inputs=[baseline/'protocol.json',baseline/'status.json',baseline/'comparison.json',
            root/'refined-source.npz',root/'initial_states.npz',root/'frozen_reference.npz']
    for job in p['jobs']:
        inputs.extend(baseline/moving.name(job)/f for f in ('history.json','result.json'))
    p['input_sha256'].update({str(f):digest(f) for f in inputs})
    for f in (Path(__file__),Path(__file__).with_name('feedback_survival_validation.py')):
        p['source_sha256'][str(f.resolve())]=digest(f)
    write_json(root/'protocol.json',p)
    write_json(root/'status.json',dict(state='prepared',protocol_sha256=digest(root/'protocol.json'),completed=0,total=len(p['jobs'])))
    return p


def response_error(coarse, coarse_control, fine, fine_control, factor):
    """Compare pulse-minus-control responses, normalized by initial log pulse."""
    arrays=[np.asarray(x,float) for x in (coarse,coarse_control,fine,fine_control)]
    if len({x.shape for x in arrays})!=1 or any(not np.isfinite(x).all() or np.any(x<=0) for x in arrays):
        raise ValueError('Positive, finite, aligned trajectories required')
    if factor<=0 or factor==1: raise ValueError('Nonzero positive pulse required')
    a,b,c,d=arrays
    return float(np.max(np.abs(np.log(a/b)-np.log(c/d)))/abs(np.log(factor)))


def assess(root):
    root=Path(root);p=verify(root);base=Path(p['baseline']);criteria=p['refinement_criteria']
    ph=digest(root/'protocol.json');bh=digest(base/'protocol.json')
    data={};raw=[]
    for folder,hash_value in ((root,ph),(base,bh)):
        for job in p['jobs']:
            key=moving.name(job);r=json.loads((folder/key/'result.json').read_text())
            if not r['quality_pass'] or r['protocol_sha256']!=hash_value: raise ValueError('Invalid result')
            h=json.loads((folder/key/'history.json').read_text())
            if len(h)!=round(p['duration']/p['interval'])+1: raise ValueError('Incomplete history')
            if not np.allclose([x['elapsed'] for x in h],np.arange(len(h))*p['interval'],atol=1e-9,rtol=0): raise ValueError('Time mismatch')
            with np.load(root/'initial_states.npz') as initial:
                if any(not np.array_equal(x['ids'],initial['ids']) for x in h): raise ValueError('ID mismatch')
            data[(folder,key)]=np.array([x['chemistry'] for x in h])
    reports={f:{moving.name(r['job']):r for r in json.loads((f/'comparison.json').read_text())['trials']} for f in (root,base)}
    rows=[]
    for job in p['jobs']:
        key=moving.name(job)
        raw.append(float(np.max(abs(np.log(data[(root,key)]/data[(base,key)])))))
        if job['target'] is None: continue
        control=job['family']+'_control'
        error=response_error(data[(base,key)],data[(base,control)],data[(root,key)],data[(root,control)],job['factor'])
        a=reports[base][key]['moving'];b=reports[root][key]['moving']
        metric='target_activator_log_auc_per_log_pulse'
        auc=abs(b[metric]/a[metric]-1)
        recovery={}
        for metric in ('target_recovery_time','network_recovery_time'):
            x,y=a[metric],b[metric]
            recovery[metric]=0. if x is None and y is None else (None if x is None or y is None else abs(x-y))
        passed=error<=criteria['max_normalized_response_error'] and auc<=criteria['max_relative_auc_error'] and all(x is not None and x<=criteria['max_recovery_time_error']+1e-9 for x in recovery.values())
        rows.append(dict(job=job,normalized_response_error=error,relative_auc_error=auc,recovery_time_errors=recovery,passed=bool(passed)))
    result=dict(passed=bool(all(r['passed'] for r in rows) and max(raw)<=criteria['max_raw_chemical_log_error']),
                criteria=criteria,max_raw_chemical_log_error=max(raw),trials=rows,
                scope='One timestep halving on one mature history. Tests chemical responses; does not validate spatial convergence, tiny mechanical residuals, or autonomous cell identity.')
    write_json(root/'refinement.json',result)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['prepare','run','assess'])
    parser.add_argument('--output',type=Path,default=Path('outputs/cell-response-moving-refined'))
    parser.add_argument('--baseline',type=Path,default=Path('outputs/cell-response-moving'))
    parser.add_argument('--workers',type=int,choices=range(1,5),default=2)
    a=parser.parse_args()
    if a.command=='prepare': prepare(a.output,a.baseline)
    elif a.command=='assess': assess(a.output)
    else:
        moving.run(a.output,a.workers)
        if json.loads((a.output/'status.json').read_text())['state']!='completed': raise RuntimeError('Simulation failed')
        assess(a.output)


if __name__=='__main__': main()
