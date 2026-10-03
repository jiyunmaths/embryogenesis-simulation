"""Full-horizon opt-in kernel gate against existing control/pulse trajectories."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from .fast_mechanics import FastAttributeSimulation,kernel
from .cell_response import pulse,response_metrics
from .cell_response_moving import observe
from .feedback_long import digest,save_checkpoint
from .resolution import write_json,_steps


def prepare(root):
    root=Path(root).resolve()
    if root.exists():raise FileExistsError(root)
    source=Path('outputs/cell-exchange-response-moving/unexchanged').resolve()
    p=json.loads((source/'protocol.json').read_text())
    for f,h in {**p['source_sha256'],**p['input_sha256']}.items():
        if digest(f)!=h:raise ValueError('Changed original input')
    jobs=[dict(key='unexchanged_control',cell=None,factor=1.),dict(key='unexchanged_cell-27_factor-0.9',cell=27,factor=.9)]
    files=[source/'protocol.json',source/'source.npz',source/'initial_states.npz']
    for job in jobs:
        folder=source/job['key'];r=json.loads((folder/'result.json').read_text())
        if not r['quality_pass'] or r['protocol_sha256']!=digest(source/'protocol.json'):raise ValueError('Requires accepted originals')
        files.extend(folder/f for f in ('result.json','history.json','latest_state.npz'))
    code=[Path(__file__),*[Path(__file__).with_name(f) for f in ('fast_mechanics.py','fast_mechanics.c','attribute_development.py','model.py','polarity.py','signaling.py','transport.py','cell_response.py','cell_response_moving.py','feedback_long.py','resolution.py')]]
    root.mkdir(parents=True)
    protocol=dict(source=str(source),checkpoint=str(source/'source.npz'),start=p['start'],duration=p['duration'],interval=p['interval'],dt=p['dt'],jobs=jobs,
        criteria=dict(chemical_log_max=1e-5,polarity_abs_max=1e-5,relative_axis_max=1e-4,relative_volume_max=1e-5,final_phi_abs_max=2e-5,normalized_response_max=.001,relative_auc_max=.002,recovery_time_error_max=.15),
        scope='Two complete 60-unit mature moving continuations, compared to stored original-kernel results at every recorded sample. Production runs unchanged. This checks one control and one negative pulse; not full developmental or all-regime equivalence.',
        source_sha256={str(f.resolve()):digest(f) for f in code},input_sha256={str(f):digest(f) for f in files})
    write_json(root/'protocol.json',protocol);write_json(root/'status.json',dict(state='prepared',completed=0,total=2))


def run(root):
    root=Path(root);p=json.loads((root/'protocol.json').read_text());ph=digest(root/'protocol.json');source=Path(p['source']);criteria=p['criteria'];kernel()
    for f,h in {**p['source_sha256'],**p['input_sha256']}.items():
        if digest(f)!=h:raise ValueError('Changed input '+f)
    histories={};done=0
    try:
        for job in p['jobs']:
            folder=root/job['key'];folder.mkdir(exist_ok=True);ck=folder/'latest_state.npz';reference=json.loads((source/job['key']/'history.json').read_text())
            write_json(root/'status.json',dict(state='running',completed=done,total=2,current=job['key']))
            if (folder/'result.json').exists():
                result=json.loads((folder/'result.json').read_text())
                if not result['passed'] or result['protocol_sha256']!=ph:raise ValueError('Saved result invalid')
                histories[job['key']]=json.loads((folder/'history.json').read_text());done+=1;continue
            if ck.exists():
                sim=FastAttributeSimulation.restore(ck)
                with np.load(ck) as d:meta=json.loads(str(d['long_experiment']))
                if meta['protocol_hash']!=ph or meta['job']!=job:raise ValueError('Restart mismatch')
                history=meta['history'];audit=meta['audit']
            else:
                sim=FastAttributeSimulation.restore(p['checkpoint'])
                if job['cell'] is not None:
                    x=pulse(np.array([sim.activator,sim.inhibitor]),list(sim.ids).index(job['cell']),job['factor']);sim.activator,sim.inhibitor=x
                history=[];audit=dict(chemical_log_max=0.,polarity_abs_max=0.,relative_axis_max=0.,relative_volume_max=0.,max_volume_error=0.,min_radius=1e100,max_clipping=0.,max_boundary=0.)
            step0=_steps(p['start'],p['dt']);stop=step0+_steps(p['duration'],p['dt']);every=_steps(p['interval'],p['dt']);last=history[-1]['time'] if history else -1
            while sim.step_number<=stop:
                v=sim.volumes();audit['max_volume_error']=max(audit['max_volume_error'],float(abs(v/sim.target-1).max()));audit['min_radius']=min(audit['min_radius'],float(((3*v/(4*np.pi))**(1/3)/sim.dx).min()));audit['max_clipping']=max(audit['max_clipping'],sim.clipped_fraction)
                if not np.isfinite(sim.activator).all() or not np.isfinite(sim.inhibitor).all() or np.any(sim.activator<=0) or np.any(sim.inhibitor<=0):raise RuntimeError('Chemical validity failure')
                if audit['max_volume_error']>=.05 or audit['min_radius']<4 or audit['max_clipping']>0:raise RuntimeError('Numerical quality failure')
                if (sim.step_number-step0)%every==0 and sim.time>last+1e-9:
                    row=observe(sim,sim.time-p['start']);ref=reference[(sim.step_number-step0)//every]
                    if row['ids']!=ref['ids'] or abs(row['time']-ref['time'])>1e-9:raise ValueError('Reference mismatch')
                    values=dict(chemical_log_max=float(abs(np.log(np.array(row['chemistry'])/ref['chemistry'])).max()),polarity_abs_max=float(abs(np.array(row['polarity'])-ref['polarity']).max()),relative_axis_max=abs(row['axis_ratio']/ref['axis_ratio']-1),relative_volume_max=float(abs(np.array(row['volumes'])/ref['volumes']-1).max()))
                    for k,val in values.items():audit[k]=max(audit[k],val)
                    audit['max_boundary']=max(audit['max_boundary'],row['boundary_occupancy'])
                    if any(audit[k]>criteria[k] for k in values) or audit['max_boundary']>=.01:raise RuntimeError('Reference agreement gate exceeded '+str(audit))
                    history.append(row);last=row['time'];write_json(folder/'history.json',history);write_json(folder/'status.json',dict(state='running',elapsed=row['elapsed'],audit=audit))
                    if (sim.step_number-step0)%_steps(3.,p['dt'])==0 or sim.step_number==stop:
                        save_checkpoint(sim,ck,audit,history,ph,job);print(job['key'],row['elapsed'],audit['chemical_log_max'],flush=True)
                if sim.step_number==stop:break
                sim.step()
            with np.load(source/job['key']/'latest_state.npz') as d:phi_error=float(abs(sim.phi-d['phi']).max())
            if phi_error>criteria['final_phi_abs_max']:raise RuntimeError('Endpoint field mismatch')
            write_json(folder/'result.json',dict(passed=True,protocol_sha256=ph,audit=audit,final_phi_abs_max=phi_error));histories[job['key']]=history;done+=1
            write_json(folder/'status.json',dict(state='completed',elapsed=p['duration'],audit=audit))
        def evaluate(collection):
            control,pulsed=[np.array([r['chemistry'] for r in collection[j['key']]]) for j in p['jobs']]
            times=np.array([r['elapsed'] for r in collection[p['jobs'][0]['key']]])
            with np.load(source/'initial_states.npz') as d:m=d['masses'];ids=d['ids']
            metrics=response_metrics(pulsed,control,times,m,list(ids).index(p['jobs'][1]['cell']),p['jobs'][1]['factor'])
            return np.log(pulsed/control)/abs(np.log(p['jobs'][1]['factor'])),metrics
        originals={j['key']:json.loads((source/j['key']/'history.json').read_text()) for j in p['jobs']}
        x,a=evaluate(histories);y,b=evaluate(originals);error=float(abs(x-y).max());auc=abs(a['target_activator_log_auc_per_log_pulse']/b['target_activator_log_auc_per_log_pulse']-1)
        recovery={}
        for key in ('target_recovery_time','network_recovery_time'):
            recovery[key]=0. if a[key] is None and b[key] is None else (None if a[key] is None or b[key] is None else abs(a[key]-b[key]))
        passed=error<=criteria['normalized_response_max'] and auc<=criteria['relative_auc_max'] and all(v is not None and v<=criteria['recovery_time_error_max']+1e-9 for v in recovery.values())
        write_json(root/'comparison.json',dict(passed=bool(passed),normalized_response_error=error,relative_auc_error=auc,recovery_time_errors=recovery,optimized=a,original=b,scope=p['scope']))
        write_json(root/'status.json',dict(state='completed',passed=bool(passed),completed=2,total=2))
    except Exception as exc:
        write_json(root/'status.json',dict(state='failed',completed=done,total=2,error=str(exc)));raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('command',choices=('prepare','run'));parser.add_argument('--output',type=Path,default=Path('outputs/fast-mechanics-validation'));a=parser.parse_args()
    prepare(a.output) if a.command=='prepare' else run(a.output)
