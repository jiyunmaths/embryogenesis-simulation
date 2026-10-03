"""Audited optimized-backend continuation of the moving response-transfer study.

Original scientific protocols stay immutable. A separate transition manifest
records accepted backend validation, code hashes, and restart checkpoints.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from . import cell_response_moving as original
from . import cell_exchange_response_moving as study
from .fast_mechanics import FastAttributeSimulation
from .feedback_long import digest,save_checkpoint,restore_checkpoint as original_restore
from .resolution import write_json,_steps
from .cell_response_moving import name,observe


def initialize(*args):
    sim=original.initialize(*args)
    sim.__class__=FastAttributeSimulation
    return sim


def restore_checkpoint(*args):
    sim,audit,history=original_restore(*args)
    sim.__class__=FastAttributeSimulation
    return sim,audit,history


def worker(args):
    root,job=args;root=Path(root);p=json.loads((root/'protocol.json').read_text());ph=digest(root/'protocol.json')
    folder=root/name(job);folder.mkdir(exist_ok=True)
    if (folder/'result.json').exists():
        result=json.loads((folder/'result.json').read_text())
        if result['protocol_sha256']!=ph:raise ValueError('Result protocol mismatch')
        return result
    checkpoint=folder/'latest_state.npz'
    if checkpoint.exists():sim,audit,history=restore_checkpoint(checkpoint,ph,job)
    else:
        with np.load(root/'initial_states.npz') as d:initial=d[job['family']]
        sim=initialize(p['checkpoint'],initial,job['target'],job['factor']);history=[]
        audit=dict(max_volume_error=0.,min_radius=1e100,max_clipping=0.,max_sampled_boundary=0.,elapsed_seconds=0.)
    step0=round(p['start']/p['dt']);stop=step0+_steps(p['duration'],p['dt']);every=_steps(p['interval'],p['dt'])
    checkpoint_every=_steps(3.,p['dt']);old=audit['elapsed_seconds'];started=time.monotonic()
    last=history[-1]['time'] if history else -1.
    try:
        while sim.step_number<=stop:
            v=sim.volumes();audit['max_volume_error']=max(audit['max_volume_error'],float(np.max(abs(v/sim.target-1))))
            audit['min_radius']=min(audit['min_radius'],float(np.min((3*v/(4*np.pi))**(1/3))/sim.dx))
            audit['max_clipping']=max(audit['max_clipping'],sim.clipped_fraction)
            if not np.isfinite(v).all() or not np.isfinite(sim.activator).all() or not np.isfinite(sim.inhibitor).all() or np.any(sim.activator<=0) or np.any(sim.inhibitor<=0):
                raise RuntimeError('Nonfinite or nonpositive state')
            if audit['max_volume_error']>=.05 or audit['min_radius']<4 or audit['max_clipping']>0:
                raise RuntimeError('Numerical quality limit exceeded')
            if (sim.step_number-step0)%every==0 and sim.time>last+1e-9:
                row=observe(sim,sim.time-p['start']);history.append(row);last=sim.time
                audit['max_sampled_boundary']=max(audit['max_sampled_boundary'],row['boundary_occupancy'])
                if audit['max_sampled_boundary']>=.01:raise RuntimeError('Boundary quality limit exceeded')
                audit['elapsed_seconds']=old+time.monotonic()-started
                write_json(folder/'history.json',history)
                write_json(folder/'status.json',dict(state='running',elapsed=row['elapsed'],until=p['duration'],audit=audit))
                if (sim.step_number-step0)%checkpoint_every==0 or sim.step_number==stop:
                    save_checkpoint(sim,checkpoint,audit,history,ph,job)
                    print(f'{name(job)}: elapsed model time {row["elapsed"]:.2f}; wall {audit["elapsed_seconds"]:.1f}s',flush=True)
            if sim.step_number==stop:break
            sim.step()
        result=dict(job=job,quality_pass=True,protocol_sha256=ph,audit=audit)
        write_json(folder/'result.json',result);write_json(folder/'status.json',dict(state='completed',elapsed=p['duration'],audit=audit))
        return result
    except Exception as exc:
        write_json(folder/'status.json',dict(state='failed',time=sim.time,error=str(exc),audit=audit));raise


def run(root):
    root=Path(root);p=json.loads((root/'protocol.json').read_text());study.verify(p)
    transition=root/'backend-transition.json';manifest=json.loads(transition.read_text())
    if digest(root/'protocol.json')!=manifest['original_protocol_sha256']:raise ValueError('Protocol changed')
    for f,h in {**manifest['source_sha256'],**manifest['validation_sha256']}.items():
        if digest(f)!=h:raise ValueError('Backend evidence changed: '+f)
    done=0
    try:
        for key in p['backgrounds']:
            while not study.source_ready(p['sources'][key]):
                write_json(root/'status.json',dict(state='waiting',completed=done,total=15,source=key,backend='optimized'));time.sleep(30)
            study.materialize(root,p,key)
            cp=json.loads((root/key/'protocol.json').read_text());study.verify(cp)
            for job in cp['jobs']:
                folder=root/key/name(job)
                write_json(root/'status.json',dict(state='running',completed=done,total=15,background=key,current=name(job),workers=1,backend='optimized',backend_transition_sha256=digest(transition)))
                if not (folder/'result.json').exists():
                    folder.mkdir(exist_ok=True)
                    execution=folder/'backend-execution.json'
                    if not execution.exists():
                        checkpoint=folder/'latest_state.npz'
                        write_json(execution,dict(backend_transition_sha256=digest(transition),resume_checkpoint_sha256=digest(checkpoint) if checkpoint.exists() else None,job=job))
                worker((str(root/key),job));done+=1
        study.assess(root,p)
        write_json(root/'status.json',dict(state='completed',completed=done,total=15,backend='optimized',backend_transition_sha256=digest(transition)))
    except Exception as exc:
        write_json(root/'status.json',dict(state='failed',completed=done,total=15,error=str(exc)));raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,default=Path('outputs/cell-exchange-response-moving'));a=parser.parse_args();run(a.output)
