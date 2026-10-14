"""Qualify faster mature clocks on actual fresh-history developmental endpoints.

Keep the accepted DevelopmentGpu solver and complete checkpoint schema. This
six-unit numerical screen neither resets chemistry nor launches long science.
"""
import argparse
from copy import deepcopy
from dataclasses import asdict, replace
from pathlib import Path
import fcntl
import json
import os
import time

import numpy as np
import torch

from . import fresh_zygote_pilot as pilot
from . import development_time_diagnosis as td
from . import qualify_repeated_development as q
from .model import Config
from .qualify_development_gpu import digest, read, write, pin, row, exact
from .resolution import _steps

ASSESSMENT=Path('docs/fresh_zygote_control_completion_assessment.json')
PARENT=Path('output/studies/fresh-zygote-control-completion')
ORIGINAL=Path('output/studies/fresh-zygote-pilot')
DESIGN=Path('docs/fresh_zygote_mature_qualification.md')
LAUNCH=Path('archive/study-launches/2026-10-12-fresh-zygote-mature-qualification')
PARENT_SHA='44af9b2bfdd77176891e015e5024db2aeca6acca413f7d8dd63e5b8a86ffd926'
HISTORIES=(1101,1102,1103)
BRANCHES=('off','tension_adhesion','polarity','full')
LEVELS=(('coarse',.00375),('fine',.001875))
START=18.;DURATION=6.;SNAPSHOTS=(0,3,6)
HOMOLOGOUS=dict(chemical_log_max=.01,polarity_abs_max=.01,
    relative_volume_max=.005,relative_transport_max=.01,unmatched_observations_max=0)


def candidate_config(base,dt):
    if base['seed'] not in HISTORIES or base['grid']!=120 or base['dt']!=.000234375:
        raise ValueError('Accepted primary fresh-history endpoint required')
    if dt not in dict(LEVELS).values():raise ValueError('Unregistered mature timestep')
    return asdict(replace(Config(**base),dt=dt,steps=_steps(START+DURATION,dt),save_every=_steps(.15,dt)))


def retime_checkpoint(source,destination,config):
    """Transform clocks only; preserve all array bytes and other metadata."""
    destination=Path(destination)
    if destination.exists():raise FileExistsError('Never overwrite a continuation start')
    with np.load(source,allow_pickle=False) as z:data={k:z[k].copy() for k in z.files}
    metadata=json.loads(str(data['metadata']));base=metadata['config']
    if (metadata['time']!=START or metadata['step_number']!=base['steps'] or
            metadata['divisions'] or len(data['ids'])!=16 or len(json.loads(str(data['event_audits'])))!=15 or
            str(data['chemical_method'])!=td.gpu.METHOD or int(data['development_schema'])!=td.gpu.SCHEMA):
        raise ValueError('Complete qualified t18 developmental endpoint required')
    if config!=candidate_config(base,config['dt']):raise ValueError('Only mature clocks may change')
    metadata['config']=config;metadata['step_number']=_steps(START,config['dt'])
    data['metadata']=np.array(json.dumps(metadata))
    np.savez_compressed(destination,**data)


def check_start(job):
    with np.load(job['source_checkpoint'],allow_pickle=False) as a,np.load(job['source'],allow_pickle=False) as b:
        if set(a.files)!=set(b.files):raise ValueError('Checkpoint fields changed during retiming')
        x,y=[json.loads(str(z['metadata'])) for z in (a,b)]
        base=x.pop('config');config=y.pop('config')
        if config!=job['config'] or config!=candidate_config(base,config['dt']):raise ValueError('Mature clock/model drift')
        if x.pop('step_number')!=base['steps'] or y.pop('step_number')!=_steps(START,config['dt']) or x!=y:
            raise ValueError('State, physical time, lineage or random metadata changed')
        for name in set(a.files)-{'metadata'}:
            if a[name].dtype!=b[name].dtype or a[name].shape!=b[name].shape or a[name].tobytes()!=b[name].tobytes():
                raise ValueError('Retiming changed saved array: '+name)


def design(root,assessment):
    jobs=[]
    for history in HISTORIES:
        for branch in BRANCHES:
            endpoint=next(e for e in assessment['combined_primary_endpoints'] if e['history']==history and e['branch']==branch)
            case=endpoint['case'];folder=(PARENT if history==1103 else ORIGINAL)/case
            with np.load(folder/'gpu-t18.npz',allow_pickle=False) as z:base=json.loads(str(z['metadata']))['config']
            for level,dt in LEVELS:
                key=f'h{history}_{branch}_{level}'
                jobs.append(dict(key=key,seed=history,branch=branch,level=level,reused=False,
                    source_checkpoint=str((folder/'gpu-t18.npz').resolve()),source_case=case,
                    source=str(Path(root)/key/'source.npz'),config=candidate_config(base,dt)))
    return jobs


def verify(p):
    q.verify(p)
    a=read(ASSESSMENT)
    if (not a['review_completed'] or not a['scoped_combined_developmental_preparations_accepted'] or
            digest(a['independent_review']['path'])!=a['independent_review']['sha256']):
        raise ValueError('Independently accepted fresh developmental preparations required')
    if (p['new_scientific_histories']!=0 or p['new_numerical_paths']!=24 or
            p['histories']!=list(HISTORIES) or p['branches']!=list(BRANCHES) or
            p['start_time']!=START or p['duration']!=DURATION or p['snapshots']!=list(SNAPSHOTS)):
        raise ValueError('Registered mature qualification scope changed')
    if p['jobs']!=design(Path(p['study_root']),a):raise ValueError('Registered endpoint design changed')
    parent=read(PARENT/'protocol.json')
    if digest(PARENT/'protocol.json')!=PARENT_SHA:raise ValueError('Accepted parent protocol changed')
    for key in ('physical_limits','time_limits','dense_time_spectrum_limit','prefix_steps','prefix_limits',
                'method','gpu_uuid','runtime_binary_sha256'):
        if p[key]!=parent[key]:raise ValueError('Inherited qualification contract changed: '+key)
    if (p['levels']!=[dict(level=l,dt=d) for l,d in LEVELS] or p['homologous_limits']!=HOMOLOGOUS or
            p['interval']!=.15 or p['diagnostic_interval']!=.015 or p['restart_elapsed']!=3. or p['restart_steps']!=8 or
            p['time_pairs']!=[dict(seed=h,branch=b,cases=[f'h{h}_{b}_{l}' for l,_ in LEVELS]) for h in HISTORIES for b in BRANCHES]):
        raise ValueError('Mature clock, pair or correspondence contract changed')
    for job in p['jobs']:check_start(job)


def prepare(root):
    root=Path(root).resolve()
    if root.exists():raise FileExistsError('Fresh qualification directory required')
    a=read(ASSESSMENT);review=read(a['independent_review']['path'])
    if (not a['review_completed'] or not a['scoped_combined_developmental_preparations_accepted'] or
            digest(a['independent_review']['path'])!=a['independent_review']['sha256'] or
            not review['scoped_combined_developmental_preparations_accepted']):
        raise ValueError('Independent combined endpoint acceptance required')
    parent=read(PARENT/'protocol.json')
    if digest(PARENT/'protocol.json')!=PARENT_SHA:raise ValueError('Accepted parent protocol changed')
    q.verify(parent)
    root.mkdir(parents=True);jobs=design(root,a);inputs=dict(parent['input_sha256'])
    inputs.update(pin([ASSESSMENT,a['independent_review']['path'],PARENT/'protocol.json',PARENT/'status.json',PARENT/'summary.json']))
    for job in jobs:
        folder=root/job['key'];folder.mkdir();retime_checkpoint(job['source_checkpoint'],job['source'],job['config']);check_start(job)
        inputs.update(pin([job['source_checkpoint'],job['source'],Path(job['source_checkpoint']).parent/'gpu-result.json']))
    sources=dict(parent['source_sha256'])
    sources.update(pin([Path(__file__),DESIGN,Path('tests/test_fresh_zygote_mature_qualification.py'),LAUNCH/'audit_prelaunch.py',LAUNCH/'launch.py']))
    p=dict(jobs=jobs,histories=list(HISTORIES),branches=list(BRANCHES),study_root=str(root),
        levels=[dict(level=l,dt=d) for l,d in LEVELS],new_scientific_histories=0,new_numerical_paths=24,
        start_time=START,duration=DURATION,interval=.15,diagnostic_interval=.015,snapshots=list(SNAPSHOTS),
        physical_limits=parent['physical_limits'],time_limits=parent['time_limits'],
        dense_time_spectrum_limit=parent['dense_time_spectrum_limit'],homologous_limits=HOMOLOGOUS,
        prefix_steps=8,prefix_limits=parent['prefix_limits'],restart_elapsed=3.,restart_steps=8,
        method=parent['method'],gpu_uuid=parent['gpu_uuid'],runtime_binary_sha256=parent['runtime_binary_sha256'],
        source_sha256=sources,input_sha256=inputs,parent_protocol_sha256=PARENT_SHA,
        time_pairs=[dict(seed=h,branch=b,cases=[f'h{h}_{b}_{l}' for l,_ in LEVELS]) for h in HISTORIES for b in BRANCHES],
        stop_policy='Stop on any prefix, physical, restart or pair failure; preserve failed and unexecuted cases. '
            'No automatic dt0.0009375 fallback, threshold change, rerun or 60/240-unit extension.',
        scope='Actual mature endpoints from three existing fresh developmental histories, four matched contexts, '
            'two larger timesteps, six elapsed units. Numerical clock/backend/transport qualification only; '
            'no new histories, reset chemistry, backend handoff, further division or long-window identity claim.')
    write(root/'protocol.json',p);verify(p)
    write(root/'status.json',dict(state='prepared',protocol_sha256=digest(root/'protocol.json')))
    return p


def invariant(sim):
    return deepcopy(dict(ids=sim.ids.tolist(),parents=sim.host.parents.tolist(),due=sim.host.due.tolist(),
        next_id=sim.host.next_id,lineage=sim.host.lineage,events=sim.event_audits,
        random_states={name:getattr(sim.host,name).bit_generator.state for name in ('rng','fate_rng','signal_rng')}))


def run_path(root,job,p):
    folder=Path(root)/job['key'];output=folder/'gpu-result.json'
    if output.exists():raise FileExistsError('Never rerun a recorded mature path')
    sim=td.gpu.DevelopmentGpu.restore(job['source']);start_step=sim.step_number
    initial=invariant(sim);initial_target=float(sim.target.sum());history=[];checkpoints=[];failure=None
    maxima=dict(max_volume_error=0.,min_radius=1e100,max_clipping=0.,dilution_amount_error=0.)
    clone=None;remaining=0;restart_passed=False;started=time.perf_counter();trace_rows=0
    try:
        with (folder/'network.jsonl').open('x') as network:
            for k in range(_steps(DURATION,job['config']['dt'])+1):
                if sim.step_number!=start_step+k or sim.time!=(start_step+k)*sim.config.dt:
                    raise ValueError('Mature continuation clock drift')
                now=sim.audit()
                for name in maxima:maxima[name]=min(maxima[name],now[name]) if name=='min_radius' else max(maxima[name],now[name])
                if not q.path_quality(maxima,p['physical_limits']):raise FloatingPointError('Every-step mature physical gate failed')
                u=torch.stack((sim.activator,sim.inhibitor))
                if not bool(torch.isfinite(u).all()) or float(u.min())<=0:raise FloatingPointError('Mature chemical positivity failed')
                if sim.host.divisions or len(sim.ids)!=16 or len(sim.event_audits)!=15:
                    raise ValueError('Unexpected division during mature continuation')
                if abs(float(sim.target.sum())/initial_target-1)>2e-14:raise ValueError('Mature target inheritance drift')
                if k%_steps(p['diagnostic_interval'],sim.config.dt)==0:
                    td.json_line(network,td.network_trace(sim));trace_rows+=1
                if k%sim.config.save_every==0:
                    if invariant(sim)!=initial:raise ValueError('Lineage, cycle or random state changed without division')
                    observed=row(sim,True);observed.update(elapsed=k*sim.config.dt,parents=sim.host.parents.tolist())
                    observed['spectrum']=q.spectrum(observed['volumes'],observed['delta']).tolist();history.append(observed)
                    if observed['boundary_occupancy']>=p['physical_limits']['boundary_max'] or any(v!=1 for v in observed['nondividing_components']):
                        raise FloatingPointError('Mature boundary/connectivity gate failed')
                    write(folder/'gpu-history.json',history);network.flush()
                    write(folder/'gpu-progress.json',dict(state='running',time=sim.time,elapsed=k*sim.config.dt,cells=16,
                        seconds=time.perf_counter()-started,audit=maxima))
                if k in [_steps(x,sim.config.dt) for x in SNAPSHOTS]:
                    elapsed=k*sim.config.dt;cp=folder/f'gpu-e{int(elapsed):02d}.npz';sim.checkpoint(cp);checkpoints.append(cp)
                if k==_steps(p['restart_elapsed'],sim.config.dt):
                    cp=folder/'restart.npz';sim.checkpoint(cp);clone=td.gpu.DevelopmentGpu.restore(cp);exact(sim,clone);remaining=p['restart_steps']
                if k==_steps(DURATION,sim.config.dt):break
                sim.step()
                if clone is not None:
                    clone.step();exact(sim,clone);remaining-=1
                    if remaining==0:clone=None;restart_passed=True
        if not restart_passed or sim.time!=START+DURATION or invariant(sim)!=initial:
            raise ValueError('Mature horizon or exact restart incomplete')
    except Exception as exc:
        failure=dict(type=type(exc).__name__,message=str(exc),time=sim.time,elapsed=sim.time-START)
        sim.checkpoint(folder/'failure-state.npz')
    write(folder/'gpu-history.json',history)
    report=dict(passed=failure is None,job=job,protocol_sha256=digest(Path(root)/'protocol.json'),failure=failure,
        final_time=sim.time,final_cells=len(sim.ids),events=sim.event_audits,audit=maxima,
        exact_restart_steps=8 if restart_passed else 0,lineage_cycle_rng_unchanged=invariant(sim)==initial,
        network_trace_rows=trace_rows,seconds=time.perf_counter()-started,
        evidence_sha256=pin([folder/'gpu-history.json',folder/'network.jsonl',*checkpoints,*folder.glob('*state.npz'),*folder.glob('restart.npz')]))
    write(output,report);write(folder/'gpu-progress.json',dict(state='completed' if report['passed'] else 'failed',
        time=sim.time,elapsed=sim.time-START,cells=len(sim.ids),failure=failure,seconds=report['seconds']))
    del sim,clone;torch.cuda.empty_cache();return report


def compare_pair(root,pair,p):
    root=Path(root);output=root/f"time-h{pair['seed']}-{pair['branch']}.json"
    if output.exists():raise FileExistsError('Never overwrite a mature comparison')
    a,b=pair['cases'];history=[read(root/key/'gpu-history.json') for key in (a,b)]
    errors,unmatched=q.compare_observations(*history);homologous,missing=q.homologous_errors(*history)
    if unmatched or missing:raise ValueError('Mature cell correspondence must remain exact')
    fields=[]
    for elapsed in SNAPSHOTS:
        densities=[q.checkpoint_density(root/key/f'gpu-e{elapsed:02d}.npz',None,3) for key in (a,b)]
        fields.append(dict(elapsed=elapsed,time=START+elapsed,errors=q.density_errors(densities[0][0],densities[1][0]),phase_ranges=[x[1] for x in densities]))
    field_max=max(v for item in fields for v in item['errors'].values())
    streams=[]
    for key in (a,b):
        with (root/key/'network.jsonl').open() as stream:streams.append([json.loads(x) for x in stream])
    dense=[];topology=[]
    if len(streams[0])!=401 or len(streams[1])!=401:raise ValueError('Incomplete dense mature record')
    for x,y in zip(*streams):
        if abs(x['time']-y['time'])>1e-9 or x['ids']!=y['ids'] or x['parents']!=y['parents']:
            raise ValueError('Dense mature correspondence changed')
        sa,sb=np.array(x['spectrum']),np.array(y['spectrum'])
        dense.append(dict(time=x['time'],elapsed=x['time']-START,relative_spectrum=float(np.linalg.norm(sa-sb)/(np.linalg.norm(sb) or 1.)),edges=[x['edges'],y['edges']]))
        if not np.array_equal(np.array(x['conductance'])>0,np.array(y['conductance'])>0):topology.append(x['time'])
    peak=max(dense,key=lambda x:x['relative_spectrum'])
    failed={k:v for k,v in errors.items() if v>p['time_limits'][k]}
    failed.update({'homologous_'+k:v for k,v in homologous.items() if v>p['homologous_limits'][k]})
    if field_max>p['time_limits']['field_max']:failed['field_max']=field_max
    if peak['relative_spectrum']>p['dense_time_spectrum_limit']:failed['dense_spectrum']=peak['relative_spectrum']
    report=dict(pair,passed=not failed,errors=errors,homologous_errors=homologous,fields=fields,field_max=field_max,
        dense_peak=peak,topology_mismatch_times=topology,failed_gates=failed,
        protocol_sha256=digest(root/'protocol.json'),evidence_sha256=pin(root/key/'gpu-result.json' for key in (a,b)))
    write(output,report);return report


def run(root):
    root=Path(root).resolve()
    with (root/'coordinator.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);p=read(root/'protocol.json');verify(p)
        if read(root/'status.json')!=dict(state='prepared',protocol_sha256=digest(root/'protocol.json')):
            raise ValueError('Unchanged prepared state required')
        if os.environ.get('CUDA_VISIBLE_DEVICES')!=p['gpu_uuid']:raise ValueError('Dedicated GPU selection required')
        torch.set_num_threads(1);lib=td.gpu.library()
        if digest(lib.binary_path)!=p['runtime_binary_sha256']:raise ValueError('Accepted runtime binary changed')
        write(root/'runtime.json',dict(gpu=torch.cuda.get_device_name(0),torch=torch.__version__,binary=lib.binary_path,binary_sha256=digest(lib.binary_path)))
        started=time.perf_counter();paths=[];pairs=[];prefixes=[];active=None
        def status(state,**extra):
            write(root/'status.json',dict(state=state,active=active,completed=len(paths),total=24,
                pairs_completed=len(pairs),pairs_total=12,new_scientific_histories=0,
                qualification_accepted=False,long_scientific_continuation_started=False,
                seconds=time.perf_counter()-started,protocol_sha256=digest(root/'protocol.json'),**extra))
        try:
            for pair in p['time_pairs']:
                for key in pair['cases']:
                    active=key;job=next(j for j in p['jobs'] if j['key']==key);status('running_prefix')
                    prefix=pilot.prefix(root,job,p);prefixes.append(prefix)
                    if not prefix['passed']:raise RuntimeError('Mature actual-start prefix failed: '+key)
                    status('running_mature');result=run_path(root,job,p);paths.append(dict(case=key,passed=result['passed'],seconds=result['seconds']))
                    if not result['passed']:raise RuntimeError('Mature path failed: '+key)
                status('comparing');result=compare_pair(root,pair,p);pairs.append(result)
                if not result['passed']:raise RuntimeError('Mature timestep gate failed: '+pair['branch'])
            verify(p)
            write(root/'summary.json',dict(execution_completed=True,mature_screens_passed=True,
                qualification_accepted=False,independent_review_pending=True,new_scientific_histories=0,
                new_numerical_paths=24,long_scientific_continuation_started=False,paths=paths,pairs=pairs,prefixes=prefixes,
                seconds=time.perf_counter()-started,protocol_sha256=digest(root/'protocol.json'),scope=p['scope']))
            status('completed_provisional',mature_screens_passed=True,independent_review_pending=True)
        except Exception as exc:
            status('failed',failure=dict(type=type(exc).__name__,message=str(exc)),paths=paths,pairs=pairs,prefixes=prefixes);raise


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=('prepare','run'));parser.add_argument('root',type=Path)
    args=parser.parse_args();(prepare if args.action=='prepare' else run)(args.root)


if __name__=='__main__':main()
