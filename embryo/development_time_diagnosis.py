"""Unchanged developmental solver, feedback-off timestep ladder and read-only traces."""
import argparse
from contextlib import contextmanager
from dataclasses import asdict, replace
from datetime import datetime, timezone
from pathlib import Path
import fcntl
import json
import os
import threading
import time
import numpy as np
import torch

from . import development_gpu as gpu
from .qualify_development_gpu import digest, read, write, pin, save_cpu, row, exact
from .qualify_repeated_development import (initial, spectrum, verify, event_passed,
    path_quality, comparison, TIME_LIMITS, PHYSICAL, FIELDS)
from .compensated_chemistry import NumpyStepper
from .compensated_moving import numpy_scale
from .model import Config
from .resolution import _steps

PARENT = Path('output/studies/repeated-development-gpu-qualification-v2')
PRESERVED = Path('archive/study-document-snapshots/2026-10-10-before-development-time-diagnosis/manifest.json')
LEVELS = (('coarse', .00375), ('fine', .001875), ('half', .0009375), ('quarter', .00046875))
_PROJECTION_LOCK = threading.Lock()


def json_line(stream, value):
    stream.write(json.dumps(value, allow_nan=False, separators=(',', ':'))+'\n')


@contextmanager
def projection_trace(sim, stream, extrema):
    """Forward original pressure routine and return its original tensors exactly.

    No arithmetic inputs, outputs or solver files are changed. Single-threaded
    instrumentation is guarded, and restored even when the original raises.
    """
    if not _PROJECTION_LOCK.acquire(blocking=False):
        raise RuntimeError('Projection tracing is serial only')
    original = gpu.project_phase
    parents = [int(cell) for cell in sim.ids if int(cell) in sim.host.divisions]
    called = []

    def traced(field, carry, volume, dx):
        parent = parents[len(called)]
        before = float(gpu.occupancy64(field.double()+carry).sum()*dx**3)
        result = original(field, carry, volume, dx)
        after = float(gpu.occupancy64(result[0].double()+result[1]).sum()*dx**3)
        item = dict(step=sim.step_number, time_start=sim.time, time_end=sim.time+sim.config.dt,
            parent=parent, capacity=volume, before_relative_volume_error=before/volume-1,
            after_relative_volume_error=after/volume-1, correction_magnitude=result[2])
        called.append(parent)
        extrema['correction_magnitude_max'] = max(extrema['correction_magnitude_max'], result[2])
        extrema['before_relative_volume_error_max'] = max(extrema['before_relative_volume_error_max'], abs(before/volume-1))
        extrema['after_relative_volume_error_max'] = max(extrema['after_relative_volume_error_max'], abs(after/volume-1))
        extrema['calls'] += 1
        json_line(stream, item)
        return result

    gpu.project_phase = traced
    try:
        yield
        assert called == parents
    finally:
        gpu.project_phase = original
        _PROJECTION_LOCK.release()


def network_trace(sim):
    """Only temporary arrays and host scalars; no simulation state writes."""
    contact, delta, _, conductance = [x.cpu().numpy() for x in sim.matrices()]
    volumes = sim.geometry[:, 0].cpu().numpy()
    active = []
    for parent, ev in sim.host.divisions.items():
        i = int(np.flatnonzero(sim.ids == parent)[0])
        _, _, neck, overlap, fraction, _ = sim._cleavage(i, ev)
        progress = min(1., max(0., (sim.time-ev['start'])/sim.config.cytokinesis_duration))
        active.append(dict(parent=int(parent), start=ev['start'], progress=progress,
            smooth_ramp=progress**2*(3-2*progress), neck=neck, overlap=overlap,
            fraction=fraction, ready_by_duration=sim.time-ev['start']>=sim.config.cytokinesis_duration,
            neck_margin=sim.config.neck_threshold-neck,
            overlap_margin=sim.config.division_overlap_tolerance-overlap))
    threshold = sim.config.graph_contact_cutoff * float(contact.max())
    n = len(volumes); upper = np.triu_indices(n, 1)
    return dict(time=sim.time, step=sim.step_number, ids=sim.ids.tolist(), parents=sim.host.parents.tolist(),
        volumes=volumes.tolist(), centers=sim.geometry[:, 1:4].cpu().numpy().tolist(),
        contact=contact.tolist(), conductance=conductance.tolist(), delta=delta.tolist(),
        spectrum=spectrum(volumes,delta).tolist(), cutoff=threshold,
        contact_cutoff_margins=(contact[upper]-threshold).tolist(),
        edges=int(np.count_nonzero(conductance[upper])), active_divisions=active,
        completed_events=len(sim.event_audits), last_step_projection_max=sim.volume_projection_max)


def prepare(root):
    root = Path(root).resolve()
    if root.exists(): raise FileExistsError('Fresh study directory required')
    parent = read(PARENT/'protocol.json'); summary = read(PARENT/'summary.json')
    assessment = read('docs/repeated_development_gpu_assessment.json')
    review_path = Path(assessment['independent_review'])
    if digest(review_path) != assessment['independent_review_sha256']: raise ValueError('Parent review changed')
    review = read(review_path)
    assert review['review_completed'] and not review['qualification_accepted'] and not summary['passed']
    assert review['protocol_sha256'] == digest(PARENT/'protocol.json')
    preserved = read(PRESERVED)
    for file, expected in preserved['preserved_sha256'].items():
        if digest(file) != expected: raise ValueError('Preserved evidence changed: '+file)
    root.mkdir(parents=True); jobs=[]
    inputs=[PARENT/'protocol.json',PARENT/'summary.json',review_path,
        Path('docs/repeated_development_gpu_assessment.json'),PRESERVED]
    base=next(j for j in parent['jobs'] if j['key']=='off_coarse')['config']
    for level,dt in LEVELS:
        config=replace(Config(**base),dt=dt,steps=_steps(18.,dt),save_every=_steps(.15,dt))
        host=initial(config);key='off_'+level;folder=root/key;folder.mkdir();source=folder/'source.npz'
        save_cpu(host,source);inputs.append(source)
        jobs.append(dict(key=key,level=level,branch='off',source=str(source),config=asdict(config)))
        if level in ('coarse','fine'):
            inputs.extend(PARENT/key/p for p in ('source.npz','gpu-history.json','gpu-result.json','gpu-lineage.json','gpu-events.json'))
            inputs.extend(PARENT/key/f'gpu-t{s:02d}.npz' for s in (0,3,6,9,12,15,18))
    inputs.extend(PARENT/'off_coarse'/p for p in ('cpu-result.json','cpu-history.json'))
    inputs.extend(PARENT/'off_fine'/f'restart-stage-{stage}.npz' for stage in (2,4,8))
    sources=dict(parent['source_sha256'])
    sources.update(pin([Path(__file__),Path('docs/development_time_diagnosis.md'),Path('tests/test_development_time_diagnosis.py')]))
    protocol=dict(jobs=jobs,duration=18.,interval=.15,checkpoint_interval=3.,diagnostic_interval=.015,
        time_limits=dict(TIME_LIMITS),physical_limits=dict(PHYSICAL),restart_steps=8,restart_stages=[2,4,8],
        minimum_mature_window=3.,diagnostic_gate_steps=8,old_trajectory_reproduction='exact history, events, lineage and every saved checkpoint array',
        source_sha256=sources,input_sha256=pin(inputs),method=gpu.METHOD,gpu_uuid=parent['gpu_uuid'],
        new_scientific_histories=0,preserved_manifest=str(PRESERVED.resolve()),
        parent_protocol_sha256=digest(PARENT/'protocol.json'),parent_study=str(PARENT.resolve()),
        scope='Feedback-off numerical diagnosis only. Four time levels at grid 72; old two are diagnostic reruns, two finer levels are new numerical paths. No new biology, model changes, full CPU fine ladder or spatial acceptance.')
    write(root/'protocol.json',protocol)
    write(root/'status.json',dict(state='prepared',protocol_sha256=digest(root/'protocol.json')))
    return protocol


def diagnostic_gate(root):
    root=Path(root);p=read(root/'protocol.json');verify(p)
    if os.environ.get('CUDA_VISIBLE_DEVICES')!=p['gpu_uuid']:raise ValueError('Explicit designated GPU required')
    target=root/'instrumentation-gate.json'
    if target.exists():raise FileExistsError('Never overwrite an instrumentation gate')
    torch.set_num_threads(1);records=[]
    for stage in p['restart_stages']:
        source=PARENT/'off_fine'/f'restart-stage-{stage}.npz'
        a=gpu.DevelopmentGpu.restore(source);b=gpu.DevelopmentGpu.restore(source);exact(a,b)
        extrema=dict(calls=0,correction_magnitude_max=0.,before_relative_volume_error_max=0.,after_relative_volume_error_max=0.)
        with (root/f'gate-stage-{stage}-projection.jsonl').open('w') as stream:
            for _ in range(p['diagnostic_gate_steps']):
                network_trace(b);exact(a,b)
                a.step()
                with projection_trace(b,stream,extrema):b.step()
                network_trace(b);exact(a,b)
        records.append(dict(stage=stage,steps=p['diagnostic_gate_steps'],bitwise_equal=True,projection=extrema))
        del a,b;torch.cuda.empty_cache()
    result=dict(passed=True,rows=records,protocol_sha256=digest(root/'protocol.json'),
        evidence_sha256=pin(root.glob('gate-stage-*-projection.jsonl')))
    write(target,result);return result


def reproduce(folder, parent_folder):
    """Instrumentation cannot reinterpret an old failure as a new trajectory."""
    for name in ('gpu-history.json','gpu-lineage.json','gpu-events.json'):
        if read(folder/name)!=read(parent_folder/name):raise FloatingPointError('Instrumented trajectory changed: '+name)
    for stamp in (0,3,6,9,12,15,18):
        name=f'gpu-t{stamp:02d}.npz'
        with np.load(folder/name,allow_pickle=False) as a,np.load(parent_folder/name,allow_pickle=False) as b:
            assert set(a.files)==set(b.files)
            for key in a.files:
                if not np.array_equal(a[key],b[key]):raise FloatingPointError('Instrumented checkpoint changed: '+name+':'+key)
    return dict(passed=True,checkpoints=7,exact_history_events_lineage=True)


def run_path(root,job):
    root=Path(root);p=read(root/'protocol.json');verify(p);folder=root/job['key']
    output=folder/'gpu-result.json'
    if output.exists():raise FileExistsError('Never rerun a recorded case')
    sim=gpu.DevelopmentGpu.restore(job['source']);clone=None;remaining=0;stage=None
    history=[];checkpoints=[];restart={};observed=0;initial_target=float(sim.target.sum());independent=0.
    extrema=dict(calls=0,correction_magnitude_max=0.,before_relative_volume_error_max=0.,after_relative_volume_error_max=0.)
    audit=dict(max_volume_error=0.,min_radius=1e100,max_clipping=0.,dilution_amount_error=0.)
    started=time.perf_counter();failure=None;reproduction=None;trace_rows=0;previous_active=set();previous_events=0
    try:
        with (folder/'projection.jsonl').open('w') as pressure,(folder/'network.jsonl').open('w') as network:
            for k in range(job['config']['steps']+1):
                now=sim.audit()
                for key,value in now.items():
                    if key in audit:audit[key]=min(audit[key],value) if key=='min_radius' else max(audit[key],value)
                if not path_quality(audit,p['physical_limits']):raise FloatingPointError('Every-step physical gate failed')
                u=torch.stack((sim.activator,sim.inhibitor)).cpu().numpy()
                if not np.isfinite(u).all() or u.min()<=0:raise FloatingPointError('Chemical positivity failed')
                if abs(float(sim.target.sum())/initial_target-1)>2e-14:raise FloatingPointError('Target inheritance drift')
                for ev in sim.event_audits[observed:]:
                    if not event_passed(ev,p['physical_limits']):raise FloatingPointError('Division event failed')
                observed=len(sim.event_audits);active=set(sim.host.divisions)
                if k%_steps(p['diagnostic_interval'],job['config']['dt'])==0 or active!=previous_active or observed!=previous_events:
                    json_line(network,network_trace(sim));trace_rows+=1
                previous_active=active;previous_events=observed
                if k%job['config']['save_every']==0:
                    record=row(sim,True);record['parents']=sim.host.parents.tolist()
                    record['spectrum']=spectrum(record['volumes'],record['delta']).tolist();history.append(record)
                    if record['boundary_occupancy']>=.01 or any(c!=1 for c in record['nondividing_components']):raise FloatingPointError('Boundary/connectivity gate failed')
                    write(folder/'gpu-history.json',history);pressure.flush();network.flush()
                    write(folder/'gpu-progress.json',dict(state='running',time=sim.time,cells=len(sim.ids),events=observed,
                        active_divisions=len(active),seconds=time.perf_counter()-started,audit=audit,projection=extrema))
                if k%_steps(p['checkpoint_interval'],job['config']['dt'])==0:
                    cp=folder/f'gpu-t{int(round(sim.time)):02d}.npz';sim.checkpoint(cp);checkpoints.append(cp)
                if clone is None:
                    for candidate in p['restart_stages']:
                        if candidate not in restart and candidate<len(sim.ids)<2*candidate and active:
                            stage=candidate;cp=folder/f'restart-stage-{stage}.npz';sim.checkpoint(cp)
                            clone=gpu.DevelopmentGpu.restore(cp);exact(sim,clone);remaining=p['restart_steps'];break
                if k==job['config']['steps']:break
                if k<8:
                    expected,correction=NumpyStepper(sim.matrices()[1].cpu().numpy(),job['config']['dt'],sim.config.signal_beta,
                        sim.config.signal_da,sim.config.signal_dh)(u,sim.chemical_correction.cpu().numpy())
                    old=sim.geometry[:,0].cpu().numpy().copy()
                with projection_trace(sim,pressure,extrema):sim.step()
                if k<8:
                    expected,correction=numpy_scale(expected,correction,old/sim.geometry[:,0].cpu().numpy())
                    actual=torch.stack((sim.activator,sim.inhibitor)).cpu().numpy()
                    independent=max(independent,float(abs(np.log(actual/expected)).max()))
                    if independent>1e-11:raise FloatingPointError('Independent chemical/dilution prefix failed')
                if clone is not None:
                    clone.step();exact(sim,clone);remaining-=1
                    if remaining==0:restart[stage]=p['restart_steps'];clone=None;stage=None
        if len(sim.ids)!=16 or len(sim.event_audits)!=15 or sim.host.divisions:raise RuntimeError('Unfinished development')
        if 18-max(e['time'] for e in sim.event_audits)<3:raise RuntimeError('Insufficient mature window')
        if set(restart)!=set(p['restart_stages']):raise RuntimeError('Missing exact restart window')
    except Exception as exc:
        failure=dict(type=type(exc).__name__,message=str(exc),time=sim.time,cells=len(sim.ids))
        sim.checkpoint(folder/'diagnostic-failure-state.npz')
    write(folder/'gpu-history.json',history);write(folder/'gpu-lineage.json',sim.host.lineage);write(folder/'gpu-events.json',sim.event_audits)
    if failure is None and job['level'] in ('coarse','fine'):
        try:reproduction=reproduce(folder,PARENT/job['key'])
        except Exception as exc:failure=dict(type=type(exc).__name__,message=str(exc),time=sim.time,cells=len(sim.ids))
    files=[folder/'gpu-history.json',folder/'gpu-lineage.json',folder/'gpu-events.json',folder/'projection.jsonl',folder/'network.jsonl',
        *checkpoints,*folder.glob('restart-stage-*.npz')]
    if (folder/'diagnostic-failure-state.npz').exists():files.append(folder/'diagnostic-failure-state.npz')
    result=dict(passed=failure is None,job=job,arm='gpu',failure=failure,audit=audit,events=sim.event_audits,
        exact_restart_steps={str(k):v for k,v in restart.items()},independent_chemical_steps=8,
        independent_chemical_log_max=independent,final_cells=len(sim.ids),final_time=sim.time,
        projection=extrema,network_trace_rows=trace_rows,old_trajectory_reproduction=reproduction,
        seconds=time.perf_counter()-started,protocol_sha256=digest(root/'protocol.json'),evidence_sha256=pin(files))
    write(output,result);write(folder/'gpu-progress.json',dict(state='completed' if result['passed'] else 'failed',
        time=sim.time,cells=len(sim.ids),events=len(sim.event_audits),seconds=result['seconds'],failure=failure,audit=audit))
    del sim,clone;torch.cuda.empty_cache();return result


def assess(root):
    root=Path(root);p=read(root/'protocol.json');verify(p);reports=[]
    for job in p['jobs']:
        result=read(root/job['key']/'gpu-result.json')
        assert result['protocol_sha256']==digest(root/'protocol.json')
        for file,expected in result['evidence_sha256'].items():
            if digest(file)!=expected:raise ValueError('Recorded evidence changed: '+file)
        reports.append(result)
    pairs=[comparison(root,p['jobs'][i]['key'],p['jobs'][i+1]['key'],p['time_limits']) for i in range(3)]
    trend=None
    if all('errors' in r for r in pairs):
        keys=('relative_spectrum_max','concentration_wasserstein_max','log_wasserstein_max','relative_axis_max')
        trend={key:[r['errors'][key] for r in pairs] for key in keys}
        trend.update({key:[max(f['errors'][key] for f in r['fields']) for r in pairs] for key in FIELDS})
    result=dict(execution_completed=True,all_paths_passed=all(r['passed'] for r in reports),
        finest_pair_time_screen_passed=all(r['passed'] for r in reports) and pairs[-1]['passed'],
        qualification_accepted=False,independent_review_pending=True,new_scientific_histories=0,
        protocol_sha256=digest(root/'protocol.json'),pairs=pairs,adjacent_error_sequences=trend,
        old_trajectory_reproduction=[r['old_trajectory_reproduction'] for r in reports[:2]],
        evidence_sha256=pin(root/job['key']/'gpu-result.json' for job in p['jobs']),
        scope=p['scope'])
    write(root/'summary.json',result);return result


def run(root):
    root=Path(root).resolve()
    with (root/'coordinator.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);p=read(root/'protocol.json');verify(p)
        status=read(root/'status.json');gate=read(root/'instrumentation-gate.json')
        assert status['state']=='prepared' and status['protocol_sha256']==digest(root/'protocol.json')
        if not gate['passed'] or gate['protocol_sha256']!=digest(root/'protocol.json'):raise ValueError('Instrumentation gate required')
        for file,expected in gate['evidence_sha256'].items():
            if digest(file)!=expected:raise ValueError('Instrumentation evidence changed')
        if os.environ.get('CUDA_VISIBLE_DEVICES')!=p['gpu_uuid']:raise ValueError('Explicit designated GPU required')
        torch.set_num_threads(1);started=time.perf_counter();completed=[]
        lib=gpu.library();write(root/'runtime.json',dict(gpu=torch.cuda.get_device_name(0),torch=torch.__version__,
            binary=lib.binary_path,binary_sha256=digest(lib.binary_path)))
        try:
            for job in p['jobs']:
                write(root/'status.json',dict(state='running',case=job['key'],completed=completed,total=4))
                result=run_path(root,job);completed.append(job['key'])
                print(job['key'],result['passed'],result['seconds'],result['failure'],flush=True)
                if not result['passed'] and job['level'] in ('coarse','fine'):
                    raise RuntimeError('Old trajectory reproduction/quality failed; do not run finer candidates')
            write(root/'status.json',dict(state='running',stage='fixed_time_assessment',completed=completed,total=4))
            result=assess(root);write(root/'timing.json',dict(seconds=time.perf_counter()-started))
            write(root/'status.json',dict(state='completed_diagnostic',finest_pair_time_screen_passed=result['finest_pair_time_screen_passed'],
                qualification_accepted=False,independent_review_pending=True,total=4));return result
        except Exception as exc:
            write(root/'status.json',dict(state='failed',error=str(exc)));raise


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('prepare','gate','run','assess'))
    parser.add_argument('--output',type=Path,default=Path('output/studies/development-time-diagnosis'))
    args=parser.parse_args()
    result={'prepare':prepare,'gate':diagnostic_gate,'run':run,'assess':assess}[args.action](args.output)
    print(json.dumps({k:result[k] for k in ('passed','finest_pair_time_screen_passed','qualification_accepted') if k in result}),flush=True)


if __name__=='__main__':main()
