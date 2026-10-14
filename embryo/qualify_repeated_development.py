"""Prospective naturally scheduled repeated-cleavage qualification, not science."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from multiprocessing import get_context
from pathlib import Path
import fcntl
import json
import os
import time

import numpy as np
from scipy.integrate import quad
from scipy.special import expit
from scipy.stats import wasserstein_distance
import torch

from .model import Config, occupancy
from .development_reference import ReferenceDevelopment
from .development_gpu import DevelopmentGpu, METHOD, library
from .qualify_development_gpu import (digest, read, write, pin, save_cpu,
    reference, row, exact, audit_cpu)
from .compensated_chemistry import NumpyStepper
from .compensated_moving import numpy_scale
from .cleavage_measurement import reconstruct, analytic
from .resolution import _steps

BRANCHES = ('off', 'tension_adhesion', 'polarity', 'full')
LEVELS = dict(coarse=(72, .00375), fine=(72, .001875),
    space88=(88, .001875), space104=(104, .001875))
FIELDS = ('occupancy', 'activator_density', 'inhibitor_density', 'polarity_density')
BACKEND_LIMITS = dict(concentration_wasserstein_max=1e-5, log_wasserstein_max=1e-5,
    contrast_abs_max=1e-5, relative_axis_max=1e-4, relative_total_volume_max=1e-5,
    whole_amount_relative_max=1e-5, mean_polarity_abs_max=1e-5,
    relative_spectrum_max=1e-5, cell_count_difference_max=1,
    unmatched_observations_max=30, field_max=2e-5, maturity_time_difference_max=.0075)
TIME_LIMITS = dict(concentration_wasserstein_max=.01, log_wasserstein_max=.01,
    contrast_abs_max=.005, relative_axis_max=.01, relative_total_volume_max=.005,
    whole_amount_relative_max=.005, mean_polarity_abs_max=.01,
    relative_spectrum_max=.01, cell_count_difference_max=1,
    unmatched_observations_max=30, field_max=.01, maturity_time_difference_max=.05)
SPACE_LIMITS = dict(TIME_LIMITS, relative_spectrum_max=.05,
    maturity_time_difference_max=.6)
PHYSICAL = dict(volume_error_max=.05, radius_min=4., clipping_max=0., boundary_max=.01,
    dilution_amount_error_max=2e-14, event_amount_jump_max=2e-14,
    event_target_jump_max=2e-14, event_volume_jump_max=1e-6,
    event_occupancy_jump_max=1e-14)
HOMOLOGOUS_LIMITS = dict(chemical_log_max=1e-5, polarity_abs_max=1e-5,
    relative_volume_max=1e-5, relative_transport_max=1e-5, phase_abs_max=2e-5,
    unmatched_observations_max=30)


def verify(protocol):
    for file, value in {**protocol['source_sha256'], **protocol['input_sha256']}.items():
        if digest(file) != value: raise ValueError('Qualification pin changed: ' + file)


def spectrum(volumes, delta):
    v, d = np.asarray(volumes), np.asarray(delta)
    if d.shape != (len(v), len(v)) or np.any(v <= 0): raise ValueError('Invalid transport arrays')
    mass = v[:, None] * d; scale = max(1., np.linalg.norm(d))
    if (not np.isfinite(d).all() or abs(d.sum(1)).max() > 1e-12 * scale or
        abs(v @ d).max() > 1e-12 * scale or abs(mass - mass.T).max() > 1e-12 * scale):
        raise FloatingPointError('Nonconservative or asymmetric weighted transport')
    off = d.copy(); np.fill_diagonal(off, 0)
    if off.min() < 0 or np.diag(d).max() > 0: raise FloatingPointError('Invalid exchange signs')
    weighted = -np.sqrt(v)[:, None] * d / np.sqrt(v)[None, :]
    return np.linalg.eigvalsh(.5 * (weighted + weighted.T))


def compare_observations(a, b):
    """Label-invariant comparisons at unchanged physical observation times."""
    keys = tuple(k for k in TIME_LIMITS if k not in ('field_max', 'maturity_time_difference_max'))
    error = dict.fromkeys(keys, 0.); unmatched = []
    if len(a) != len(b): raise ValueError('Incomplete observation histories')
    for x, y in zip(a, b):
        if abs(x['time'] - y['time']) > 1e-9: raise ValueError('Observation clocks differ')
        u, w = np.array(x['chemistry']), np.array(y['chemistry'])
        v, z = np.array(x['volumes']), np.array(y['volumes'])
        if min(u.min(), w.min(), v.min(), z.min()) <= 0: raise ValueError('Nonpositive observation')
        count_difference = abs(len(v) - len(z))
        values = dict(concentration_wasserstein_max=max(wasserstein_distance(i,j,v,z) for i,j in zip(u,w)),
            log_wasserstein_max=max(wasserstein_distance(np.log(i),np.log(j),v,z) for i,j in zip(u,w)),
            contrast_abs_max=float(abs(u.std(1)-w.std(1)).max()),
            relative_axis_max=abs(x['axis_ratio']/y['axis_ratio']-1),
            relative_total_volume_max=abs(v.sum()/z.sum()-1),
            whole_amount_relative_max=float(abs((u@v)/(w@z)-1).max()),
            mean_polarity_abs_max=float(abs(v@np.array(x['polarity'])/v.sum()-z@np.array(y['polarity'])/z.sum()).max()),
            cell_count_difference_max=count_difference)
        sx, sy = spectrum(v, x['delta']), spectrum(z, y['delta'])
        if count_difference:
            unmatched.append(x['time'])
        else:
            denominator = np.linalg.norm(sy)
            values['relative_spectrum_max'] = float(np.linalg.norm(sx-sy)/denominator) if denominator else float(np.linalg.norm(sx-sy))
        for key, value in values.items(): error[key] = max(error[key], float(value))
    error['unmatched_observations_max'] = len(unmatched)
    return error, unmatched


def densities(phi, chemistry, polarity):
    h = occupancy(np.asarray(phi, dtype=np.float64))
    u, p = np.asarray(chemistry), np.asarray(polarity)
    return np.concatenate((h.sum(0)[None], np.einsum('sn,nijk->sijk', u, h),
        np.einsum('nd,nijk->dijk', p, h)), axis=0)


def density_errors(a, b):
    """First three are relative L2; polarity uses occupancy as its denominator."""
    if a.shape != b.shape or a.shape[0] != 6: raise ValueError('Wrong density fields')
    norm = lambda x: float(np.sqrt(np.sum(x*x, dtype=np.float64)))
    errors = [norm(a[k]-b[k])/max(norm(b[k]), 1e-30) for k in range(3)]
    errors.append(norm(a[3:]-b[3:])/max(norm(b[0]), 1e-30))
    return dict(zip(FIELDS, errors))


def homologous_errors(a,b):
    """Only call after completed ancestry and partition correspondence passes."""
    errors=dict.fromkeys(('chemical_log_max','polarity_abs_max','relative_volume_max','relative_transport_max'),0.)
    unmatched=[]
    for x,y in zip(a,b):
        if x['ids']!=y['ids'] or x['parents']!=y['parents']:
            unmatched.append(x['time']);continue
        u,v=np.array(x['chemistry']),np.array(y['chemistry'])
        d,e=np.array(x['delta']),np.array(y['delta']);norm=np.linalg.norm(e)
        values=dict(chemical_log_max=float(abs(np.log(u/v)).max()),
            polarity_abs_max=float(abs(np.array(x['polarity'])-np.array(y['polarity'])).max()),
            relative_volume_max=float(abs(np.array(x['volumes'])/np.array(y['volumes'])-1).max()),
            relative_transport_max=float(np.linalg.norm(d-e)/norm) if norm else float(np.linalg.norm(d-e)))
        for key,value in values.items():errors[key]=max(errors[key],value)
    errors['unmatched_observations_max']=len(unmatched)
    return errors,unmatched


def initial(config):
    host = ReferenceDevelopment(config)
    radius = np.sqrt(np.sum(host.xyz**2, axis=0))
    host.phi = expit(np.sqrt(2)*(.8-radius)/config.interface_width)[None].astype(np.float32)
    host.phase_carry = np.zeros_like(host.phi, dtype=np.float64)
    target = quad(lambda r:4*np.pi*r*r*occupancy(expit(np.sqrt(2)*(.8-r)/config.interface_width)),
        0, .8+30*config.interface_width, epsabs=1e-12, epsrel=1e-12)[0]
    host.target[:] = target; host.initial_volume = target
    host.activator[:] = 1.03; host.inhibitor[:] = .98; host.polarity[:] = 0
    return host


def prepare(root):
    root = Path(root).resolve()
    if root.exists(): raise FileExistsError('Choose a fresh study directory')
    accepted = Path('docs/development_gpu_qualification_assessment.json')
    acceptance = read(accepted)
    if not acceptance['passed']: raise ValueError('Single-cleavage method not accepted')
    review = Path(acceptance['independent_review'])
    if digest(review) != acceptance['independent_review_sha256']: raise ValueError('Acceptance record changed')
    parent = read(review)
    for file, value in parent['consumed_sha256'].items():
        if digest(file) != value: raise ValueError('Accepted evidence changed: '+file)
    root.mkdir(parents=True); jobs=[]; inputs=[accepted, review]
    for branch in BRANCHES:
        for level, (grid, dt) in LEVELS.items():
            config = Config(grid=grid, extent=2.24, dt=dt, steps=_steps(18.,dt),
                save_every=_steps(.15,dt), max_cells=16, division_interval=2., signal_dh=.55,
                signal_partition_noise=.001, feedback=branch!='off',
                fate_tension=.25 if branch in ('tension_adhesion','full') else 0.,
                fate_adhesion=.35 if branch in ('tension_adhesion','full') else 0.,
                polarity_tension=.35 if branch in ('polarity','full') else 0.)
            host=initial(config); name=branch+'_'+level; folder=root/name; folder.mkdir()
            source=folder/'source.npz'; save_cpu(host,source); inputs.append(source)
            jobs.append(dict(key=name,branch=branch,level=level,source=str(source),config=asdict(host.config)))
    names=('qualify_repeated_development.py','qualify_development_gpu.py','development_gpu.py',
        'development_reference.py','gpu_development_kernels.cu','gpu_backend.py','gpu_precision_control.py',
        'gpu_precision_kernels.cu','gpu_kernels.cu','cuda_spatial_bench.cu','cuda_mechanics.cu',
        'compensated_moving.py','compensated_chemistry.py','attribute_development.py','model.py',
        'transport.py','signaling.py','polarity.py','cell_response_moving.py','shape.py','cleavage_measurement.py','resolution.py')
    sources=[Path(__file__).with_name(name) for name in names]
    sources += [Path('docs/repeated_development_gpu.md'), Path('tests/test_repeated_development_gpu.py')]
    protocol=dict(jobs=jobs,duration=18.,interval=.15,checkpoint_interval=3.,snapshots=[0.,3.,6.,9.,12.,15.,18.],
        cpu_jobs=[job['key'] for job in jobs if job['level']=='coarse'],
        backend_limits=BACKEND_LIMITS,time_limits=TIME_LIMITS,space_limits=SPACE_LIMITS,physical_limits=PHYSICAL,
        homologous_limits=HOMOLOGOUS_LIMITS,
        restart_steps=8,restart_stages=[2,4,8],minimum_mature_window=3.,
        probes=[112,128],orders=[3,5],primary_probe=112,primary_order=3,
        calibration_limit=.0025,diagnostic_spread_limit=.001,trend_floor=1e-10,
        source_sha256=pin(sources),input_sha256=pin(inputs),method=METHOD,
        gpu_uuid='GPU-2991ab3d-dcbd-41f3-eac1-133be278f56e',seed=7,new_scientific_histories=0,
        preserved_manifest=str(Path('archive/study-document-snapshots/2026-10-10-before-repeated-development-gpu/manifest.json').resolve()),
        scope='Naturally scheduled synthetic repeated cleavage; zero new scientific histories. Temporal/spatial observable and spatial-density screens, no identity or continuum inference.')
    write(root/'protocol.json',protocol); write(root/'status.json',dict(state='prepared',protocol_sha256=digest(root/'protocol.json')))
    return protocol


def calibrate(root):
    root=Path(root); p=read(root/'protocol.json'); verify(p)
    if (root/'calibration.json').exists(): raise FileExistsError('Calibration already recorded')
    geometries=[dict(name='initial_zygote',radius=.8,center=[0,0,0],stretch=[1,1,1],angle=0.),
        dict(name='translated_small',radius=.36,center=[.017,-.023,.031],stretch=[1,1,1],angle=0.),
        dict(name='rotated_ellipsoid',radius=.44,center=[-.029,.013,.019],stretch=[.9,1.05,1.1],angle=.43)]
    rows=[]
    for g in geometries:
        args={key:value for key,value in g.items() if key!='name'}
        for probe in p['probes']:
            expected=occupancy(analytic(probe,2.24,.085,**args)).astype(float)
            for grid in (72,88,104):
                native=analytic(grid,2.24,.085,**args)
                for order in p['orders']:
                    projected=reconstruct(native,2.24,probe,order)
                    err=float(np.linalg.norm((occupancy(projected).astype(float)-expected).ravel())/np.linalg.norm(expected.ravel()))
                    rows.append(dict(geometry=g['name'],grid=grid,probe=probe,order=order,relative_l2=err,
                        phase_range=[float(projected.min()),float(projected.max())]))
    result=dict(passed=all(row['relative_l2']<=p['calibration_limit'] for row in rows),rows=rows,
        protocol_sha256=digest(root/'protocol.json'),source_pins_verified=True)
    write(root/'calibration.json',result)
    write(root/'calibration-manifest.json',dict(protocol_sha256=digest(root/'protocol.json'),calibration_sha256=digest(root/'calibration.json')))
    return result


def event_passed(event, limits):
    return (event['neck']<=.15 and event['overlap']<=.002 and .1<event['fraction']<.9 and
        event['daughter_components']==[1,1] and all(event[key]<=limits[limit] for key,limit in (
            ('amount_relative_jump','event_amount_jump_max'),('target_relative_jump','event_target_jump_max'),
            ('volume_relative_jump','event_volume_jump_max'),('effective_occupancy_absolute_jump','event_occupancy_jump_max'))))


def path_quality(audit, limits):
    return (audit['max_volume_error']<limits['volume_error_max'] and audit['min_radius']>=limits['radius_min'] and
        audit['max_clipping']==0 and audit['dilution_amount_error']<=limits['dilution_amount_error_max'])


def run_path(root, job, arm):
    root=Path(root); p=read(root/'protocol.json'); verify(p); folder=root/job['key']
    output=folder/(arm+'-result.json')
    if output.exists(): raise FileExistsError('Never rerun a recorded path')
    gpu=arm=='gpu'; sim=DevelopmentGpu.restore(job['source']) if gpu else reference(job['source'])
    history=[]; audit=dict(max_volume_error=0.,min_radius=1e100,max_clipping=0.,dilution_amount_error=0.)
    started=time.perf_counter(); checkpoints=[]; observed_events=0; restart_checks={}; clone=None; remaining=0; stage=None
    independent_error=0.; initial_target=float(sim.target.sum()); failure=None
    try:
        for k in range(job['config']['steps']+1):
            now=sim.audit() if gpu else audit_cpu(sim)
            for key,value in now.items():
                if key in audit: audit[key]=min(audit[key],value) if key=='min_radius' else max(audit[key],value)
            if not path_quality(audit,p['physical_limits']): raise FloatingPointError('Every-step physical gate failed')
            u=torch.stack((sim.activator,sim.inhibitor)).cpu().numpy() if gpu else np.stack((sim.activator,sim.inhibitor))
            if not np.isfinite(u).all() or u.min()<=0: raise FloatingPointError('Chemical positivity/finite gate failed')
            if abs(float(sim.target.sum())/initial_target-1)>2e-14: raise FloatingPointError('Target inheritance drift')
            for ev in sim.event_audits[observed_events:]:
                if not event_passed(ev,p['physical_limits']): raise FloatingPointError('Event gate failed')
            observed_events=len(sim.event_audits)
            active=sim.host.divisions if gpu else sim.divisions
            if k%job['config']['save_every']==0:
                record=row(sim,gpu); record['parents']=(sim.host.parents if gpu else sim.parents).tolist()
                record['spectrum']=spectrum(record['volumes'],record['delta']).tolist(); history.append(record)
                write(folder/(arm+'-history.json'),history)
                if record['boundary_occupancy']>=.01 or any(c!=1 for c in record['nondividing_components']):
                    raise FloatingPointError('Boundary/connectivity observation gate failed: '+str(record['nondividing_components']))
                write(folder/(arm+'-progress.json'),dict(state='running',time=sim.time,cells=len(sim.ids),
                    active_divisions=len(active),events=observed_events,audit=audit,
                    seconds=time.perf_counter()-started,seconds_per_step=(time.perf_counter()-started)/max(k,1),
                    estimated_remaining_seconds=(time.perf_counter()-started)/max(k,1)*(job['config']['steps']-k)))
            if k%_steps(p['checkpoint_interval'],job['config']['dt'])==0:
                cp=folder/(arm+f'-t{int(round(sim.time)):02d}.npz')
                sim.checkpoint(cp) if gpu else save_cpu(sim,cp); checkpoints.append(cp)
            if gpu and clone is None:
                for candidate in p['restart_stages']:
                    if candidate not in restart_checks and candidate<len(sim.ids)<2*candidate and active:
                        stage=candidate; cp=folder/f'restart-stage-{stage}.npz'; sim.checkpoint(cp)
                        clone=DevelopmentGpu.restore(cp); exact(sim,clone); remaining=p['restart_steps']; break
            if k==job['config']['steps']: break
            if gpu and k<8:
                expected, correction=NumpyStepper(sim.matrices()[1].cpu().numpy(),job['config']['dt'],sim.config.signal_beta,
                    sim.config.signal_da,sim.config.signal_dh)(u,sim.chemical_correction.cpu().numpy())
                old=sim.geometry[:,0].cpu().numpy().copy()
            sim.step()
            if gpu and k<8:
                expected,correction=numpy_scale(expected,correction,old/sim.geometry[:,0].cpu().numpy())
                actual=torch.stack((sim.activator,sim.inhibitor)).cpu().numpy()
                independent_error=max(independent_error,float(abs(np.log(actual/expected)).max()))
                if independent_error>1e-11: raise FloatingPointError('Independent chemical/dilution prefix failed')
            if clone is not None:
                clone.step(); exact(sim,clone); remaining-=1
                if remaining==0: restart_checks[stage]=p['restart_steps']; clone=None; stage=None
        if len(sim.ids)!=16 or len(sim.event_audits)!=15 or active: raise RuntimeError('Repeated development unfinished')
        if p['duration']-max(ev['time'] for ev in sim.event_audits)<p['minimum_mature_window']:
            raise RuntimeError('Insufficient post-cleavage window')
        if gpu and set(restart_checks)!=set(p['restart_stages']): raise RuntimeError('Missing mixed-count restart window')
    except Exception as exc:
        failure=dict(type=type(exc).__name__,message=str(exc),time=sim.time,cells=len(sim.ids))
        def host_array(value):return value.cpu().numpy() if torch.is_tensor(value) else np.asarray(value)
        np.savez_compressed(folder/(arm+'-failure-arrays.npz'),
            **{key:host_array(getattr(sim,key)) for key in ('phi','phase_carry','activator','inhibitor','polarity','target','chemical_correction')},
            ids=sim.ids,diagnostic_only=np.array(True),failure=np.array(json.dumps(failure)))
    lineage=sim.host.lineage if gpu else sim.lineage
    write(folder/(arm+'-history.json'),history)
    write(folder/(arm+'-lineage.json'),lineage); write(folder/(arm+'-events.json'),sim.event_audits)
    evidence=[folder/(arm+'-history.json'),folder/(arm+'-lineage.json'),folder/(arm+'-events.json'),*checkpoints]
    if failure:evidence.append(folder/(arm+'-failure-arrays.npz'))
    if gpu:evidence += list(folder.glob('restart-stage-*.npz'))
    result=dict(passed=failure is None,job=job,arm=arm,failure=failure,audit=audit,events=sim.event_audits,
        exact_restart_steps={str(key):value for key,value in restart_checks.items()},independent_chemical_steps=8 if gpu else None,
        independent_chemical_log_max=independent_error if gpu else None,final_cells=len(sim.ids),final_time=sim.time,
        seconds=time.perf_counter()-started,protocol_sha256=digest(root/'protocol.json'),evidence_sha256=pin(evidence))
    write(output,result); write(folder/(arm+'-progress.json'),dict(state='completed' if result['passed'] else 'failed',
        time=sim.time,cells=len(sim.ids),events=len(sim.event_audits),seconds=result['seconds'],audit=audit,failure=failure))
    del sim,clone
    if gpu:torch.cuda.empty_cache()
    return result


def cpu_worker(args):
    return run_path(*args,'cpu')


def checkpoint_density(path,probe=None,order=3):
    with np.load(path,allow_pickle=False) as z:
        phi=z['phi']; u=np.stack((z['activator'],z['inhibitor'])); polarity=z['polarity']
        extent=json.loads(str(z['metadata']))['config']['extent']
        if probe is not None: phi=reconstruct(phi,extent,probe,order)
        result=densities(phi,u,polarity)
        phase_range=[float(phi.min()),float(phi.max())]
    return result,phase_range


def field_comparison(root,first,second,arms=('gpu','gpu'),spatial=False):
    rows=[]
    for stamp in (0,3,6,9,12,15,18):
        variants=[(112,3)] if spatial else [(None,3)]
        if spatial and stamp==18: variants=[(112,3),(112,5),(128,3),(128,5)]
        for probe,order in variants:
            a,ar=checkpoint_density(root/first/(arms[0]+f'-t{stamp:02d}.npz'),probe,order)
            b,br=checkpoint_density(root/second/(arms[1]+f'-t{stamp:02d}.npz'),probe,order)
            rows.append(dict(time=stamp,probe=probe,order=order,errors=density_errors(a,b),phase_ranges=[ar,br]))
    return rows


def comparison(root,first,second,limits,arms=('gpu','gpu'),spatial=False):
    reports=[read(root/name/(arm+'-result.json')) for name,arm in zip((first,second),arms)]
    if not all(r['passed'] for r in reports): return dict(passed=False,cases=[first,second],arms=arms,incomplete_path=True)
    history=[read(root/name/(arm+'-history.json')) for name,arm in zip((first,second),arms)]
    error,unmatched=compare_observations(*history)
    maturity=abs(max(e['time'] for e in reports[0]['events'])-max(e['time'] for e in reports[1]['events']))
    fields=field_comparison(root,first,second,arms,spatial)
    field_max=max(value for item in fields for value in item['errors'].values())
    passed=all(value<=limits[key] for key,value in error.items()) and maturity<=limits['maturity_time_difference_max'] and field_max<=limits['field_max']
    spread={key:max(r['errors'][key] for r in fields if r['time']==18)-min(r['errors'][key] for r in fields if r['time']==18) for key in FIELDS}
    if spatial:passed=passed and max(spread.values())<=.001
    extra={}
    if arms==('gpu','cpu'):
        trees=[[(r['id'],r['parent']) for r in read(root/first/(arm+'-lineage.json'))] for arm in arms]
        events=[{e['parent']:e for e in report['events']} for report in reports]
        correspondence=trees[0]==trees[1] and set(events[0])==set(events[1])
        event_error=0.
        if correspondence:
            for parent,e in events[0].items():
                other=events[1][parent];correspondence=correspondence and e['children']==other['children'] and e['draws']==other['draws']
                event_error=max(event_error,abs(e['time']-other['time']))
        extra=dict(lineage_and_partition_correspondence=bool(correspondence),parent_event_time_max=event_error)
        passed=passed and correspondence and event_error<=.0075+1e-12
        if correspondence:
            raw,raw_unmatched=homologous_errors(*history);phase_max=0.;unmatched_snapshots=[]
            for stamp in (0,3,6,9,12,15,18):
                with np.load(root/first/f'gpu-t{stamp:02d}.npz') as a,np.load(root/first/f'cpu-t{stamp:02d}.npz') as b:
                    if not np.array_equal(a['ids'],b['ids']) or not np.array_equal(a['parents'],b['parents']):
                        unmatched_snapshots.append(stamp);continue
                    phase_max=max(phase_max,float(abs(a['phi']-b['phi']).max()))
            raw['phase_abs_max']=phase_max
            passed=passed and all(value<=HOMOLOGOUS_LIMITS[key] for key,value in raw.items()) and 18 not in unmatched_snapshots
            extra.update(homologous_errors=raw,homologous_unmatched_times=raw_unmatched,homologous_unmatched_snapshots=unmatched_snapshots)
    return dict(passed=bool(passed),cases=[first,second],arms=arms,errors=error,
        unmatched_times=unmatched,maturity_time_error=maturity,fields=fields,field_max=field_max,
        diagnostic_spread=spread,**extra)


def assess(root):
    root=Path(root);p=read(root/'protocol.json');verify(p);rows=[];backends=[];temporal=[];spatial=[];trends=[];evidence=[]
    for job in p['jobs']:
        for arm in (('gpu','cpu') if job['key'] in p['cpu_jobs'] else ('gpu',)):
            report=read(root/job['key']/(arm+'-result.json'))
            if report['protocol_sha256']!=digest(root/'protocol.json'): raise ValueError('Path protocol changed')
            for file,value in report['evidence_sha256'].items():
                if digest(file)!=value: raise ValueError('Path evidence changed: '+file)
            rows.append(dict(case=job['key'],arm=arm,passed=report['passed'],failure=report['failure']))
            evidence.append(root/job['key']/(arm+'-result.json'))
    for branch in BRANCHES:
        backends.append(comparison(root,branch+'_coarse',branch+'_coarse',p['backend_limits'],('gpu','cpu')))
        temporal.append(comparison(root,branch+'_coarse',branch+'_fine',p['time_limits']))
        pairs=[comparison(root,branch+'_'+a,branch+'_'+b,p['space_limits'],spatial=True) for a,b in (('fine','space88'),('space88','space104'))]
        spatial.extend(pairs)
        trend=False
        if all('errors' in pair for pair in pairs):
            keys=('concentration_wasserstein_max','log_wasserstein_max','relative_axis_max')
            trend=all(pairs[1]['errors'][key]<=pairs[0]['errors'][key]+p['trend_floor'] for key in keys)
            for key in FIELDS:
                maxima=[max(r['errors'][key] for r in pair['fields']) for pair in pairs]
                trend=trend and maxima[1]<=maxima[0]+p['trend_floor']
        trends.append(dict(branch=branch,passed=bool(trend)))
    verify(p)
    result=dict(passed=all(row['passed'] for row in rows+backends+temporal+spatial+trends),
        independent_review_pending=True,protocol_sha256=digest(root/'protocol.json'),new_scientific_histories=0,
        paths=rows,backend_rows=backends,time_rows=temporal,space_rows=spatial,spatial_trend_rows=trends,
        evidence_sha256=pin(evidence),scope=p['scope'])
    write(root/'summary.json',result);return result


def run(root,workers=4):
    root=Path(root).resolve()
    if workers not in (1,2,3,4):raise ValueError('Use one to four independent CPU workers')
    with (root/'coordinator.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        p=read(root/'protocol.json');verify(p);status=read(root/'status.json')
        if status['state']!='prepared' or status['protocol_sha256']!=digest(root/'protocol.json'): raise ValueError('Unchanged prepared study required')
        calibration=read(root/'calibration.json');manifest=read(root/'calibration-manifest.json')
        if (not calibration['passed'] or calibration['protocol_sha256']!=digest(root/'protocol.json') or
                manifest['protocol_sha256']!=digest(root/'protocol.json') or manifest['calibration_sha256']!=digest(root/'calibration.json') or
                len(calibration['rows'])!=36 or any(r['relative_l2']>p['calibration_limit'] for r in calibration['rows'])):
            raise ValueError('Reconstruction calibration failed or changed')
        if os.environ.get('CUDA_VISIBLE_DEVICES')!=p['gpu_uuid']:raise ValueError('Explicit qualified GPU UUID required')
        torch.set_num_threads(1);started=time.perf_counter()
        try:
            lib=library();write(root/'runtime.json',dict(gpu=torch.cuda.get_device_name(0),torch=torch.__version__,
                binary=lib.binary_path,binary_sha256=digest(lib.binary_path)))
            with ProcessPoolExecutor(max_workers=workers,mp_context=get_context('spawn')) as pool:
                futures=[pool.submit(cpu_worker,(str(root),job)) for job in p['jobs'] if job['key'] in p['cpu_jobs']]
                completed=[]
                for job in p['jobs']:
                    write(root/'status.json',dict(state='running',case=job['key'],completed_gpu=completed,total_gpu=16,total_cpu=4))
                    result=run_path(root,job,'gpu');completed.append(job['key'])
                    print(job['key'],result['passed'],result['seconds'],result['failure'],flush=True)
                    for future in futures:
                        if future.done():future.result()
                write(root/'status.json',dict(state='running',stage='waiting_cpu_references',completed_gpu=completed,total_gpu=16,total_cpu=4))
                for future in futures:future.result()
            write(root/'status.json',dict(state='running',stage='field_and_refinement_assessment',completed_gpu=completed,total_gpu=16,total_cpu=4))
            result=assess(root);write(root/'timing.json',dict(seconds=time.perf_counter()-started))
            write(root/'status.json',dict(state='completed' if result['passed'] else 'completed_with_unresolved_checks',
                passed=result['passed'],independent_review_pending=True,total_gpu=16,total_cpu=4))
            return result
        except Exception as exc:
            write(root/'status.json',dict(state='failed',error=str(exc)));raise


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('prepare','calibrate','run','assess'))
    parser.add_argument('--output',type=Path,default=Path('output/studies/repeated-development-gpu-qualification'))
    parser.add_argument('--workers',type=int,default=4);args=parser.parse_args()
    result=run(args.output,args.workers) if args.action=='run' else globals()[args.action](args.output)
    print(json.dumps(dict(passed=result.get('passed'),action=args.action),indent=2))


if __name__=='__main__':main()
