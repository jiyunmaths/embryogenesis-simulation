"""Prospective feedback-off spatial ladder with finest-grid timestep check.

Reuse the pinned diagnostic runner/comparisons; do not change solver dynamics.
The single reused path is read-only and retains its original provenance.
"""
import argparse
from dataclasses import asdict, replace
from pathlib import Path
import fcntl
import json
import os
import time

import numpy as np
import torch

from . import development_time_diagnosis as td
from . import qualify_repeated_development as q
from .model import Config, occupancy
from .resolution import _steps
from .qualify_development_gpu import digest, read, write, pin, save_cpu
from .cleavage_measurement import analytic, reconstruct

PARENT = Path('output/studies/development-time-diagnosis')
PRESERVED = Path('archive/study-document-snapshots/2026-10-10-before-development-space-refinement/manifest.json')
LEVELS = (('off72',72,.00046875,True),('off88',88,.00046875,False),
    ('off104',104,.00046875,False),('off120',120,.00046875,False),
    ('off120_half',120,.000234375,False))


def config_changes(base, candidate):
    changed = {k for k,v in candidate.items() if v != base[k]}
    if not changed <= {'grid','dt','steps','save_every'}:
        raise ValueError('Physical model changed: '+str(changed))
    if candidate['steps'] != _steps(18.,candidate['dt']) or candidate['save_every'] != _steps(.15,candidate['dt']):
        raise ValueError('Physical horizon or observation clock changed')
    if candidate['feedback'] or any(candidate[k] for k in ('fate_tension','fate_adhesion','polarity_tension')):
        raise ValueError('Only feedback-off geometry is authorized')
    return sorted(changed)


def prepare(root):
    root = Path(root).resolve()
    if root.exists(): raise FileExistsError('Fresh study directory required')
    assessment_path = Path('docs/development_time_diagnosis_assessment.json')
    accepted = read(assessment_path); review_path = Path(accepted['independent_review'])
    if not accepted['review_completed'] or not accepted['temporal_screen_passed']:
        raise ValueError('Reviewed feedback-off time screen required')
    if digest(review_path) != accepted['independent_review_sha256']: raise ValueError('Time review changed')
    parent = read(PARENT/'protocol.json'); review = read(review_path)
    assert review['protocol_sha256'] == digest(PARENT/'protocol.json') == accepted['protocol_sha256']
    preserved = read(PRESERVED)
    for file,expected in preserved['preserved_sha256'].items():
        if digest(file) != expected: raise ValueError('Preserved evidence changed: '+file)
    baseline = next(j for j in parent['jobs'] if j['key'] == 'off_quarter')
    baseline_folder = (PARENT/baseline['key']).resolve()
    baseline_report = read(baseline_folder/'gpu-result.json')
    assert baseline_report['passed'] and baseline_report['job'] == baseline
    assert baseline_report['protocol_sha256'] == digest(PARENT/'protocol.json')
    root.mkdir(parents=True); jobs = []
    inputs = [assessment_path, review_path, PARENT/'protocol.json', PARENT/'summary.json', PRESERVED,
        PARENT/'instrumentation-gate.json', PARENT/'runtime.json',
        Path('output/studies/repeated-development-gpu-qualification-v2/calibration.json')]
    inputs += [Path(f) for f in baseline_report['evidence_sha256']]
    inputs += [baseline_folder/'gpu-result.json', Path(baseline['source'])]
    for key,grid,dt,reused in LEVELS:
        cfg = replace(Config(**baseline['config']),grid=grid,dt=dt,steps=_steps(18.,dt),save_every=_steps(.15,dt))
        config_changes(baseline['config'],asdict(cfg))
        folder = root/key
        if reused:
            folder.symlink_to(baseline_folder,target_is_directory=True)
            source = Path(baseline['source'])
        else:
            folder.mkdir(); source = folder/'source.npz'; save_cpu(q.initial(cfg),source); inputs.append(source)
        jobs.append(dict(key=key,level='reused' if reused else key,branch='off',
            source=str(source.resolve()),config=asdict(cfg),reused=reused))
    sources = dict(parent['source_sha256'])
    sources.update(pin([Path(__file__),Path('docs/development_space_refinement.md'),Path('tests/test_development_space_refinement.py')]))
    p = dict(jobs=jobs,duration=18.,interval=.15,checkpoint_interval=3.,diagnostic_interval=.015,
        snapshots=[0.,3.,6.,9.,12.,15.,18.],time_limits=dict(parent['time_limits']),
        space_limits=dict(q.SPACE_LIMITS),physical_limits=dict(parent['physical_limits']),
        restart_steps=8,restart_stages=[2,4,8],minimum_mature_window=3.,
        space_pairs=[['off72','off88'],['off88','off104'],['off104','off120']],
        time_pair=['off120','off120_half'],primary_probe=112,primary_order=3,
        probes=[112,128],orders=[3,5],calibration_limit=.0025,diagnostic_spread_limit=.001,trend_floor=1e-10,
        dense_time_spectrum_limit=.01,source_sha256=sources,input_sha256=pin(inputs),method=parent['method'],
        gpu_uuid=parent['gpu_uuid'],new_scientific_histories=0,new_numerical_paths=4,reused_numerical_paths=1,
        preserved_manifest=str(PRESERVED.resolve()),
        reuse_provenance=dict(folder=str(baseline_folder),job=baseline,report_sha256=digest(baseline_folder/'gpu-result.json'),
            protocol_sha256=digest(PARENT/'protocol.json')),
        instrumentation_gate=str((PARENT/'instrumentation-gate.json').resolve()),
        scope='Feedback-off spatial/time numerical screen, same synthetic seed7, t0–18. Finest spatial pair104→120 and grid120 time pair; zero new independent scientific histories. No mechanical-control/full-restart/fine-CPU or biological qualification.')
    write(root/'protocol.json',p); write(root/'status.json',dict(state='prepared',protocol_sha256=digest(root/'protocol.json')))
    return p


def calibrate(root):
    root = Path(root); p = read(root/'protocol.json'); q.verify(p)
    if (root/'calibration.json').exists(): raise FileExistsError('Calibration already recorded')
    geometries = [dict(name='initial_zygote',radius=.8,center=[0,0,0],stretch=[1,1,1],angle=0.),
        dict(name='translated_small',radius=.36,center=[.017,-.023,.031],stretch=[1,1,1],angle=0.),
        dict(name='rotated_ellipsoid',radius=.44,center=[-.029,.013,.019],stretch=[.9,1.05,1.1],angle=.43)]
    rows = []
    for g in geometries:
        args = {k:v for k,v in g.items() if k != 'name'}
        for probe in p['probes']:
            expected = occupancy(analytic(probe,2.24,.085,**args)).astype(float)
            for grid in sorted({j['config']['grid'] for j in p['jobs']}):
                native = analytic(grid,2.24,.085,**args)
                for order in p['orders']:
                    projected = reconstruct(native,2.24,probe,order)
                    err = float(np.linalg.norm((occupancy(projected).astype(float)-expected).ravel())/np.linalg.norm(expected.ravel()))
                    rows.append(dict(geometry=g['name'],grid=grid,probe=probe,order=order,
                        relative_l2=err,phase_range=[float(projected.min()),float(projected.max())]))
    result = dict(passed=len(rows)==48 and all(i['relative_l2']<=p['calibration_limit'] for i in rows),
        rows=rows,protocol_sha256=digest(root/'protocol.json'))
    write(root/'calibration.json',result)
    write(root/'calibration-manifest.json',dict(protocol_sha256=digest(root/'protocol.json'),calibration_sha256=digest(root/'calibration.json')))
    return result


def case_report(root,job,p):
    folder = Path(root)/job['key']; report = read(folder/'gpu-result.json')
    if job['reused']:
        provenance = p['reuse_provenance']
        if folder.resolve() != Path(provenance['folder']) or digest(folder/'gpu-result.json') != provenance['report_sha256']:
            raise ValueError('Reused path or report changed')
        if report['job'] != provenance['job'] or report['protocol_sha256'] != provenance['protocol_sha256']:
            raise ValueError('Reused report lost original provenance')
        if report['job']['config'] != job['config']: raise ValueError('Reuse configuration mismatch')
    elif report['job'] != job or report['protocol_sha256'] != digest(Path(root)/'protocol.json'):
        raise ValueError('New case protocol or job changed')
    for file,expected in report['evidence_sha256'].items():
        if digest(file) != expected: raise ValueError('Recorded case evidence changed: '+file)
    return report


def execute_path(root,job):
    if job['reused'] or (Path(root)/job['key']).is_symlink():
        raise ValueError('Reused paths are read-only; never execute or overwrite')
    return td.run_path(root,job)


def dense_pair(root,first,second):
    def records(name):
        rows = [json.loads(l) for l in (Path(root)/name/'network.jsonl').open()]
        regular = {round(r['time'],9):r for r in rows if abs(r['time']/.015-round(r['time']/.015)) < 1e-7}
        expected = [round(i*.015,9) for i in range(1201)]
        if sorted(regular) != expected: raise ValueError('Incomplete dense physical clocks')
        return regular
    a,b = records(first),records(second); peak = dict(relative_spectrum=0.,time=0.)
    unmatched = []; topology = []
    for t,x in a.items():
        y = b[t]
        sx,sy = q.spectrum(x['volumes'],x['delta']),q.spectrum(y['volumes'],y['delta'])
        if len(sx) != len(sy): unmatched.append(t); continue
        er = float(np.linalg.norm(sx-sy)/(np.linalg.norm(sy) or 1.))
        if er > peak['relative_spectrum']:
            peak = dict(relative_spectrum=er,time=t,cells=len(sx),edges=[x['edges'],y['edges']])
        if x['ids'] == y['ids'] and x['parents'] == y['parents']:
            if not np.array_equal(np.array(x['conductance'])>0,np.array(y['conductance'])>0): topology.append(t)
    return dict(cases=[first,second],peak=peak,unmatched_times=unmatched,topology_mismatch_times=topology)


def spatial_trend(pairs,floor):
    if not all('errors' in r and r.get('fields') for r in pairs): return False
    sequences = {k:[r['errors'][k] for r in pairs] for k in
        ('concentration_wasserstein_max','log_wasserstein_max','relative_axis_max')}
    sequences.update({k:[max(f['errors'][k] for f in r['fields']) for r in pairs] for k in q.FIELDS})
    return all(all(b <= a+floor for a,b in zip(values,values[1:])) for values in sequences.values())


def scoped_decision(paths,space,time_row,dense_time,trend,limit):
    return bool(all(r['passed'] for r in paths) and space[-1]['passed'] and time_row['passed']
        and dense_time['peak']['relative_spectrum'] <= limit and trend)


def assess(root):
    root = Path(root); p = read(root/'protocol.json'); q.verify(p)
    if (root/'summary.json').exists(): raise FileExistsError('Never overwrite a recorded assessment')
    reports = [case_report(root,job,p) for job in p['jobs']]
    spaces = [q.comparison(root,a,b,p['space_limits'],spatial=True) for a,b in p['space_pairs']]
    times = q.comparison(root,*p['time_pair'],p['time_limits'])
    dense_time = dense_pair(root,*p['time_pair'])
    dense_spaces = [dense_pair(root,a,b) for a,b in p['space_pairs']]
    trend = spatial_trend(spaces,p['trend_floor'])
    decision = scoped_decision(reports,spaces,times,dense_time,trend,p['dense_time_spectrum_limit'])
    q.verify(p)
    result = dict(execution_completed=True,all_paths_passed=all(r['passed'] for r in reports),
        finest_space_screen_passed=spaces[-1]['passed'],all_space_pairs_passed=all(r['passed'] for r in spaces),
        fine_grid_time_screen_passed=times['passed'],dense_time_screen_passed=dense_time['peak']['relative_spectrum']<=p['dense_time_spectrum_limit'],
        spatial_error_decay_passed=trend,scoped_screens_passed=decision,qualification_accepted=False,
        independent_review_pending=True,new_scientific_histories=0,protocol_sha256=digest(root/'protocol.json'),
        space_rows=spaces,time_row=times,dense_time=dense_time,dense_space_rows=dense_spaces,
        paths=[dict(case=j['key'],reused=j['reused'],passed=r['passed'],seconds=r['seconds']) for j,r in zip(p['jobs'],reports)],
        evidence_sha256=pin(root/j['key']/'gpu-result.json' for j in p['jobs']),scope=p['scope'])
    write(root/'summary.json',result); return result


def run(root):
    root = Path(root).resolve()
    with (root/'coordinator.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        p = read(root/'protocol.json'); q.verify(p)
        status = read(root/'status.json')
        if status['state'] != 'prepared' or status['protocol_sha256'] != digest(root/'protocol.json'):
            raise ValueError('Unchanged prepared study required')
        calibration = read(root/'calibration.json'); manifest = read(root/'calibration-manifest.json')
        if not calibration['passed'] or len(calibration['rows']) != 48 or any(i['relative_l2']>p['calibration_limit'] for i in calibration['rows']):
            raise ValueError('Spatial reconstruction calibration failed')
        if calibration['protocol_sha256'] != digest(root/'protocol.json') or manifest['protocol_sha256'] != digest(root/'protocol.json') or manifest['calibration_sha256'] != digest(root/'calibration.json'):
            raise ValueError('Calibration changed')
        gate = read(p['instrumentation_gate'])
        if not gate['passed'] or gate['protocol_sha256'] != p['reuse_provenance']['protocol_sha256']:
            raise ValueError('Pinned instrumentation evidence required')
        for file,expected in gate['evidence_sha256'].items():
            if digest(file) != expected: raise ValueError('Instrumentation evidence changed')
        if os.environ.get('CUDA_VISIBLE_DEVICES') != p['gpu_uuid']: raise ValueError('Explicit designated GPU required')
        torch.set_num_threads(1); started = time.perf_counter(); completed = []
        try:
            lib = td.gpu.library()
            write(root/'runtime.json',dict(gpu=torch.cuda.get_device_name(0),torch=torch.__version__,binary=lib.binary_path,binary_sha256=digest(lib.binary_path)))
            for job in p['jobs']:
                if job['reused']:
                    result = case_report(root,job,p)
                    if not result['passed']: raise ValueError('Reused baseline failed')
                else:
                    write(root/'status.json',dict(state='running',case=job['key'],completed=completed,total=5,new_total=4))
                    result = execute_path(root,job)
                    if not result['passed']: raise RuntimeError('Path physical/event gate failed; preserve result before further work')
                completed.append(job['key'])
                print(job['key'],'reused' if job['reused'] else 'new',result['passed'],result['seconds'],flush=True)
            write(root/'status.json',dict(state='running',stage='fixed_time_refinement_assessment',completed=completed,total=5))
            result = assess(root); write(root/'timing.json',dict(seconds=time.perf_counter()-started))
            write(root/'status.json',dict(state='completed_provisional' if result['scoped_screens_passed'] else 'completed_with_unresolved_checks',
                scoped_screens_passed=result['scoped_screens_passed'],qualification_accepted=False,independent_review_pending=True,
                total=5,new_total=4,reused_total=1)); return result
        except Exception as exc:
            write(root/'status.json',dict(state='failed',error=str(exc),completed=completed)); raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('prepare','calibrate','run','assess'))
    parser.add_argument('--output',type=Path,default=Path('output/studies/development-space-refinement'))
    args = parser.parse_args(); result = globals()[args.action](args.output)
    print(json.dumps({k:result[k] for k in ('passed','scoped_screens_passed','qualification_accepted') if k in result}),flush=True)


if __name__ == '__main__': main()
