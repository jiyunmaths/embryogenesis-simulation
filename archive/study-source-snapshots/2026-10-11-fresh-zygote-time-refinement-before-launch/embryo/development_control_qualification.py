"""Matched mechanical-control refinement and full live-cleavage restart replays.

Uses the pinned developmental solver and comparisons without changing dynamics.
Production and replay evidence are written only to a fresh study directory.
"""
import argparse
from dataclasses import asdict, replace
from pathlib import Path
import fcntl
import os
import time

import numpy as np
import torch

from . import development_time_diagnosis as td
from . import development_space_refinement as space
from . import qualify_repeated_development as q
from .model import Config
from .qualify_development_gpu import digest, read, write, pin, save_cpu, row, exact
from .resolution import _steps

PARENT = Path('output/studies/development-space-refinement')
ASSESSMENT = Path('docs/development_space_refinement_assessment.json')
PRESERVED = Path('archive/study-document-snapshots/2026-10-10-before-development-controls/manifest.json')
BRANCHES = ('tension_adhesion', 'polarity', 'full')
LEVELS = (('space88', 88, .00046875), ('space104', 104, .00046875),
          ('space120', 120, .00046875), ('time120_half', 120, .000234375))
COUPLINGS = {'off': (False, 0., 0., 0.),
    'tension_adhesion': (True, .25, .35, 0.),
    'polarity': (True, 0., 0., .35), 'full': (True, .25, .35, .35)}


def config_changes(base, candidate, branch):
    if branch not in COUPLINGS:
        raise ValueError('Unknown matched mechanical branch')
    changed = {key for key, value in candidate.items() if value != base[key]}
    allowed = {'grid', 'dt', 'steps', 'save_every', 'feedback',
               'fate_tension', 'fate_adhesion', 'polarity_tension'}
    if not changed <= allowed:
        raise ValueError('Physical model drift: '+str(changed-allowed))
    expected = COUPLINGS[branch]
    actual = tuple(candidate[key] for key in
        ('feedback', 'fate_tension', 'fate_adhesion', 'polarity_tension'))
    if actual != expected:
        raise ValueError('Mechanical branch is not the prescribed matched control')
    if (candidate['steps'] != _steps(18., candidate['dt']) or
            candidate['save_every'] != _steps(.15, candidate['dt'])):
        raise ValueError('Physical horizon or observation clock drift')
    return sorted(changed)


def verify_preserved(p):
    q.verify(p)
    for name, expected in read(p['preserved_manifest'])['preserved_sha256'].items():
        if digest(name) != expected:
            raise ValueError('Preserved evidence changed: '+name)


def prepare(root):
    root = Path(root).resolve()
    if root.exists():
        raise FileExistsError('Fresh qualification directory required')
    accepted = read(ASSESSMENT)
    review_path = Path(accepted['independent_review']['path'])
    if not accepted['review_completed'] or not accepted['scoped_screen_accepted']:
        raise ValueError('Independently accepted feedback-off screen required')
    if digest(review_path) != accepted['independent_review']['sha256']:
        raise ValueError('Accepted feedback-off review changed')
    parent = read(PARENT/'protocol.json')
    if accepted['protocol_sha256'] != digest(PARENT/'protocol.json'):
        raise ValueError('Accepted parent protocol changed')
    for name, expected in read(PRESERVED)['preserved_sha256'].items():
        if digest(name) != expected:
            raise ValueError('Preserved evidence changed: '+name)
    baseline = next(j for j in parent['jobs'] if j['key'] == 'off120')
    original_folder = (PARENT/'off120').resolve()
    original_report = read(original_folder/'gpu-result.json')
    if (not original_report['passed'] or original_report['job'] != baseline or
            original_report['protocol_sha256'] != digest(PARENT/'protocol.json')):
        raise ValueError('Accepted grid-120 baseline report required')
    root.mkdir(parents=True)
    jobs = []; inputs = [ASSESSMENT, review_path, PARENT/'protocol.json', PARENT/'summary.json',
        PRESERVED, PARENT/'runtime.json', PARENT/'calibration.json', PARENT/'calibration-manifest.json',
        Path(parent['instrumentation_gate'])]
    for branch in BRANCHES:
        feedback, tension, adhesion, polarity = COUPLINGS[branch]
        for level, grid, dt in LEVELS:
            cfg = replace(Config(**baseline['config']), grid=grid, dt=dt,
                steps=_steps(18., dt), save_every=_steps(.15, dt), feedback=feedback,
                fate_tension=tension, fate_adhesion=adhesion, polarity_tension=polarity)
            host = q.initial(cfg)
            config_changes(baseline['config'], asdict(host.config), branch)
            key = branch+'_'+level; folder = root/key; folder.mkdir()
            source = folder/'source.npz'; save_cpu(host, source); inputs.append(source)
            jobs.append(dict(key=key, branch=branch, level=level, source=str(source),
                config=asdict(host.config), reused=False))
    (root/'off120').symlink_to(original_folder, target_is_directory=True)
    jobs.append(dict(baseline, reused=True))
    inputs += [original_folder/'gpu-result.json', Path(baseline['source'])]
    inputs += [Path(name) for name in original_report['evidence_sha256']]
    sources = dict(parent['source_sha256'])
    sources.update(pin([Path(__file__), Path('docs/development_control_qualification.md'),
        Path('tests/test_development_control_qualification.py')]))
    p = {key:parent[key] for key in ('duration', 'interval', 'checkpoint_interval',
        'diagnostic_interval', 'snapshots', 'time_limits', 'space_limits', 'physical_limits',
        'restart_steps', 'restart_stages', 'minimum_mature_window', 'primary_probe', 'primary_order',
        'probes', 'orders', 'calibration_limit', 'diagnostic_spread_limit', 'trend_floor',
        'dense_time_spectrum_limit', 'method', 'gpu_uuid')}
    p.update(jobs=jobs, branches=list(BRANCHES), new_scientific_histories=0,
        new_numerical_paths=12, reused_numerical_paths=1,
        source_sha256=sources, input_sha256=pin(inputs), preserved_manifest=str(PRESERVED.resolve()),
        parent_protocol_sha256=digest(PARENT/'protocol.json'),
        instrumentation_gate=parent['instrumentation_gate'],
        instrumentation_protocol_sha256=parent['reuse_provenance']['protocol_sha256'],
        calibration_reference=str((PARENT/'calibration.json').resolve()),
        calibration_manifest=str((PARENT/'calibration-manifest.json').resolve()),
        runtime_binary_sha256=read(PARENT/'runtime.json')['binary_sha256'],
        space_pairs={branch:[[branch+'_space88', branch+'_space104'],
            [branch+'_space104', branch+'_space120']] for branch in BRANCHES},
        time_pairs={branch:[branch+'_space120', branch+'_time120_half'] for branch in BRANCHES},
        replay_cases=['off120']+[branch+'_space120' for branch in BRANCHES],
        full_restart_replays=12,
        reuse_provenance={'off120':dict(folder=str(original_folder), job=baseline,
            report_sha256=digest(original_folder/'gpu-result.json'),
            protocol_sha256=original_report['protocol_sha256'])},
        scope='Matched mechanical-control numerical qualification of one synthetic seed7 zygote, '
            't0–18. Three branches, spatial/time screens and full remaining-horizon restart replays '
            'at stages2/4/8 for all four candidate branches at grid120,dt0.00046875. '
            'Zero new independent histories; no biological, continuum or fine CPU qualification.')
    write(root/'protocol.json', p)
    write(root/'status.json', dict(state='prepared', protocol_sha256=digest(root/'protocol.json')))
    return p


def case_report(root, job, p):
    folder = Path(root)/job['key']; report = read(folder/'gpu-result.json')
    if job['reused']:
        old = p['reuse_provenance'][job['key']]
        if (not folder.is_symlink() or folder.resolve() != Path(old['folder']) or
                digest(folder/'gpu-result.json') != old['report_sha256']):
            raise ValueError('Reused baseline path or report changed')
        if report['job'] != old['job'] or report['protocol_sha256'] != old['protocol_sha256']:
            raise ValueError('Reused baseline lost original provenance')
        if job['config'] != report['job']['config']:
            raise ValueError('Reused baseline configuration drift')
    elif report['job'] != job or report['protocol_sha256'] != digest(Path(root)/'protocol.json'):
        raise ValueError('New path protocol changed')
    for name, expected in report['evidence_sha256'].items():
        if digest(name) != expected:
            raise ValueError('Case evidence changed: '+name)
    return report


def execute_path(root, job):
    if job['reused'] or (Path(root)/job['key']).is_symlink():
        raise ValueError('Reused baseline is read-only; never execute it')
    return td.run_path(root, job)


def replay_live_state(root, job, stage, p):
    """Replay one saved mixed-count live-furrow state to every later observation.

    No input is overwritten. Exact states are compared at all later regular
    checkpoints and history records, with final arrays saved as separate evidence.
    """
    root = Path(root); source_folder = root/job['key']
    target = root/'replays'/job['key']/f'stage-{stage}'
    if target.exists():
        raise FileExistsError('Never overwrite a recorded full restart replay')
    target.mkdir(parents=True)
    source = source_folder/f'restart-stage-{stage}.npz'
    evidence = pin([source, source_folder/'gpu-result.json', source_folder/'gpu-history.json',
        *[source_folder/f'gpu-t{stamp:02d}.npz' for stamp in (0,3,6,9,12,15,18)]])
    sim = td.gpu.DevelopmentGpu.restore(source)
    if asdict(sim.config) != job['config']:
        raise ValueError('Replay source configuration differs from its candidate')
    if not stage < len(sim.ids) < 2*stage or not sim.host.divisions:
        raise ValueError('Replay must start from the prescribed live mixed-count window')
    start = sim.time; first_step = sim.step_number; started = time.perf_counter()
    expected_rows = {int(round(x['time']/job['config']['dt'])):x
        for x in read(source_folder/'gpu-history.json')}
    checked = []; checked_history = []; observed = len(sim.event_audits)
    initial_target = float(sim.target.sum()); failure = None
    audit = dict(max_volume_error=0., min_radius=1e100, max_clipping=0., dilution_amount_error=0.)
    try:
        while True:
            now = sim.audit()
            for key in audit:
                audit[key] = min(audit[key], now[key]) if key=='min_radius' else max(audit[key], now[key])
            if not q.path_quality(audit, p['physical_limits']):
                raise FloatingPointError('Replay physical gate failed')
            chemistry = torch.stack((sim.activator, sim.inhibitor))
            if not bool(torch.isfinite(chemistry).all()) or float(chemistry.min()) <= 0:
                raise FloatingPointError('Replay chemical positivity failed')
            if abs(float(sim.target.sum())/initial_target-1)>2e-14:
                raise FloatingPointError('Replay target inheritance drift')
            for ev in sim.event_audits[observed:]:
                if not q.event_passed(ev, p['physical_limits']):
                    raise FloatingPointError('Replay event gate failed')
            observed = len(sim.event_audits)
            k = sim.step_number
            if k in expected_rows:
                record = row(sim, True); record['parents'] = sim.host.parents.tolist()
                record['spectrum'] = q.spectrum(record['volumes'], record['delta']).tolist()
                if record != expected_rows[k]:
                    raise RuntimeError('Full restart observation changed at '+str(sim.time))
                checked_history.append(sim.time)
                write(target/'progress.json', dict(state='running', time=sim.time, start=start,
                    checked_observations=len(checked_history), seconds=time.perf_counter()-started))
            if k % _steps(p['checkpoint_interval'], job['config']['dt']) == 0:
                expected = td.gpu.DevelopmentGpu.restore(source_folder/f'gpu-t{int(round(sim.time)):02d}.npz')
                exact(sim, expected); checked.append(sim.time); del expected
            if k == job['config']['steps']:
                break
            if k > job['config']['steps']:
                raise RuntimeError('Replay clock passed the fixed horizon')
            sim.step()
        if len(sim.ids)!=16 or len(sim.event_audits)!=15 or sim.host.divisions:
            raise RuntimeError('Replay did not finish development')
        sim.checkpoint(target/'final.npz')
    except Exception as exc:
        failure = dict(type=type(exc).__name__, message=str(exc), time=sim.time)
        sim.checkpoint(target/'failure-state.npz')
    for name, expected in evidence.items():
        if digest(name) != expected:
            raise ValueError('Full replay changed an input: '+name)
    saved = target/('final.npz' if failure is None else 'failure-state.npz')
    result = dict(passed=failure is None, failure=failure, case=job['key'], stage=stage,
        start=start, end=sim.time, steps=sim.step_number-first_step,
        checkpoint_times=checked, exact_observation_times=checked_history,
        complete_state_bitwise_equal=failure is None,
        final_cells=len(sim.ids), final_events=len(sim.event_audits), audit=audit,
        seconds=time.perf_counter()-started, protocol_sha256=digest(root/'protocol.json'),
        input_sha256=evidence, output_sha256=pin([saved]))
    write(target/'result.json', result)
    write(target/'progress.json', dict(state='completed' if failure is None else 'failed',
        time=sim.time, seconds=result['seconds'], failure=failure))
    del sim; torch.cuda.empty_cache()
    return result


def branch_decision(paths, spatial, temporal, dense, trend, limit):
    return space.scoped_decision(paths, spatial, temporal, dense, trend, limit)


def assess_branch(root, branch, p):
    jobs = [j for j in p['jobs'] if j['branch']==branch]
    reports = [case_report(root, job, p) for job in jobs]
    spatial = [q.comparison(root, *pair, p['space_limits'], spatial=True)
        for pair in p['space_pairs'][branch]]
    temporal = q.comparison(root, *p['time_pairs'][branch], p['time_limits'])
    dense = space.dense_pair(root, *p['time_pairs'][branch])
    dense_spatial = [space.dense_pair(root, *pair) for pair in p['space_pairs'][branch]]
    trend = space.spatial_trend(spatial, p['trend_floor'])
    return dict(branch=branch, space_rows=spatial, time_row=temporal, dense_time=dense,
        dense_space_rows=dense_spatial, spatial_error_decay_passed=trend,
        scoped_screens_passed=branch_decision(reports, spatial, temporal, dense, trend,
            p['dense_time_spectrum_limit']))


def final_decision(branches, replays, replay_count):
    required = {(case, stage) for case in
        ['off120']+[branch+'_space120' for branch in BRANCHES] for stage in (2,4,8)}
    return bool(len(branches)==3 and {r['branch'] for r in branches}==set(BRANCHES)
        and all(r['scoped_screens_passed'] for r in branches)
        and len(replays)==replay_count==len(required)
        and {(r['case'], r['stage']) for r in replays}==required
        and all(r['passed'] for r in replays))


def validate_prelaunch(root, p):
    verify_preserved(p)
    calibration = read(p['calibration_reference']); manifest = read(p['calibration_manifest'])
    if (not calibration['passed'] or len(calibration['rows'])!=48 or
            calibration['protocol_sha256'] != p['parent_protocol_sha256'] or
            manifest['protocol_sha256'] != p['parent_protocol_sha256'] or
            manifest['calibration_sha256'] != digest(p['calibration_reference']) or
            any(x['relative_l2']>p['calibration_limit'] for x in calibration['rows'])):
        raise ValueError('Pinned reconstruction calibration changed or failed')
    combinations = {(x['grid'], x['probe'], x['order']) for x in calibration['rows']}
    if not {(g, pr, order) for g in (88,104,120) for pr in p['probes'] for order in p['orders']} <= combinations:
        raise ValueError('Calibration does not cover selected grids')
    gate = read(p['instrumentation_gate'])
    if not gate['passed'] or gate['protocol_sha256'] != p['instrumentation_protocol_sha256']:
        raise ValueError('Unchanged instrumentation gate required')
    for name, expected in gate['evidence_sha256'].items():
        if digest(name) != expected:
            raise ValueError('Instrumentation evidence changed')
    if os.environ.get('CUDA_VISIBLE_DEVICES') != p['gpu_uuid']:
        raise ValueError('Explicit designated GPU required')
    status = read(Path(root)/'status.json')
    if status != dict(state='prepared', protocol_sha256=digest(Path(root)/'protocol.json')):
        raise ValueError('Unchanged prepared status required')


def run(root):
    root = Path(root).resolve()
    with (root/'coordinator.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
        p = read(root/'protocol.json'); validate_prelaunch(root, p)
        torch.set_num_threads(1); started = time.perf_counter()
        completed = []; branches = []; replays = []
        lib = td.gpu.library()
        if digest(lib.binary_path) != p['runtime_binary_sha256']:
            raise ValueError('Runtime binary changed')
        write(root/'runtime.json', dict(gpu=torch.cuda.get_device_name(0), torch=torch.__version__,
            binary=lib.binary_path, binary_sha256=digest(lib.binary_path)))
        try:
            for job in p['jobs']:
                if job['reused']:
                    report = case_report(root, job, p)
                else:
                    write(root/'status.json', dict(state='running', stage='production', case=job['key'],
                        completed=completed, new_total=12, replay_total=12))
                    report = execute_path(root, job)
                if not report['passed']:
                    raise RuntimeError('Physical/event path failed; preserve it before further work')
                completed.append(job['key'])
                print(job['key'], 'reused' if job['reused'] else 'new', report['seconds'], flush=True)
            for branch in BRANCHES:
                write(root/'status.json', dict(state='running', stage='refinement_assessment',
                    branch=branch, completed=completed, new_total=12, replay_total=12))
                result = assess_branch(root, branch, p); branches.append(result)
                write(root/(branch+'-refinement.json'), result)
                print(branch, 'refinement', result['scoped_screens_passed'], flush=True)
            # Replays are numerical restart checks, even when a refinement gate fails.
            # They cannot promote a failed refinement or authorize biological science.
            for key in p['replay_cases']:
                job = next(j for j in p['jobs'] if j['key']==key)
                for stage in p['restart_stages']:
                    write(root/'status.json', dict(state='running', stage='full_restart_replay',
                        case=key, live_stage=stage, completed=completed,
                        replays_completed=len(replays), replay_total=12))
                    result = replay_live_state(root, job, stage, p); replays.append(result)
                    print(key, 'full restart', stage, result['passed'], result['seconds'], flush=True)
                    if not result['passed']:
                        raise RuntimeError('Full live-event restart replay failed')
            verify_preserved(p)
            passed = final_decision(branches, replays, p['full_restart_replays'])
            result = dict(execution_completed=True, scoped_screens_passed=passed,
                qualification_accepted=False, independent_review_pending=True,
                branch_rows=branches, replay_rows=replays, new_scientific_histories=0,
                paths=[dict(case=j['key'], reused=j['reused'], passed=True) for j in p['jobs']],
                protocol_sha256=digest(root/'protocol.json'), scope=p['scope'],
                evidence_sha256=pin([*[root/j['key']/'gpu-result.json' for j in p['jobs']],
                    *root.glob('*-refinement.json'), *root.glob('replays/*/stage-*/result.json')]))
            write(root/'summary.json', result)
            write(root/'timing.json', dict(seconds=time.perf_counter()-started))
            write(root/'status.json', dict(state='completed_provisional' if passed else 'completed_with_unresolved_checks',
                scoped_screens_passed=passed, qualification_accepted=False, independent_review_pending=True,
                new_total=12, reused_total=1, replays_completed=len(replays)))
            return result
        except Exception as exc:
            write(root/'status.json', dict(state='failed', error=str(exc), completed=completed,
                replays_completed=len(replays)))
            raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'run'))
    parser.add_argument('--output', type=Path, default=Path('output/studies/development-control-qualification'))
    args = parser.parse_args(); result = globals()[args.action](args.output)
    print({key:result[key] for key in ('scoped_screens_passed', 'qualification_accepted') if key in result}, flush=True)


if __name__ == '__main__':
    main()
