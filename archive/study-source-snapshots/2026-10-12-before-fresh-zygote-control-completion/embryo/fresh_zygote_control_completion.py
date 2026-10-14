"""Complete history 1103's matched controls at the accepted finer clocks.

Register six new full-start paths separately; retain the failed original pilot.
Reuse its independently accepted feedback-off refinement without executing it.
"""
import argparse
from dataclasses import asdict, replace
from pathlib import Path
import fcntl
import os
import time

import torch

from . import fresh_zygote_pilot as pilot
from . import fresh_zygote_time_refinement as refinement
from . import development_control_qualification as controls
from . import development_time_diagnosis as td
from . import qualify_repeated_development as q
from .model import Config
from .qualify_development_gpu import digest, read, write, pin, save_cpu
from .resolution import _steps

PARENT = Path('output/studies/fresh-zygote-time-refinement')
ASSESSMENT = Path('docs/fresh_zygote_time_refinement_assessment.json')
DESIGN = Path('docs/fresh_zygote_control_completion.md')
LAUNCH = Path('archive/study-launches/2026-10-12-fresh-zygote-control-completion')
PARENT_SHA = '54632a1f284affcbcec93cd0f1036d2d9c3788bd93ef5480df5fa1f14336c6c3'
BRANCHES = ('tension_adhesion', 'polarity', 'full')
LEVELS = (('development_half', .000234375), ('development_quarter', .0001171875))
REFERENCES = (refinement.REFERENCE, refinement.CANDIDATE)


def design(base):
    if base['seed'] != 1103 or base['grid'] != 120 or base['dt'] != LEVELS[0][1]:
        raise ValueError('Accepted history-1103 feedback-off reference is required')
    controls.config_changes(base, base, 'off')
    jobs = []
    for branch in BRANCHES:
        feedback, tension, adhesion, polarity = controls.COUPLINGS[branch]
        for level, dt in LEVELS:
            config = asdict(replace(Config(**base), dt=dt, steps=_steps(18., dt),
                save_every=_steps(.15, dt), feedback=feedback,
                fate_tension=tension, fate_adhesion=adhesion, polarity_tension=polarity))
            jobs.append(dict(key=f'h1103_{branch}_{level}', seed=1103,
                branch=branch, level=level, config=config, reused=False))
    return jobs


def validate_jobs(root, p):
    reference = next(j for j in p['jobs'] if j['key'] == REFERENCES[0])
    expected = [dict(j, source=str(Path(root)/j['key']/'source.npz'))
        for j in design(reference['config'])]
    actual = [j for j in p['jobs'] if not j['reused']]
    pairs = [dict(seed=1103, branch=branch,
        cases=[f'h1103_{branch}_{level}' for level, _ in LEVELS]) for branch in BRANCHES]
    if actual != expected or p['time_pairs'] != pairs:
        raise ValueError('Registered controls, model or clocks changed')
    if ([j['key'] for j in p['jobs'] if j['reused']] != list(REFERENCES) or
            len(p['jobs']) != 8 or p['histories'] != [1103] or
            p['branches'] != list(BRANCHES) or p['new_scientific_histories'] != 0 or
            p['new_numerical_paths'] != 6 or p['reused_numerical_paths'] != 2):
        raise ValueError('Matched controls are nested within one existing history')


def check_references(root, p):
    reports = []
    for key in REFERENCES:
        job = next(j for j in p['jobs'] if j['key'] == key)
        report = controls.case_report(root, job, p)
        if not report['passed']:
            raise ValueError('Reused feedback-off path failed its physical gates')
        reports.append(report)
    return reports


def verify(root, p):
    q.verify(p)
    parent = read(PARENT/'protocol.json')
    if digest(PARENT/'protocol.json') != PARENT_SHA:
        raise ValueError('Accepted refinement protocol changed')
    if any(p[key] != parent[key] for key in refinement.INHERITED):
        raise ValueError('Inherited limits or instrumentation changed')
    validate_jobs(root, p)
    starts = [refinement.signature(j['source'], j) for j in p['jobs']]
    if starts != [p['matched_start_signature']]*8:
        raise ValueError('Zygote arrays, cell-cycle schedule or random streams differ')
    accepted = read(ASSESSMENT)
    review = read(accepted['independent_review']['path'])
    if (not accepted['review_completed'] or not accepted['scoped_refinement_accepted'] or
            accepted['original_cohort_accepted'] or accepted['protocol_sha256'] != PARENT_SHA or
            digest(accepted['independent_review']['path']) != accepted['independent_review']['sha256'] or
            not review['review_completed'] or not review['scoped_refinement_accepted']):
        raise ValueError('Independent feedback-off refinement acceptance is required')
    original = refinement.PARENT
    if (read(original/'status.json')['state'] != 'failed' or
            read(original/'time-h1103-off.json')['passed']):
        raise ValueError('Original failed pilot must remain failed')
    gate = read(p['instrumentation_gate'])
    if not gate['passed'] or gate['protocol_sha256'] != p['instrumentation_protocol_sha256']:
        raise ValueError('Accepted unchanged instrumentation is required')
    for name, expected in gate['evidence_sha256'].items():
        if digest(name) != expected:
            raise ValueError('Instrumentation evidence changed: '+name)


def prepare(root):
    root = Path(root).resolve()
    if root.exists():
        raise FileExistsError('A fresh matched-control directory is required')
    if digest(PARENT/'protocol.json') != PARENT_SHA:
        raise ValueError('Accepted refinement protocol mismatch')
    parent = read(PARENT/'protocol.json')
    refinement.verify(parent)
    refinement.check_reference(PARENT, parent)
    accepted = read(ASSESSMENT)
    review_path = Path(accepted['independent_review']['path'])
    if (not accepted['review_completed'] or not accepted['scoped_refinement_accepted'] or
            accepted['original_cohort_accepted'] or
            digest(review_path) != accepted['independent_review']['sha256']):
        raise ValueError('Independently accepted finer feedback-off pair required')
    base = next(j for j in parent['jobs'] if j['key'] == REFERENCES[0])['config']
    new_jobs = design(base)
    root.mkdir(parents=True)
    jobs = []; reuse = {}; inputs = dict(parent['input_sha256'])
    inputs.update(pin([PARENT/'protocol.json', PARENT/'status.json', PARENT/'summary.json',
        PARENT/'runtime.json', PARENT/'time-h1103-off.json', ASSESSMENT, review_path]))
    for original in parent['jobs']:
        folder = (PARENT/original['key']).resolve()
        report = read(folder/'gpu-result.json')
        if not report['passed'] or report['job']['config'] != original['config']:
            raise ValueError('Accepted feedback-off report changed')
        (root/original['key']).symlink_to(folder, target_is_directory=True)
        jobs.append(dict(original, reused=True))
        reuse[original['key']] = dict(folder=str(folder), job=report['job'],
            report_sha256=digest(folder/'gpu-result.json'), protocol_sha256=report['protocol_sha256'])
        inputs.update(report['evidence_sha256'])
        inputs.update(pin([folder/'gpu-result.json', folder/'actual-start-prefix.json']))
    for job in new_jobs:
        folder = root/job['key']; folder.mkdir()
        source = folder/'source.npz'
        save_cpu(q.initial(Config(**job['config'])), source)
        job['source'] = str(source); jobs.append(job)
        inputs.update(pin([source]))
    start = refinement.signature(jobs[0]['source'], jobs[0])
    if any(refinement.signature(j['source'], j) != start for j in jobs):
        raise ValueError('Complete initial state and RNG matching failed')
    sources = dict(parent['source_sha256'])
    sources.update(pin([Path(__file__), DESIGN, Path('tests/test_fresh_zygote_control_completion.py'),
        LAUNCH/'audit_prelaunch.py', LAUNCH/'launch.py']))
    p = {key:parent[key] for key in refinement.INHERITED}
    p.update(jobs=jobs, histories=[1103], branches=list(BRANCHES), grid=120,
        levels=[dict(level=level, dt=dt) for level, dt in LEVELS],
        new_scientific_histories=0, new_numerical_paths=6, reused_numerical_paths=2,
        matched_start_signature=start, source_sha256=sources, input_sha256=inputs,
        parent_protocol_sha256=PARENT_SHA, accepted_assessment=str(ASSESSMENT.resolve()),
        reuse_provenance=reuse,
        time_pairs=[dict(seed=1103, branch=branch,
            cases=[f'h1103_{branch}_{level}' for level, _ in LEVELS]) for branch in BRANCHES],
        stop_policy='Stop on any actual-start, physical, live restart or timestep gate failure. '
            'No reruns, seed substitutions, limit relaxation, contact-law changes or long extension.',
        scope='Six new full-start matched-control paths within existing history1103 at grid120 '
            'through t18, two finer timesteps. Two accepted feedback-off references reused read-only. '
            'No new independent history, original pilot acceptance, full fine CPU qualification, '
            'long formation/maintenance continuation or biological identity claim.')
    write(root/'protocol.json', p)
    verify(root, p); check_references(root, p)
    write(root/'status.json', dict(state='prepared', protocol_sha256=digest(root/'protocol.json')))
    return p


def validate_prelaunch(root, p):
    verify(root, p); check_references(root, p)
    if os.environ.get('CUDA_VISIBLE_DEVICES') != p['gpu_uuid']:
        raise ValueError('Select the dedicated GPU explicitly')
    if read(Path(root)/'status.json') != dict(state='prepared', protocol_sha256=digest(Path(root)/'protocol.json')):
        raise ValueError('Unchanged prepared status required')
    for job in p['jobs']:
        if not job['reused'] and (Path(root)/job['key']/'gpu-result.json').exists():
            raise FileExistsError('Never rerun a recorded trajectory')


def complete_design(pairs):
    return (len(pairs) == 3 and {(r['seed'], r['branch']) for r in pairs} ==
        {(1103, b) for b in BRANCHES} and all(r['passed'] for r in pairs))


def run(root):
    root = Path(root).resolve()
    with (root/'coordinator.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
        p = read(root/'protocol.json'); validate_prelaunch(root, p)
        torch.set_num_threads(1); lib = td.gpu.library()
        if digest(lib.binary_path) != p['runtime_binary_sha256']:
            raise ValueError('Qualified runtime binary changed')
        write(root/'runtime.json', dict(gpu=torch.cuda.get_device_name(0), torch=torch.__version__,
            binary=lib.binary_path, binary_sha256=digest(lib.binary_path)))
        started = time.perf_counter(); paths = []; pairs = []; prefixes = []; active = None
        def status(state, **extra):
            write(root/'status.json', dict(state=state, active=active, completed=len(paths), total=6,
                pairs_completed=len(pairs), pairs_total=3, new_scientific_histories=0, reused_paths=2,
                post_cleavage_started=False, qualification_accepted=False,
                original_pilot_remains_failed=True, seconds=time.perf_counter()-started,
                protocol_sha256=digest(root/'protocol.json'), **extra))
        try:
            for pair in p['time_pairs']:
                for case in pair['cases']:
                    job = next(j for j in p['jobs'] if j['key'] == case); active = case
                    if job['reused']:
                        raise ValueError('Read-only references cannot be executed')
                    status('running_prefix')
                    prefix = pilot.prefix(root, job, p); prefixes.append(prefix)
                    if not prefix['passed']:
                        raise RuntimeError('Actual-start prefix failed: '+case)
                    status('running_development')
                    report = controls.execute_path(root, job)
                    paths.append(dict(case=case, passed=report['passed'], seconds=report['seconds']))
                    if not report['passed']:
                        raise RuntimeError('Developmental path failed: '+case)
                status('comparing')
                result = pilot.assess_pair(root, pair, p); pairs.append(result)
                if not result['passed']:
                    raise RuntimeError('Registered timestep comparison failed: '+pair['branch'])
            verify(root, p); check_references(root, p)
            if not complete_design(pairs):
                raise RuntimeError('Incomplete matched-control design')
            result = dict(execution_completed=True, developmental_screens_passed=True,
                qualification_accepted=False, independent_review_pending=True,
                original_pilot_remains_failed=True, post_cleavage_started=False,
                new_scientific_histories=0, new_numerical_paths=6, reused_numerical_paths=2,
                paths=paths, pairs=pairs, prefixes=prefixes, seconds=time.perf_counter()-started,
                protocol_sha256=digest(root/'protocol.json'), scope=p['scope'])
            write(root/'summary.json', result)
            status('completed_provisional', developmental_screens_passed=True, independent_review_pending=True)
            return result
        except Exception as exc:
            status('failed', failure=dict(type=type(exc).__name__, message=str(exc)),
                paths=paths, pairs=pairs, prefixes=prefixes)
            raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'run'))
    parser.add_argument('root', type=Path)
    args = parser.parse_args()
    (prepare if args.action == 'prepare' else run)(args.root)


if __name__ == '__main__':
    main()
