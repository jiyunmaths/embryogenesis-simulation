"""A separately registered timestep diagnosis of the stopped history-1103 pilot.

Reuse one immutable saved reference; execute one new full zygote trajectory.
The developmental solver, contact law, instrumentation and gates are unchanged.
"""
import argparse
from dataclasses import asdict, replace
from pathlib import Path
import fcntl
import hashlib
import json
import os
import time

import torch

from . import fresh_zygote_pilot as pilot
from . import development_control_qualification as controls
from . import development_time_diagnosis as td
from . import qualify_repeated_development as q
from .model import Config
from .qualify_development_gpu import digest, read, write, pin, save_cpu
from .resolution import _steps

PARENT = Path('output/studies/fresh-zygote-pilot')
ASSESSMENT = Path('docs/fresh_zygote_pilot_assessment.json')
DESIGN = Path('docs/fresh_zygote_time_refinement.md')
REFERENCE = 'h1103_off_development_half'
CANDIDATE = 'h1103_off_development_quarter'
REFERENCE_DT = .000234375
CANDIDATE_DT = .0001171875
PARENT_SHA = '728f584676ea4eeeee530656015aec61931e9296ff095eff8b7c94f9832a13ee'
INHERITED = ('duration', 'interval', 'checkpoint_interval', 'diagnostic_interval',
    'snapshots', 'time_limits', 'physical_limits', 'restart_steps', 'restart_stages',
    'minimum_mature_window', 'method', 'gpu_uuid', 'dense_time_spectrum_limit',
    'instrumentation_gate', 'instrumentation_protocol_sha256', 'runtime_binary_sha256',
    'prefix_steps', 'prefix_limits')


def candidate_config(reference):
    if (reference['seed'] != 1103 or reference['grid'] != 120 or
            reference['dt'] != REFERENCE_DT):
        raise ValueError('The recorded history-1103 reference is required')
    controls.config_changes(reference, reference, 'off')
    return asdict(replace(Config(**reference), dt=CANDIDATE_DT,
        steps=_steps(18., CANDIDATE_DT), save_every=_steps(.15, CANDIDATE_DT)))


def validate_config(reference, candidate):
    if candidate != candidate_config(reference):
        raise ValueError('Only the registered timestep and its integer clocks may change')


def signature(source, job):
    state = pilot.initial_signature(source, job)
    return hashlib.sha256(json.dumps(state, sort_keys=True).encode()).hexdigest()


def check_reference(root, p):
    job = next(j for j in p['jobs'] if j['key'] == REFERENCE)
    result = controls.case_report(root, job, p)
    if not result['passed']:
        raise ValueError('The immutable reference path must pass its physical gates')
    return result


def verify(p):
    q.verify(p)
    parent = read(p['parent_protocol'])
    if digest(p['parent_protocol']) != PARENT_SHA:
        raise ValueError('The original failed pilot protocol changed')
    if any(p[key] != parent[key] for key in INHERITED):
        raise ValueError('Inherited numerical limits or instrumentation changed')
    if p['new_scientific_histories'] != 0 or p['new_numerical_paths'] != 1:
        raise ValueError('Numerical repeats must not count as new independent histories')
    reference, candidate = p['jobs']
    if (reference['key'] != REFERENCE or candidate['key'] != CANDIDATE or
            not reference['reused'] or candidate['reused'] or
            candidate['seed'] != 1103 or candidate['branch'] != 'off'):
        raise ValueError('The registered single-history diagnostic design changed')
    validate_config(reference['config'], candidate['config'])
    starts = [signature(j['source'], j) for j in p['jobs']]
    if starts[0] != starts[1] or starts != [p['matched_start_signature']]*2:
        raise ValueError('The zygote state or random streams do not match the reference')
    gate = read(p['instrumentation_gate'])
    if not gate['passed'] or gate['protocol_sha256'] != p['instrumentation_protocol_sha256']:
        raise ValueError('Unchanged accepted instrumentation is required')
    for file, expected in gate['evidence_sha256'].items():
        if digest(file) != expected:
            raise ValueError('Instrumentation evidence changed: '+file)


def prepare(root):
    root = Path(root).resolve()
    if root.exists():
        raise FileExistsError('A fresh diagnostic directory is required')
    parent_path = PARENT/'protocol.json'
    if digest(parent_path) != PARENT_SHA:
        raise ValueError('Original failed-pilot protocol mismatch')
    parent = read(parent_path)
    q.verify(parent)
    assessment = read(ASSESSMENT)
    review_path = Path(assessment['independent_saved_evidence_review']['path'])
    if (not assessment['review_completed'] or assessment['qualification_accepted'] or
            assessment['protocol_sha256'] != PARENT_SHA or
            digest(review_path) != assessment['independent_saved_evidence_review']['sha256']):
        raise ValueError('Independently reviewed failed pilot is required')
    old_status = read(PARENT/'status.json')
    failed_pair = read(PARENT/'time-h1103-off.json')
    if (old_status['state'] != 'failed' or old_status['completed'] != 18 or
            failed_pair['passed'] or
            failed_pair['dense']['peak']['relative_spectrum'] <= parent['dense_time_spectrum_limit']):
        raise ValueError('The original stopped design must retain its numerical failure')
    original = next(j for j in parent['jobs'] if j['key'] == REFERENCE)
    folder = (PARENT/REFERENCE).resolve()
    report = read(folder/'gpu-result.json')
    if (not report['passed'] or report['job'] != original or
            report['protocol_sha256'] != PARENT_SHA):
        raise ValueError('Original reference provenance is inconsistent')
    for file, expected in report['evidence_sha256'].items():
        if digest(file) != expected:
            raise ValueError('Original reference evidence changed: '+file)
    config = candidate_config(original['config'])
    root.mkdir(parents=True)
    (root/REFERENCE).symlink_to(folder, target_is_directory=True)
    candidate_folder = root/CANDIDATE
    candidate_folder.mkdir()
    source = candidate_folder/'source.npz'
    host = q.initial(Config(**config))
    save_cpu(host, source)
    candidate = dict(key=CANDIDATE, seed=1103, branch='off', level='development_quarter',
        config=config, reused=False, source=str(source))
    reference = dict(original, reused=True)
    start = signature(source, candidate)
    if start != signature(original['source'], original):
        raise ValueError('New full-start arrays and RNG differ from the original zygote')
    inputs = dict(parent['input_sha256'])
    inputs.update(report['evidence_sha256'])
    inputs.update(pin([parent_path, PARENT/'status.json', PARENT/'time-h1103-off.json',
        folder/'gpu-result.json', folder/'actual-start-prefix.json',
        ASSESSMENT, review_path, source]))
    sources = dict(parent['source_sha256'])
    sources.update(pin([Path(__file__), DESIGN,
        Path('tests/test_fresh_zygote_time_refinement.py'),
        Path('archive/study-launches/2026-10-11-fresh-zygote-time-refinement/audit_prelaunch.py'),
        Path('archive/study-launches/2026-10-11-fresh-zygote-time-refinement/launch.py')]))
    p = {key:parent[key] for key in INHERITED}
    p.update(jobs=[reference, candidate], histories=[1103], branches=['off'],
        grid=120, new_scientific_histories=0, new_numerical_paths=1, reused_numerical_paths=1,
        matched_start_signature=start, source_sha256=sources, input_sha256=inputs,
        parent_protocol=str(parent_path.resolve()), parent_protocol_sha256=PARENT_SHA,
        original_failed_status_sha256=digest(PARENT/'status.json'),
        original_failed_pair_sha256=digest(PARENT/'time-h1103-off.json'),
        parent_assessment=str(ASSESSMENT.resolve()),
        reuse_provenance={REFERENCE:dict(folder=str(folder), job=original,
            report_sha256=digest(folder/'gpu-result.json'), protocol_sha256=PARENT_SHA)},
        time_pair=dict(seed=1103, branch='off', cases=[REFERENCE,CANDIDATE]),
        stop_policy='Stop on any prefix, physical, restart or timestep failure. No reruns, '
            'seed substitutions, threshold relaxation, contact-law changes or automatic extension.',
        preservation_scope='All original pinned sources/inputs, failed status and comparison, '
            'immutable reference report and its full saved evidence are verified. The original '
            'historical document manifest remains unmodified; mutable manuscript drafts are not '
            'added as simulation pins. Archived manuscript reconciliations remain in the parent review.',
        scope='Numerical diagnosis of existing history1103, feedback off, grid120, t0–18. '
            'One new full-start path at dt0.0001171875 versus the original dt0.000234375 path. '
            'No new independent history, full fine CPU qualification, original cohort acceptance, '
            'remaining mechanical controls, long continuation or biological identity claim.')
    write(root/'protocol.json', p)
    verify(p)
    check_reference(root, p)
    write(root/'status.json', dict(state='prepared', protocol_sha256=digest(root/'protocol.json')))
    return p


def validate_prelaunch(root, p):
    verify(p)
    check_reference(root, p)
    if os.environ.get('CUDA_VISIBLE_DEVICES') != p['gpu_uuid']:
        raise ValueError('The designated GPU must be selected explicitly')
    if read(Path(root)/'status.json') != dict(state='prepared',protocol_sha256=digest(Path(root)/'protocol.json')):
        raise ValueError('Unchanged prepared status is required')
    if (Path(root)/CANDIDATE/'gpu-result.json').exists():
        raise FileExistsError('Never rerun a recorded trajectory')


def run(root):
    root = Path(root).resolve()
    with (root/'coordinator.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
        p = read(root/'protocol.json')
        validate_prelaunch(root, p)
        torch.set_num_threads(1)
        lib = td.gpu.library()
        if digest(lib.binary_path) != p['runtime_binary_sha256']:
            raise ValueError('The qualified runtime binary changed')
        write(root/'runtime.json',dict(gpu=torch.cuda.get_device_name(0),torch=torch.__version__,
            binary=lib.binary_path,binary_sha256=digest(lib.binary_path)))
        started = time.perf_counter()
        job = next(j for j in p['jobs'] if j['key'] == CANDIDATE)
        paths = []
        prefix = None
        pair = None
        def status(state, **extra):
            write(root/'status.json',dict(state=state,active=CANDIDATE,completed=len(paths),
                total=1,new_scientific_histories=0,reused_paths=1,post_cleavage_started=False,
                qualification_accepted=False,seconds=time.perf_counter()-started,
                protocol_sha256=digest(root/'protocol.json'),**extra))
        try:
            status('running_prefix')
            prefix = pilot.prefix(root,job,p)
            if not prefix['passed']:
                raise RuntimeError('Actual-start prefix failed')
            status('running_development')
            report = controls.execute_path(root,job)
            paths.append(dict(case=CANDIDATE,passed=report['passed'],seconds=report['seconds']))
            if not report['passed']:
                raise RuntimeError('Full developmental path failed')
            status('comparing')
            pair = pilot.assess_pair(root,p['time_pair'],p)
            if not pair['passed']:
                raise RuntimeError('Registered timestep comparison failed')
            verify(p)
            check_reference(root,p)
            result = dict(execution_completed=True,developmental_screens_passed=True,
                qualification_accepted=False,independent_review_pending=True,
                original_pilot_remains_failed=True,post_cleavage_started=False,
                new_scientific_histories=0,new_numerical_paths=1,reused_numerical_paths=1,
                paths=paths,prefix=prefix,pair=pair,seconds=time.perf_counter()-started,
                protocol_sha256=digest(root/'protocol.json'),scope=p['scope'])
            write(root/'summary.json',result)
            status('completed_provisional',developmental_screens_passed=True,independent_review_pending=True)
            return result
        except Exception as exc:
            status('failed',failure=dict(type=type(exc).__name__,message=str(exc)),
                paths=paths,prefix=prefix,pair=pair,original_pilot_remains_failed=True)
            raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('prepare','run'))
    parser.add_argument('root',type=Path)
    args = parser.parse_args()
    (prepare if args.action == 'prepare' else run)(args.root)


if __name__ == '__main__':
    main()
