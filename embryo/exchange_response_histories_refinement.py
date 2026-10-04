"""Timestep halving of moving exchange formation and selected history responses.

Keeps previous sources/evidence immutable. A separate full-horizon GPU gate
reuses accepted dt=0.001875 native references before any scientific continuation.
Response starts are the original t=210 states, not the new formation endpoints.
"""
import argparse
import fcntl
import json
from pathlib import Path
import shutil

import numpy as np
import torch

from .attribute_development import AttributeSimulation
from .cell_exchange_response_moving import BACKGROUNDS, verify
from .cell_exchange_response_refinement import assess as response_assessment
from .cell_response_moving import name
from .cell_response_moving_refinement import CRITERIA as RESPONSE_CRITERIA
from .exchange_response_histories import (completed_history, formation_assessment,
                                         prefix_check, read)
from .feedback_long import digest
from .feedback_survival_validation import retime
from .gpu_response_runner import (compatible, prepare as gpu_prepare,
                                 require_validation, run as gpu_run)
from .resolution import write_json
from . import validate_gpu_backend as gpu_validation


RESPONSE_BACKGROUNDS = ('unexchanged', 'fresh_exchange')
FORMATION_CRITERIA = dict(chemical_log_max=.01, polarity_abs_max=.01,
                         relative_axis_max=.01, relative_volume_max=.005,
                         relative_transport_max=.01, final_phi_abs_max=.02,
                         pair_ratio_abs_max=.05, same_outcome_required=True)


def hashes(files):
    return {str(Path(f).resolve()): digest(f) for f in files}


def require_completed(root, expected=None, passed=False):
    status = read(Path(root)/'status.json')
    if (status.get('state') != 'completed' or
            (expected is not None and status.get('completed') != expected) or
            (passed and not status.get('passed'))):
        raise ValueError('Requires completed passing evidence: '+str(root))
    return status


def prepare_gate(root, cpu, coarse_gate, sources):
    """Canonical read-only CPU reference adapter for the unchanged GPU gate."""
    root, cpu = Path(root), Path(cpu)
    original = read(cpu/'protocol.json'); verify(original)
    require_completed(cpu, 6, passed=True)
    if not read(cpu/'refinement.json')['passed']:
        raise ValueError('CPU timestep refinement did not pass')
    old = read(Path(coarse_gate)/'protocol.json'); verify(old)
    require_completed(coarse_gate, 4, passed=True)
    accepted = AttributeSimulation.restore(Path(old['baseline'])/'unexchanged/source.npz').config
    if original['dt'] != old['dt']/2:
        raise ValueError('Expected one timestep halving')
    # Validate every selected reference before creating the adapter.
    jobs = [dict(family=f, target=t, factor=1. if t is None else .9)
            for f in RESPONSE_BACKGROUNDS for t in (None, 27)]
    inputs = [cpu/'protocol.json', cpu/'status.json', cpu/'refinement.json',
              Path(coarse_gate)/'protocol.json', Path(coarse_gate)/'status.json',
              Path(coarse_gate)/'comparison.json']
    times = np.arange(round(original['duration']/original['interval'])+1)*original['interval']
    for family in RESPONSE_BACKGROUNDS:
        child = cpu/family; cp = read(child/'protocol.json'); verify(cp)
        host = AttributeSimulation.restore(cp['checkpoint'])
        config = type(host.config)(**vars(host.config)); config.dt = accepted.dt
        compatible(config, accepted)
        if host.time != original['start'] or host.config.dt != original['dt'] or host.divisions:
            raise ValueError('Invalid fine CPU physical start')
        for job in (j for j in jobs if j['family'] == family):
            completed_history(child, job, times)
            inputs.extend(child/name(job)/f for f in ('result.json', 'history.json', 'latest_state.npz'))
        inputs.extend((child/'protocol.json', Path(cp['checkpoint']), child/'initial_states.npz'))
    root.mkdir()
    baseline = root/'cpu-references'; baseline.mkdir()
    for family in RESPONSE_BACKGROUNDS:
        child = baseline/family; child.mkdir(); source = cpu/family
        cp = read(source/'protocol.json')
        shutil.copy2(cp['checkpoint'], child/'source.npz')
        shutil.copy2(source/'initial_states.npz', child/'initial_states.npz')
        # Reference directories are never written by the validation runner.
        for job in (j for j in jobs if j['family'] == family):
            (child/name(job)).symlink_to((source/name(job)).resolve(), target_is_directory=True)
        inputs.extend((child/'source.npz', child/'initial_states.npz'))
    p = dict(baseline=str(baseline.resolve()), jobs=jobs, start=original['start'],
             dt=original['dt'], duration=original['duration'], interval=original['interval'],
             checkpoint_interval=3., criteria=gpu_validation.CRITERIA,
             scope='Four full same-state dt=0.001875 GPU/native comparisons, using accepted seed-7 untouched/fresh controls and cell-27 negative pulses. CPU references are existing completed fine trajectories, not rerun or fabricated.',
             limits='Backend agreement in the tested mature regime; not time/spatial convergence, cleavage validation, or acceptance of arbitrary parameters.',
             source_sha256=sources, input_sha256=hashes(inputs))
    write_json(root/'protocol.json', p)
    write_json(root/'reference-provenance.json', dict(cpu=str(cpu.resolve()),
        cpu_protocol_sha256=digest(cpu/'protocol.json'), references=[name(j) for j in jobs],
        note='Canonical path adapter to unchanged validated CPU results.'))
    write_json(root/'status.json', dict(state='prepared', completed=0, total=4))
    return p


def prepare_source(source, output, dt, jobs, sources, extra_files=()):
    source, output = Path(source), Path(output)
    cp = read(source/'protocol.json'); verify(cp)
    output.mkdir(parents=True)
    sim = retime(cp['checkpoint'], output/'source.npz', dt)
    if sim.time != cp['start'] or sim.divisions or len(sim.ids) != sim.config.max_cells:
        raise ValueError('Retiming changed physical start or allows cleavage')
    shutil.copy2(source/'initial_states.npz', output/'initial_states.npz')
    with np.load(output/'initial_states.npz') as initial:
        if not np.array_equal(sim.ids, initial['ids']):
            raise ValueError('Initial cell IDs changed')
        for job in jobs:
            state = initial[job['family']]
            if state.shape != (2, len(sim.ids)) or not np.isfinite(state).all() or np.any(state <= 0):
                raise ValueError('Invalid prepared chemistry')
            if job['target'] is not None and job['target'] not in sim.ids:
                raise ValueError('Missing prespecified target')
    files = [source/'protocol.json', Path(cp['checkpoint']), source/'initial_states.npz',
             output/'source.npz', output/'initial_states.npz', *extra_files]
    p = dict(checkpoint=str((output/'source.npz').resolve()), start=cp['start'], dt=dt,
             duration=cp['duration'], interval=cp['interval'], jobs=jobs,
             source_sha256=sources, input_sha256=hashes(files))
    write_json(output/'protocol.json', p)
    return files+[output/'protocol.json']


def prepare(root, baseline=Path('outputs/exchange-response-histories'),
            cpu=Path('outputs/cell-exchange-response-refined'),
            coarse_gate=Path('outputs/gpu-backend-validation')):
    root, baseline, cpu, coarse_gate = [Path(x).resolve() for x in (root, baseline, cpu, coarse_gate)]
    if root.exists():
        raise FileExistsError(root)
    old = read(baseline/'protocol.json'); verify(old); require_completed(baseline, 36)
    if not read(baseline/'comparison.json')['completed'] or old['seeds'] != [8, 9]:
        raise ValueError('Requires completed prespecified history cohort')
    fine = read(cpu/'protocol.json'); verify(fine)
    dt = old['dt']/2
    if fine['dt'] != dt:
        raise ValueError('Existing CPU references have another timestep')
    sources = {**old['source_sha256'], **fine['source_sha256'],
               **read(coarse_gate/'protocol.json')['source_sha256'],
               **hashes([Path(__file__), Path(__file__).with_name('cell_exchange_response_refinement.py')])}
    verify(dict(source_sha256=sources))
    # Check all coarse results before freezing any new protocol.
    for seed in old['seeds']:
        for stage in ('formation', 'response'):
            groups = [baseline/f'seed-{seed}'/stage] if stage == 'formation' else [
                baseline/f'seed-{seed}'/stage/f for f in RESPONSE_BACKGROUNDS]
            for group in groups:
                cp = read(group/'protocol.json')
                times = np.arange(round(cp['duration']/cp['interval'])+1)*cp['interval']
                for job in cp['jobs']:
                    if stage == 'formation' or job['target'] is None or job['factor'] == .9:
                        completed_history(group, job, times)
    root.mkdir(parents=True)
    prepare_gate(root/'gpu-validation', cpu, coarse_gate, sources)
    files = [baseline/f for f in ('protocol.json', 'status.json', 'comparison.json')]
    files.append(root/'gpu-validation/protocol.json')
    contexts = []
    for seed in old['seeds']:
        original = baseline/f'seed-{seed}'; dest = root/f'seed-{seed}'
        targets = next(h['targets'] for h in old['histories'] if h['seed'] == seed)
        formation_jobs = [dict(family=f, target=None, factor=1.) for f in BACKGROUNDS]
        source = dest/'formation-source'
        frozen = original/'formation-source/frozen_reference.npz'
        files += prepare_source(original/'formation', source, dt, formation_jobs, sources, [frozen])
        shutil.copy2(frozen, source/'frozen_reference.npz'); files.append(source/'frozen_reference.npz')
        for stage, families in (('formation', BACKGROUNDS), ('response', RESPONSE_BACKGROUNDS)):
            for family in families:
                if stage == 'response':
                    source = dest/'response-source'/family
                    jobs = [dict(family=family, target=None, factor=1.)]+[
                        dict(family=family, target=cell, factor=.9) for cell in targets]
                    files += prepare_source(original/'response'/family, source, dt, jobs, sources)
                else:
                    source = dest/'formation-source'
                # Reuse the existing strict prefix checker through a path adapter.
                alias = root/'backend-checks'/f'seed-{seed}-{stage}-{family}'
                parent = alias/f'seed-{seed}'; parent.mkdir(parents=True)
                (parent/'formation-source').symlink_to(source.resolve(), target_is_directory=True)
                contexts.append(dict(seed=seed, stage=stage, family=family, root=str(alias)))
        response = dest/'response'; response.mkdir()
        selected_files = []
        for family in RESPONSE_BACKGROUNDS:
            cp = read(original/'response'/family/'protocol.json')
            selected_files += [original/'response'/family/'protocol.json',
                              original/'response'/family/'initial_states.npz']
            for job in cp['jobs']:
                if job['target'] is None or job['factor'] == .9:
                    selected_files.extend(original/'response'/family/name(job)/f for f in ('result.json', 'history.json'))
        rp = dict(baseline=str(original/'response'), backgrounds=RESPONSE_BACKGROUNDS,
                  targets=targets, factor=.9, start=old['response_start'], duration=old['response_duration'],
                  interval=old['interval'], dt=dt, refinement_criteria=RESPONSE_CRITERIA,
                  classification_criteria=dict(reference_separation_min=.01, nearest_reference_unchanged=True),
                  scope='Same-state t=210 response timestep halving in one independently developed history; untouched/fresh controls and negative pulses in both prespecified targets.',
                  limits='Excludes positive and pre-relaxed pulses. Does not propagate fine formation endpoints into responses, test identity inheritance, or establish spatial convergence.',
                  source_sha256=sources, input_sha256=hashes(selected_files))
        write_json(response/'protocol.json', rp); files += selected_files+[response/'protocol.json']
        files += [original/'formation_comparison.json']
        for job in formation_jobs:
            files.extend(original/'formation'/name(job)/f for f in ('result.json', 'history.json', 'latest_state.npz'))
    p = dict(baseline=str(baseline), cpu_references=str(cpu), validation=str(root/'gpu-validation'),
             seeds=old['seeds'], histories=old['histories'], dt=dt, coarse_dt=old['dt'],
             formation_start=old['formation_start'], formation_duration=old['formation_duration'],
             response_start=old['response_start'], response_duration=old['response_duration'],
             interval=old['interval'], criteria=old['criteria'], formation_criteria=FORMATION_CRITERIA,
             response_criteria=RESPONSE_CRITERIA, prefix_duration=.6, prefix_native_threads=4,
             prefix_criteria=gpu_validation.CRITERIA, contexts=contexts, scientific_jobs=18,
             validation_jobs=4, prefix_checks=len(contexts),
             scope='One timestep halving, 0.00375 to 0.001875: six moving formation/retention continuations and twelve same-state negative-pulse/control continuations across seeds 8 and 9. Four full CPU/GPU replays and ten short new-context checks gate execution.',
             limits='Two additional histories, no new developmental replicates. Responses start from coarse-study t=210 physical states, separately from formation refinement. No positive/pre-relaxed response, end-to-end refined development, spatial convergence, or inherited/autonomous identity claim.',
             source_sha256=sources, input_sha256=hashes(files))
    write_json(root/'protocol.json', p)
    write_json(root/'status.json', dict(state='prepared', completed=0, total=18,
                                      prefix_completed=0, prefix_total=len(contexts)))
    return p


def formation_comparison(root, seed, p):
    fine = formation_assessment(root, seed, p)
    coarse = read(Path(p['baseline'])/f'seed-{seed}'/'formation_comparison.json')
    rows = []; times = np.arange(round(p['formation_duration']/p['interval'])+1)*p['interval']
    for family in BACKGROUNDS:
        job = dict(family=family, target=None, factor=1.)
        fine_root, coarse_root = Path(root)/f'seed-{seed}'/'formation', Path(p['baseline'])/f'seed-{seed}'/'formation'
        fh, fx = completed_history(fine_root, job, times)
        ch, cx = completed_history(coarse_root, job, times)
        if not np.array_equal(fx[0], cx[0]):
            raise ValueError('Formation refinement changed starting chemistry')
        errors = dict.fromkeys(('chemical_log_max', 'polarity_abs_max', 'relative_axis_max',
                               'relative_volume_max', 'relative_transport_max'), 0.)
        for a, b in zip(fh, ch):
            for key, value in gpu_validation.discrepancies(a, b).items(): errors[key] = max(errors[key], value)
        with np.load(fine_root/name(job)/'latest_state.npz') as a, np.load(coarse_root/name(job)/'latest_state.npz') as b:
            if (a['phi'].shape != b['phi'].shape or not np.isfinite(a['phi']).all()
                    or not np.isfinite(b['phi']).all()):
                raise ValueError('Incompatible or nonfinite endpoint phase fields')
            errors['final_phi_abs_max'] = float(abs(a['phi']-b['phi']).max())
        a = next(r for r in fine['trials'] if r['family'] == family)
        b = next(r for r in coarse['trials'] if r['family'] == family)
        outcomes = ('informative_pair', 'destination_like', 'transferred_like',
                    'contrast_retained', 'late_pair_order_reversed')
        same = all(a[k] == b[k] for k in outcomes)
        ratios = [abs(a[k]-b[k]) if a[k] is not None and b[k] is not None else
                  (0. if a[k] is None and b[k] is None else None)
                  for k in ('late_destination_ratio', 'late_transferred_ratio')]
        c = p['formation_criteria']
        passed = (same and all(np.isfinite(v) and v <= c[k] for k, v in errors.items())
                  and all(v is not None and v <= c['pair_ratio_abs_max'] for v in ratios))
        rows.append(dict(family=family, errors=errors, pair_ratio_errors=ratios,
                         outcomes_unchanged=same, fine=a, coarse=b, passed=bool(passed)))
    report = dict(seed=seed, passed=all(r['passed'] for r in rows), trials=rows,
                  criteria=p['formation_criteria'])
    write_json(Path(root)/f'seed-{seed}'/'formation_refinement.json', report)
    return report


def assess(root):
    root = Path(root).resolve(); p = read(root/'protocol.json'); verify(p)
    reports = []
    for seed in p['seeds']:
        formation = formation_comparison(root, seed, p)
        # Strong complete-history checks precede the inherited response assessor.
        for parent in (root, Path(p['baseline'])):
            for family in RESPONSE_BACKGROUNDS:
                child = parent/f'seed-{seed}'/'response'/family
                jobs = read(root/f'seed-{seed}'/'response-source'/family/'protocol.json')['jobs']
                times = np.arange(round(p['response_duration']/p['interval'])+1)*p['interval']
                for job in jobs: completed_history(child, job, times)
        responses = response_assessment(root/f'seed-{seed}'/'response')
        reports.append(dict(seed=seed, formation=formation, responses=responses,
                            passed=formation['passed'] and responses['passed']))
    result = dict(completed=True, passed=all(r['passed'] for r in reports), histories=reports,
                  protocol_sha256=digest(root/'protocol.json'), scope=p['scope'], limits=p['limits'])
    write_json(root/'refinement.json', result)
    return result


def run(root):
    root = Path(root).resolve()
    with (root/'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        p = read(root/'protocol.json'); verify(p); done = prefixes = 0
        torch.set_num_threads(1)
        def status(stage, **extra):
            write_json(root/'status.json', dict(state='running', stage=stage, completed=done,
                total=18, prefix_completed=prefixes, prefix_total=len(p['contexts']), **extra))
        try:
            status('fine_timestep_gpu_validation')
            gate = Path(p['validation'])
            if read(gate/'status.json').get('state') != 'completed': gpu_validation.run(gate)
            accepted, _ = require_validation(gate)
            if accepted.dt != p['dt']: raise ValueError('Wrong accepted GPU timestep')
            for context in p['contexts']:
                status('starting_state_backend_check', seed=context['seed'],
                       background=context['family'], context_stage=context['stage'])
                prefix_check(Path(context['root']), context['seed'], context['family'], p)
                prefixes += 1
            for seed in p['seeds']:
                source, output = root/f'seed-{seed}'/'formation-source', root/f'seed-{seed}'/'formation'
                status('formation', seed=seed)
                if not (output/'protocol.json').exists(): gpu_prepare(output, source, gate)
                gpu_run(output); done += 3
                formation_comparison(root, seed, p)
            for seed in p['seeds']:
                for family in RESPONSE_BACKGROUNDS:
                    source, output = root/f'seed-{seed}'/'response-source'/family, root/f'seed-{seed}'/'response'/family
                    status('response', seed=seed, background=family)
                    if not (output/'protocol.json').exists(): gpu_prepare(output, source, gate)
                    gpu_run(output); done += 3
                response_assessment(root/f'seed-{seed}'/'response')
            result = assess(root)
            write_json(root/'status.json', dict(state='completed', passed=result['passed'], completed=done,
                                               total=18, prefix_completed=prefixes, prefix_total=len(p['contexts'])))
        except Exception as exc:
            write_json(root/'status.json', dict(state='failed', completed=done, total=18,
                prefix_completed=prefixes, prefix_total=len(p['contexts']), error=str(exc)))
            raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('prepare', 'run', 'assess'))
    parser.add_argument('--output', type=Path, default=Path('outputs/exchange-response-histories-refined'))
    args = parser.parse_args()
    if args.command == 'prepare': prepare(args.output)
    elif args.command == 'run': run(args.output)
    else: assess(args.output)
