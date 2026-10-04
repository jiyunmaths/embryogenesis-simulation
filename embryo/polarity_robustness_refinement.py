"""Targeted second timestep halving of history-9 zero-contrast initiation.

Preserve the completed parent and reuse its dt=.001875 trajectory read-only.
Gate one dt=.0009375 restart on exact-context native/GPU and 60-unit agreement.
The inherited quantitative tolerances remain unchanged, including failures.
"""
import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import fcntl
from pathlib import Path
import shutil

import numpy as np

from .attribute_development import AttributeSimulation
from .feedback_long import digest, restore_checkpoint
from .feedback_survival_validation import retime
from .gpu_response_runner import compatible, require_validation
from .neighbor_context import read
from .parameter_robustness import verify
from .parameter_robustness_moving import initialize, prefix
from .polarity_robustness import (REFINEMENT_CRITERIA, advance,
    assert_same_physical_start, checked_history, job_view, refinement_comparison)
from .resolution import _steps, write_json


def check_halving(original, refined):
    """Allow only the timestep and its physically equivalent bookkeeping."""
    if refined.config.dt != original.config.dt/2:
        raise ValueError('Requires exactly one timestep halving')
    assert_same_physical_start(original, refined)
    if not np.array_equal(original.fate, refined.fate) or original.divisions != refined.divisions:
        raise ValueError('Retiming changed auxiliary state or divisions')
    old, new = asdict(original.config), asdict(refined.config)
    if any(old[k] != new[k] for k in old if k not in ('dt', 'steps', 'save_every')):
        raise ValueError('Retiming changed physical parameters')
    if (refined.step_number != _steps(original.time, refined.config.dt) or
            refined.config.steps != _steps(old['steps']*old['dt'], new['dt']) or
            refined.config.save_every != _steps(old['save_every']*old['dt'], new['dt'])):
        raise ValueError('Retiming changed physical observation or endpoint clocks')


def parent_evidence(parent, p, job):
    """Validate a completed reference including its actual physical checkpoint."""
    folder = parent/job['key']; ph = digest(parent/'protocol.json')
    r = read(folder/'result.json'); h = read(folder/'history.json')
    if (r['job'] != job or not r['quality_pass'] or r['protocol_sha256'] != ph or
            r['history_sha256'] != digest(folder/'history.json') or
            r['checkpoint_sha256'] != digest(folder/'latest_state.npz')):
        raise ValueError('Changed completed parent reference')
    checked_history(h, p, job, p['duration'])
    host, audit, saved = restore_checkpoint(folder/'latest_state.npz', ph, job)
    if (saved != h or audit != r['audit'] or abs(host.time-p['start']-p['duration']) > 1e-9 or
            not np.array_equal(np.array([host.activator, host.inhibitor]), h[-1]['chemistry'])):
        raise ValueError('Parent checkpoint and observations disagree')
    cp = read(folder/'prefix/comparison.json')
    if cp['protocol_sha256'] != ph or cp['job'] != job or not cp['passed']:
        raise ValueError('Parent backend context gate failed')
    for f, expected in cp['evidence_sha256'].items():
        if digest(f) != expected: raise ValueError('Changed parent backend evidence')
    ep = read(folder/'endpoint/protocol.json'); verify(ep)
    e = read(folder/'endpoint/assay/result.json')
    if (not e['numerical_pass'] or not e['all_trials_settled'] or
            e['protocol_sha256'] != digest(folder/'endpoint/protocol.json') or
            e['paths_sha256'] != digest(folder/'endpoint/assay/paths.npz') or
            ep['moving_protocol_sha256'] != ph or ep['moving_history_sha256'] != r['history_sha256']):
        raise ValueError('Changed parent endpoint evidence')
    files = [folder/f for f in ('result.json', 'history.json', 'latest_state.npz',
        'prefix/comparison.json', 'endpoint/protocol.json', 'endpoint/assay/result.json',
        'endpoint/assay/paths.npz')]
    files += [Path(f) for f in cp['evidence_sha256']]
    files += [Path(f) for f in ep['input_sha256']]
    return r, h, files


def prepare(root, parent=Path('outputs/polarity-robustness')):
    root, parent = Path(root).resolve(), Path(parent).resolve()
    if root.exists(): raise FileExistsError(root)
    old = read(parent/'protocol.json'); verify(old)
    status = read(parent/'status.json')
    if (status['state'] != 'completed_with_unresolved_checks' or status['completed'] != 21 or
            old['duration'] != 240. or old['refinement_criteria'] != REFINEMENT_CRITERIA):
        raise ValueError('Requires the completed parent with its unchanged unresolved gate')
    selected = [next(j for j in old['jobs'] if j['seed'] == 9 and j['family'] == 'uniform' and
                    j['point']['chi'] == 0. and j['level'] == level) for level in ('coarse', 'fine')]
    references = [parent_evidence(parent, old, j) for j in selected]
    prior = refinement_comparison(references[0][1], references[1][1], old, old['duration'])
    stored = read(parent/'long-refinement.json')
    if (stored['protocol_sha256'] != digest(parent/'protocol.json') or
            next(r for r in stored['comparisons'] if r['chi'] == 0.) != dict(chi=0., **prior) or
            prior['passed'] or prior['errors']['chemical_log_max'] <= REFINEMENT_CRITERIA['chemical_log_max']):
        raise ValueError('Parent failure no longer matches the prespecified target')
    validation = old['validation_roots']['fine']
    accepted, gate_files = require_validation(validation)
    if asdict(accepted) != old['accepted_configs']['fine'] or accepted.dt != .001875:
        raise ValueError('Changed accepted fine backend baseline')
    original = AttributeSimulation.restore(selected[1]['source']); compatible(original.config, accepted)
    root.mkdir(parents=True); folder = root/'seed-9'; folder.mkdir()
    source, chemicals = folder/'finer-source.npz', folder/'initial_states.npz'
    refined = retime(selected[1]['source'], source, accepted.dt/2)
    check_halving(original, refined)
    shutil.copy2(selected[1]['chemical_file'], chemicals)
    job = dict(key='seed-9_finer_polarity-0_uniform', seed=9, level='finer', dt=refined.config.dt,
        point=selected[1]['point'], family='uniform', source=str(source), chemical_file=str(chemicals))
    # Same physical/chemical start, after applying the identical intervention.
    a = initialize(selected[1]['source'], selected[1]['chemical_file'], 'uniform', job['point'], asdict(accepted))
    b = initialize(source, chemicals, 'uniform', job['point'], asdict(refined.config))
    check_halving(a, b)
    start_check = folder/'physical-start-verification.json'
    write_json(start_check, dict(physical_state_and_random_streams_equal=True,
        identical_initial_chemical_file=digest(chemicals) == digest(selected[1]['chemical_file']),
        fine_dt=accepted.dt, finer_dt=refined.config.dt, physical_time=a.time,
        fine_step=a.step_number, finer_step=b.step_number,
        scope='Retiming only; no restart from an evolved patterned endpoint.'))
    inputs = [parent/'protocol.json', parent/'status.json', parent/'pilot-refinement.json',
        parent/'long-refinement.json', source, chemicals, start_check,
        Path(selected[1]['source']), Path(selected[1]['chemical_file']), *gate_files]
    for _, _, files in references: inputs += files
    sources = {**old['source_sha256'], **{
        str(Path(__file__).with_name(f).resolve()): digest(Path(__file__).with_name(f))
        for f in ('polarity_robustness_refinement.py', 'resolution.py', 'cell_response_exchange.py')}}
    p = dict(parent=str(parent), parent_protocol_sha256=digest(parent/'protocol.json'),
        reference_job=selected[1], reference_folder=str(parent/selected[1]['key']),
        previous_coarse_fine_comparison=prior, validation_root=validation,
        validated_fine_config=asdict(accepted), accepted_configs=dict(finer=asdict(refined.config)),
        jobs=[job], total_jobs=1, independent_histories=1, new_independent_histories=0,
        start=old['start'], duration=old['duration'], pilot_duration=old['pilot_duration'],
        interval=old['interval'], checkpoint_interval=old['checkpoint_interval'],
        late_window=old['late_window'], pilot_late_window=old['pilot_late_window'],
        prefix_duration=old['prefix_duration'], prefix_native_threads=old['prefix_native_threads'],
        prefix_criteria=old['prefix_criteria'], criteria=old['criteria'],
        refinement_criteria=old['refinement_criteria'], endpoint_criteria=old['endpoint_criteria'],
        endpoint_horizons=old['endpoint_horizons'], endpoint_interval=old['endpoint_interval'],
        estimated_scientific_seconds=2*references[1][0]['audit']['wall_seconds'],
        stages=['retained full baseline backend gate at dt=.001875',
                'exact-context .6-unit native/GPU gate at dt=.0009375',
                'one finer restart to 60, compared to immutable fine reference',
                'continue that same physical checkpoint to 240 only after pilot passes',
                'frozen endpoint assay and full fine/finer comparison with unchanged tolerances'],
        design='Same t=150 history-9 mature geometry, polarity, lineage, random streams and exact near-uniform chemistry; beta=2, D_a=.02, D_b=.55, chi=0; activity-tension/adhesion unchanged. Restart from the original physical start, not from elapsed 60 or the patterned endpoint. Reuse completed dt=.001875 reference read-only. One new dt=.0009375 trajectory, ending at t=390.',
        backend_scope='Retain the four full native/GPU baseline replays at dt=.001875 and their hardware/software/binary hashes. New dt=.0009375 acceptance is an explicitly scoped .6-unit exact-context check with unchanged strict tolerances, not full-horizon native/GPU agreement at that timestep. Use resident PyTorch matrices/arrays and unchanged custom CUDA mechanics, geometry and polarity.',
        interpretation='The target was chosen to resolve a recorded formation-transient numerical failure, not to add an independent success history. Fine/finer agreement cannot erase the earlier coarse/fine failure or establish continuum/spatial convergence. Quantitative failure remains a limitation even when outcomes agree.',
        source_sha256=sources, input_sha256={str(f.resolve()): digest(f) for f in inputs})
    write_json(root/'protocol.json', p)
    write_json(root/'status.json', dict(state='prepared', stage='backend_context_gate',
        completed=0, total=1, independent_histories=1, new_independent_histories=0))
    return p


def compare(root, p, horizon):
    root = Path(root); ph = digest(root/'protocol.json')
    pilot = horizon == p['pilot_duration']
    if horizon not in (p['pilot_duration'], p['duration']): raise ValueError('Undeclared comparison horizon')
    target = root/('pilot-refinement.json' if pilot else 'long-refinement.json')
    if target.exists():
        report = read(target)
        if report['protocol_sha256'] != ph or report['horizon'] != horizon:
            raise ValueError('Changed refinement protocol')
        for f, expected in report['source_histories_sha256'].items():
            if digest(f) != expected: raise ValueError('Changed refinement evidence')
        paths = list(report['source_histories_sha256'])
        fine, finer = [[r for r in read(f) if r['elapsed'] <= horizon+1e-9] for f in paths]
        checked_history(fine, p, p['reference_job'], horizon)
        checked_history(finer, p, p['jobs'][0], horizon)
        recomputed = refinement_comparison(fine, finer, p, horizon)
        if any(report[k] != recomputed[k] for k in recomputed if k != 'interpretation'):
            raise ValueError('Changed refinement decision')
        return report
    job = p['jobs'][0]; folder = root/job['key']; reference = Path(p['reference_folder'])/'history.json'
    fine = [r for r in read(reference) if r['elapsed'] <= horizon+1e-9]
    full = read(folder/'history.json'); finer = [r for r in full if r['elapsed'] <= horizon+1e-9]
    checked_history(fine, p, p['reference_job'], horizon); checked_history(finer, p, job, horizon)
    host, _, saved = restore_checkpoint(folder/'latest_state.npz', ph, job)
    if saved != full or host.time < p['start']+horizon-1e-9:
        raise ValueError('Refined checkpoint and observations disagree')
    source = folder/'pilot-history.json' if pilot else folder/'history.json'
    if pilot:
        if source.exists() and read(source) != finer: raise ValueError('Changed immutable pilot observations')
        if not source.exists(): write_json(source, finer)
    result = refinement_comparison(fine, finer, p, horizon)
    result['interpretation'] = 'History-9 chi=0 exact-start fine/finer timestep halving; frozen modal growth is not full moving-system stability.'
    report = dict(**result, protocol_sha256=ph, compared_timesteps=[p['reference_job']['dt'], job['dt']],
        source_histories_sha256={str(f.resolve()): digest(f) for f in (reference, source)})
    write_json(target, report)
    return report


def assess(root):
    root = Path(root); p = read(root/'protocol.json'); verify(p)
    job = p['jobs'][0]; folder = root/job['key']
    result = read(folder/'result.json'); history = read(folder/'history.json')
    if (not result['quality_pass'] or result['protocol_sha256'] != digest(root/'protocol.json') or
            result['history_sha256'] != digest(folder/'history.json') or
            result['checkpoint_sha256'] != digest(folder/'latest_state.npz')):
        raise ValueError('Changed completed refined evidence')
    checked_history(history, p, job, p['duration'])
    gate = prefix(root, job, job_view(p, job)); pilot = compare(root, p, p['pilot_duration'])
    long = compare(root, p, p['duration'])
    from .parameter_robustness_moving import endpoint_assay
    endpoint = endpoint_assay(root, job, job_view(p, job), history[-1])
    fine = read(Path(p['reference_folder'])/'history.json')
    t = np.array([r['elapsed'] for r in history])
    errors = abs(np.log(np.array([r['chemistry'] for r in fine])/
                        np.array([r['chemistry'] for r in history]))).max(axis=(1, 2))
    previous = p['previous_coarse_fine_comparison']['errors']['chemical_log_max']
    report = dict(passed=bool(gate['passed'] and pilot['passed'] and long['passed'] and
        endpoint['numerical_pass'] and endpoint['all_trials_settled']),
        completed=1, total=1, independent_histories=1, new_independent_histories=0,
        protocol_sha256=digest(root/'protocol.json'), assessed_utc=datetime.now(timezone.utc).isoformat(),
        prefix_pass=gate['passed'], pilot_refinement_pass=pilot['passed'], long_refinement_pass=long['passed'],
        fine_finer=long, previous_coarse_fine=p['previous_coarse_fine_comparison'],
        maximum_error_fraction_of_previous=float(errors.max()/previous),
        late_chemical_log_max=float(errors[t >= p['duration']-p['late_window']-1e-9].max()),
        peak_chemical_error_elapsed=float(t[errors.argmax()]),
        errors_above_tolerance=int(np.sum(errors > p['refinement_criteria']['chemical_log_max'])),
        final_log_sd=history[-1]['log_activator_sd'], late_min_log_sd=result['late_min_log_activator_sd'],
        moving_quality_pass=True, audit=result['audit'], endpoint_numerical_pass=endpoint['numerical_pass'],
        endpoint_phase=endpoint['phase'], endpoint_local_bistability=endpoint['local_bistability_supported'],
        scope=p['interpretation'])
    write_json(root/'summary.json', report)
    return report


def run_stages(root, p):
    job = p['jobs'][0]
    write_json(root/'status.json', dict(state='running', stage='backend_context_gate',
        current_job=job['key'], completed=0, total=1))
    try:
        prefix(root, job, job_view(p, job))
    except RuntimeError as error:
        report = dict(state='stopped_at_backend_gate', error=str(error), completed=0, total=1)
        write_json(root/'status.json', report); return report
    write_json(root/'status.json', dict(state='running', stage='timestep_pilot',
        current_job=job['key'], completed=0, total=1, prefix_pass=True))
    advance(root, job, p, p['pilot_duration'])
    pilot = compare(root, p, p['pilot_duration'])
    if not pilot['passed']:
        report = dict(state='stopped_at_numerical_gate', stage='timestep_pilot',
            completed=0, total=1, pilot_refinement_pass=False)
        write_json(root/'status.json', report); return report
    write_json(root/'status.json', dict(state='running', stage='long_continuation',
        current_job=job['key'], completed=0, total=1, prefix_pass=True, pilot_refinement_pass=True))
    advance(root, job, p, p['duration'])
    report = assess(root)
    write_json(root/'status.json', dict(state='completed' if report['passed'] else 'completed_with_unresolved_checks',
        completed=1, total=1, prefix_pass=report['prefix_pass'],
        pilot_refinement_pass=report['pilot_refinement_pass'], long_refinement_pass=report['long_refinement_pass'],
        endpoint_numerical_pass=report['endpoint_numerical_pass']))
    return report


def run(root):
    import torch
    root = Path(root).resolve()
    with (root/'coordinator.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        p = read(root/'protocol.json'); verify(p); torch.set_num_threads(1)
        try:
            accepted, _ = require_validation(p['validation_root'])
            if asdict(accepted) != p['validated_fine_config']:
                raise ValueError('Changed baseline backend validation')
            original = AttributeSimulation.restore(p['reference_job']['source'])
            refined = AttributeSimulation.restore(p['jobs'][0]['source']); check_halving(original, refined)
            return run_stages(root, p)
        except Exception as error:
            write_json(root/'status.json', dict(state='failed', error=str(error), completed=0, total=1))
            raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'run', 'assess', 'status'))
    parser.add_argument('--output', type=Path, default=Path('outputs/polarity-robustness-refined'))
    parser.add_argument('--parent', type=Path, default=Path('outputs/polarity-robustness'))
    args = parser.parse_args()
    result = (prepare(args.output, args.parent) if args.action == 'prepare' else
              run(args.output) if args.action == 'run' else
              assess(args.output) if args.action == 'assess' else read(args.output/'status.json'))
    print({k: v for k, v in result.items() if k in ('state', 'stage', 'completed', 'total',
        'total_jobs', 'independent_histories', 'new_independent_histories', 'passed',
        'prefix_pass', 'pilot_refinement_pass', 'long_refinement_pass', 'estimated_scientific_seconds', 'error')})
