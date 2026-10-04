"""Extend the validated phase-carry pair through history-9 formation.

Preserve physical state and carry from elapsed 6; keep old failures immutable.
"""
import argparse
import fcntl
import json
from pathlib import Path
import shutil

import numpy as np
import torch

from .feedback_long import digest, restore_checkpoint
from .geometry_precision import worker, checked_history, exact_state
from .gpu_backend import library as accepted_library
from .gpu_precision_control import library, PrecisionSimulation
from .neighbor_context import read
from .parameter_robustness import verify
from .parameter_robustness_moving import endpoint_assay
from .polarity_robustness import refinement_comparison, first_crossing
from .resolution import write_json


def transfer(source, destination, old_hash, new_hash, job, horizon):
    """Re-encode experiment metadata, preserving every physical/residual byte."""
    host, audit, history = restore_checkpoint(source, old_hash, job)
    if abs(host.time-150.-horizon) > 1e-9 or history[-1]['elapsed'] != horizon:
        raise ValueError('Wrong transferred physical clock')
    with np.load(source, allow_pickle=False) as z:
        payload = {key:z[key].copy() for key in z.files}
    if (str(payload['precision_arm']) != 'phase_carry' or payload['phase_carry'].dtype != np.float64
            or payload['phase_carry'].shape != host.phi.shape or not np.isfinite(payload['phase_carry']).all()):
        raise ValueError('Invalid transferred phase carry')
    payload['long_experiment'] = np.array(json.dumps(dict(audit=audit, history=history, protocol_hash=new_hash, job=job)))
    temporary = destination.with_name(destination.stem+'.tmp.npz')
    np.savez_compressed(temporary, **payload); temporary.replace(destination)
    with np.load(destination, allow_pickle=False) as z:
        for key, value in payload.items():
            if not np.array_equal(z[key], value): raise ValueError('Transfer changed physical checkpoint bytes')
    return audit, history


def prepare(root, short=Path('outputs/geometry-precision'), old=Path('outputs/polarity-robustness-refined')):
    root, short, old = [Path(x).resolve() for x in (root, short, old)]
    if root.exists(): raise FileExistsError(root)
    p = read(short/'protocol.json'); verify(p)
    long = read(old/'protocol.json'); verify(long)
    if read(short/'status.json')['state'] != 'completed' or not read(short/'summary.json')['quality_pass']:
        raise ValueError('Accepted short precision controls required')
    inputs = [short/f for f in ('protocol.json','summary.json','status.json','implementation-gate.json')]
    inputs += [old/f for f in ('protocol.json','long-refinement.json','summary.json','status.json')]
    jobs = [j for j in p['jobs'] if j['arm'] == 'phase_carry']
    for job in jobs:
        folder = short/job['key']; r = read(folder/'result.json')
        if (r['protocol_sha256'] != digest(short/'protocol.json') or r['job'] != job or not r['quality_pass']
                or digest(folder/'history.json') != r['history_sha256'] or
                digest(folder/'latest_state.npz') != r['checkpoint_sha256']):
            raise ValueError('Changed completed short carry trajectory')
        checked_history(read(folder/'history.json'), p, job, p['duration'])
        inputs += [folder/f for f in ('result.json','history.json','latest_state.npz')]
    references = dict(fine=long['reference_folder'], finer=str(old/long['jobs'][0]['key']))
    for folder in references.values():
        folder = Path(folder); r = read(folder/'result.json')
        if not r['quality_pass'] or digest(folder/'history.json') != r['history_sha256']:
            raise ValueError('Changed long baseline')
        inputs += [folder/f for f in ('history.json','result.json','latest_state.npz')]
    root.mkdir(parents=True)
    sources = {**p['source_sha256'], **long['source_sha256'], str(Path(__file__).resolve()):digest(__file__)}
    protocol = dict(p, jobs=jobs, duration=240., pilot_duration=60., pilot_late_window=12., late_window=24.,
        criteria=long['criteria'], refinement_criteria=long['refinement_criteria'],
        endpoint_horizons=long['endpoint_horizons'], endpoint_interval=long['endpoint_interval'],
        endpoint_criteria=long['endpoint_criteria'], short=str(short), baseline_folders=references,
        reuse_duration=6., estimated_new_seconds=sum(read(short/j['key']/'result.json')['audit']['wall_seconds']
            *(240./6.-1) for j in jobs),
        design='Two carry-preserving continuations of the unchanged history-9 near-uniform t=150 restart at dt=.001875 and .0009375. Reuse original elapsed 0–6 without resetting residual, chemistry, geometry, polarity, lineage or random streams; then continue to elapsed 240 (physical t=390). Contact accumulation stays float32. No new histories or zygotes.',
        interpretation='Test whether phase-update carry preserves initiation and improves full formation-time agreement under the original .01 raw chemical-error criterion and all inherited limits. Compare with both immutable no-carry baselines at matching timesteps. Passing a two-timestep pair is not full temporal/spatial convergence, new-backend scientific acceptance, or autonomous identity. Earlier failures remain failures.',
        source_sha256=sources, input_sha256={**p['input_sha256'], **{str(f.resolve()):digest(f) for f in inputs}})
    write_json(root/'protocol.json', protocol); ph = digest(root/'protocol.json')
    transfers = []
    for job in jobs:
        folder = root/job['key']; folder.mkdir()
        audit, history = transfer(short/job['key']/'latest_state.npz', folder/'latest_state.npz',
            digest(short/'protocol.json'), ph, job, p['duration'])
        write_json(folder/'history.json', history)
        write_json(folder/'status.json', dict(state='prepared', elapsed=6., until=240., audit=audit))
        transfers.append(dict(job=job['key'], residual_preserved=True, reused_samples=len(history),
            source_checkpoint_sha256=digest(short/job['key']/'latest_state.npz'),
            transferred_checkpoint_sha256=digest(folder/'latest_state.npz')))
    write_json(root/'transfer-verification.json', dict(passed=True, transfers=transfers))
    write_json(root/'status.json', dict(state='prepared', completed=0, total=2, independent_histories=1,
        new_histories=0, estimated_new_seconds=protocol['estimated_new_seconds']))
    return protocol


def assess(root):
    root = Path(root); p = read(root/'protocol.json'); verify(p); ph = digest(root/'protocol.json')
    histories = {}; rows = []; endpoints = []
    for job in p['jobs']:
        folder = root/job['key']; r = read(folder/'result.json'); h = read(folder/'history.json')
        if (r['job'] != job or r['protocol_sha256'] != ph or not r['quality_pass'] or
                r['history_sha256'] != digest(folder/'history.json') or r['checkpoint_sha256'] != digest(folder/'latest_state.npz')):
            raise ValueError('Changed completed formation evidence')
        checked_history(h, p, job, p['duration'])
        host, audit, saved = restore_checkpoint(folder/'latest_state.npz', ph, job)
        if h != saved or audit != r['audit'] or host.time != 390.: raise ValueError('Checkpoint/history mismatch')
        histories[job['level']] = h
        t = np.array([x['elapsed'] for x in h]); spread = np.array([x['log_activator_sd'] for x in h])
        late_min = float(spread[t >= 216.-1e-9].min())
        rows.append(dict(**r, late_min_log_sd=late_min, persistent_contrast=late_min > .1,
            onset_elapsed=first_crossing(t, spread, .1, 'up')))
        e = endpoint_assay(root, job, p, h[-1]); endpoints.append(dict(level=job['level'],
            numerical_pass=e['numerical_pass'], all_trials_settled=e['all_trials_settled'],
            phase=e['phase'], local_bistability=e['local_bistability_supported']))
    pair = refinement_comparison(histories['fine'], histories['finer'], p, 240.)
    pilot = refinement_comparison([h for h in histories['fine'] if h['elapsed'] <= 60.+1e-9],
        [h for h in histories['finer'] if h['elapsed'] <= 60.+1e-9], p, 60.)
    baseline = {level:read(Path(folder)/'history.json') for level, folder in p['baseline_folders'].items()}
    unchanged_pair = refinement_comparison(baseline['fine'], baseline['finer'], p, 240.)
    effects = {level:refinement_comparison(baseline[level], histories[level], p, 240.) for level in ('fine','finer')}
    endpoint_pass = all(e['numerical_pass'] and e['all_trials_settled'] for e in endpoints)
    summary = dict(protocol_sha256=ph, completed=2, total=2, independent_histories=1, new_histories=0,
        quality_pass=True, pilot_refinement=pilot, full_refinement=pair, no_carry_full_refinement=unchanged_pair,
        same_dt_interventions=effects, results=rows, endpoints=endpoints, endpoint_numerical_pass=endpoint_pass,
        passed=bool(pair['passed'] and pilot['passed'] and endpoint_pass), scope=p['interpretation'])
    write_json(root/'summary.json', summary); return summary


def run(root):
    root = Path(root).resolve()
    with (root/'coordinator.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        p = read(root/'protocol.json'); verify(p); device = p['device']; torch.set_num_threads(1)
        if (device['name'] != torch.cuda.get_device_name(0) or device['torch'] != torch.__version__ or
                device['torch_cuda'] != torch.version.cuda or
                digest(library()._name) != device['experimental_binary_sha256'] or
                digest(accepted_library()._name) != device['accepted_binary_sha256']):
            raise ValueError('Changed GPU hardware/software')
        done = 0
        try:
            gate = root/'restart-gate.json'
            if not gate.exists():
                write_json(root/'status.json', dict(state='running', stage='carry_restart_gate', completed=0, total=2))
                checks = []
                for job in p['jobs']:
                    original = PrecisionSimulation.restore(Path(p['short'])/job['key']/'latest_state.npz')
                    reused = PrecisionSimulation.restore(root/job['key']/'latest_state.npz')
                    exact_state(original, reused)
                    for _ in range(2):
                        original.step(); reused.step(); exact_state(original, reused)
                        if not torch.equal(original.phase_carry, reused.phase_carry):
                            raise ValueError('Transferred carry changes the next physical step')
                    checks.append(dict(job=job['key'], bit_exact_steps=2, residual_preserved=True))
                    del original, reused
                torch.cuda.empty_cache()
                write_json(gate, dict(protocol_sha256=digest(root/'protocol.json'), passed=True, checks=checks))
            elif not read(gate)['passed'] or read(gate)['protocol_sha256'] != digest(root/'protocol.json'):
                raise ValueError('Changed continuation restart gate')
            for job in p['jobs']:
                write_json(root/'status.json', dict(state='running', current_job=job['key'], completed=done, total=2))
                worker(root, job, p); done += 1
            write_json(root/'status.json', dict(state='running', stage='endpoint_and_refinement', completed=2, total=2))
            summary = assess(root)
            write_json(root/'status.json', dict(state='completed' if summary['passed'] else 'completed_with_unresolved_checks',
                completed=2, total=2, passed=summary['passed'], quality_pass=True))
            return summary
        except Exception as error:
            write_json(root/'status.json', dict(state='failed', completed=done, total=2, error=str(error))); raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare','run','assess'))
    parser.add_argument('--output', type=Path, default=Path('outputs/phase-carry-formation')); args = parser.parse_args()
    report = prepare(args.output) if args.action == 'prepare' else run(args.output) if args.action == 'run' else assess(args.output)
    print({k:v for k,v in report.items() if k in ('completed','total','passed','estimated_new_seconds','full_refinement')})
