"""Full-horizon resident-GPU gate against accepted native CPU trajectories."""
import argparse
import fcntl
import json
from pathlib import Path
import time

import numpy as np
import torch

from .attribute_development import AttributeSimulation
from .cell_exchange_response_moving import verify
from .cell_response import pulse, response_metrics
from .cell_response_moving import name
from .cell_response_moving_refinement import response_error
from .feedback_long import digest, save_checkpoint
from .gpu_backend import GpuSimulation, library
from .resolution import write_json, _steps


CRITERIA = dict(chemical_log_max=1e-5, polarity_abs_max=1e-5, relative_axis_max=1e-4,
                relative_volume_max=1e-5, relative_transport_max=1e-5, final_phi_abs_max=2e-5,
                normalized_response_max=.001, relative_auc_max=.002,
                recovery_time_error_max=.15, dilution_amount_error=2e-14)


def prepare(root, baseline=Path('outputs/cell-exchange-response-moving')):
    root, baseline = Path(root).resolve(), Path(baseline).resolve()
    if root.exists():
        raise FileExistsError(root)
    p = json.loads((baseline/'protocol.json').read_text())
    verify(p)
    if json.loads((baseline/'status.json').read_text())['state'] != 'completed':
        raise ValueError('Requires completed CPU study')
    manifest = baseline/'native-backend-transition.json'
    m = json.loads(manifest.read_text())
    for f, h in {**m['source_sha256'], **m['validation_sha256'], **m['completed_results']}.items():
        if digest(f) != h:
            raise ValueError('Native provenance changed: '+f)
    inputs = [baseline/'protocol.json', baseline/'status.json', manifest,
              baseline/'comparison.json']
    jobs = []
    for background in ('unexchanged', 'fresh_exchange'):
        child = baseline/background
        cp = json.loads((child/'protocol.json').read_text())
        verify(cp)
        sim = AttributeSimulation.restore(child/'source.npz')
        if sim.config.dt != p['dt'] or sim.time != p['start'] or sim.divisions:
            raise ValueError('Starting state mismatch')
        inputs.extend(child/f for f in ('protocol.json', 'source.npz', 'initial_states.npz'))
        # The same recipient samples low- and high-activity contexts before/after exchange.
        for cell, factor in ((None, 1.), (27, .9)):
            job = dict(family=background, target=cell, factor=factor)
            folder = child/name(job)
            result = json.loads((folder/'result.json').read_text())
            if not result['quality_pass'] or result['protocol_sha256'] != digest(child/'protocol.json'):
                raise ValueError('Invalid CPU reference')
            inputs.extend(folder/f for f in ('result.json', 'history.json', 'latest_state.npz'))
            jobs.append(job)
    sources = [Path(__file__), *[Path(__file__).with_name(f) for f in (
        'gpu_backend.py', 'gpu_kernels.cu', 'cuda_mechanics.cu', 'cuda_spatial_bench.cu',
        'attribute_development.py', 'model.py', 'polarity.py', 'signaling.py', 'transport.py',
        'cell_response.py', 'cell_response_moving.py', 'cell_response_moving_refinement.py',
        'cell_exchange_response_moving.py', 'feedback_long.py', 'resolution.py')]]
    root.mkdir(parents=True)
    protocol = dict(baseline=str(baseline), jobs=jobs, start=p['start'], dt=p['dt'],
        duration=60., interval=p['interval'], checkpoint_interval=3., criteria=CRITERIA,
        scope='Four complete t=210-to-270 continuations: untouched/fresh-exchange moving controls and -10% activator pulses in recipient 27. Reference mature 16-cell, 72-cubed, direct-feedback conservative model; both low/high chemical contexts. Resident float32 fields, float64 mechanics/chemistry, PyTorch matrices and custom CUDA geometry/polarity. No cleavage, no prescribed fate, no full-field per-step CPU transfers.',
        limits='Backend agreement on one developmental history, two backgrounds and one negative pulse target at dt=.00375. Does not validate cleavage, cue forcing, no-feedback/nonpolar branches, spatial convergence, all parameter regimes, or autonomous identity.',
        source_sha256={str(f.resolve()): digest(f) for f in sources},
        input_sha256={str(f): digest(f) for f in inputs})
    write_json(root/'protocol.json', protocol)
    write_json(root/'status.json', dict(state='prepared', completed=0, total=len(jobs)))
    return protocol


def discrepancies(row, reference):
    if row['ids'] != reference['ids'] or abs(row['time']-reference['time']) > 1e-9:
        raise ValueError('Reference time/cell alignment mismatch')
    for value in (row, reference):
        if not all(np.isfinite(np.asarray(value[k])).all() for k in (
                'time', 'chemistry', 'polarity', 'axis_ratio', 'volumes', 'delta')):
            raise ValueError('Nonfinite observation')
        if np.any(np.asarray(value['chemistry']) <= 0) or np.any(np.asarray(value['volumes']) <= 0) or value['axis_ratio'] <= 0:
            raise ValueError('Nonpositive observation')
    if any(np.shape(row[k]) != np.shape(reference[k]) for k in ('chemistry', 'polarity', 'volumes', 'delta')):
        raise ValueError('Observation shape mismatch')
    current = np.array(row['delta'])
    old = np.array(reference['delta'])
    return dict(chemical_log_max=float(abs(np.log(np.array(row['chemistry'])/reference['chemistry'])).max()),
                polarity_abs_max=float(abs(np.array(row['polarity'])-reference['polarity']).max()),
                relative_axis_max=abs(row['axis_ratio']/reference['axis_ratio']-1),
                relative_volume_max=float(abs(np.array(row['volumes'])/reference['volumes']-1).max()),
                relative_transport_max=float(np.linalg.norm(current-old)/max(np.linalg.norm(old), 1e-300)))


def assess(root):
    root = Path(root)
    p = json.loads((root/'protocol.json').read_text())
    verify(p)
    ph = digest(root/'protocol.json')
    trials, runs = [], []
    for family in ('unexchanged', 'fresh_exchange'):
        selected = [j for j in p['jobs'] if j['family'] == family]
        gpu, cpu = [], []
        for job in selected:
            folder = root/name(job)
            result = json.loads((folder/'result.json').read_text())
            if not result['passed'] or result['protocol_sha256'] != ph:
                raise ValueError('Invalid GPU result')
            histories = [json.loads((parent/'history.json').read_text()) for parent in (
                folder, Path(p['baseline'])/family/name(job))]
            for history in histories:
                expected = np.arange(_steps(p['duration'], p['interval'])+1)*p['interval']
                if len(history) != len(expected) or not np.allclose([r['elapsed'] for r in history], expected, rtol=0, atol=1e-9):
                    raise ValueError('Incomplete or misaligned trajectory')
            maxima = dict.fromkeys(('chemical_log_max', 'polarity_abs_max', 'relative_axis_max',
                                   'relative_volume_max', 'relative_transport_max'), 0.)
            for a, b in zip(*histories):
                for key, value in discrepancies(a, b).items():
                    maxima[key] = max(maxima[key], value)
            # Independently assess recorded trajectories, not just the worker's decision.
            field_paths = (folder/'latest_state.npz', Path(p['baseline'])/family/name(job)/'latest_state.npz')
            with np.load(field_paths[0]) as x, np.load(field_paths[1]) as y:
                phi_error = float(abs(x['phi']-y['phi']).max())
            c = p['criteria']
            audit = result['audit']
            quality = (audit['max_volume_error'] < .05 and audit['min_radius'] >= 4 and
                       audit['max_clipping'] == 0 and audit['max_boundary'] < .01 and
                       audit['dilution_amount_error'] <= c['dilution_amount_error'])
            valid = all(np.isfinite(value) and value <= c[key] for key, value in maxima.items())
            result = {**result, 'recomputed_errors': maxima, 'final_phi_abs_max': phi_error,
                      'passed': bool(quality and valid and np.isfinite(phi_error) and phi_error <= c['final_phi_abs_max'])}
            runs.append(result)
            gpu.append(np.array([r['chemistry'] for r in histories[0]]))
            cpu.append(np.array([r['chemistry'] for r in histories[1]]))
        with np.load(Path(p['baseline'])/family/'initial_states.npz') as data:
            masses, ids = data['masses'], data['ids']
        times = np.arange(len(gpu[0]))*p['interval']
        cell = list(ids).index(27)
        a = response_metrics(gpu[1], gpu[0], times, masses, cell, .9)
        b = response_metrics(cpu[1], cpu[0], times, masses, cell, .9)
        error = response_error(cpu[1], cpu[0], gpu[1], gpu[0], .9)
        auc = abs(a['target_activator_log_auc_per_log_pulse']/b['target_activator_log_auc_per_log_pulse']-1)
        recovery = {}
        for key in ('target_recovery_time', 'network_recovery_time'):
            recovery[key] = 0. if a[key] is None and b[key] is None else (
                None if a[key] is None or b[key] is None else abs(a[key]-b[key]))
        c = p['criteria']
        passed = error <= c['normalized_response_max'] and auc <= c['relative_auc_max'] and all(
            v is not None and v <= c['recovery_time_error_max']+1e-9 for v in recovery.values())
        trials.append(dict(background=family, normalized_response_error=error, relative_auc_error=auc,
                           recovery_time_errors=recovery, gpu_metrics=a, cpu_metrics=b, passed=bool(passed)))
    passed = all(r['passed'] for r in trials+runs)
    result = dict(passed=passed, scientific_ready=passed, protocol_sha256=ph, runs=runs,
                  responses=trials, criteria=p['criteria'], scope=p['scope'], limits=p['limits'])
    write_json(root/'comparison.json', result)
    return result


def run(root, device='cuda:0'):
    root = Path(root).resolve()
    with (root/'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        p = json.loads((root/'protocol.json').read_text())
        verify(p)
        ph = digest(root/'protocol.json')
        torch.set_num_threads(1)
        library()
        write_json(root/'device.json', dict(name=torch.cuda.get_device_name(device),
                    torch=torch.__version__, torch_cuda=torch.version.cuda,
                    binary_sha256=digest(library()._name)))
        done = 0
        try:
            for job in p['jobs']:
                folder = root/name(job)
                folder.mkdir(exist_ok=True)
                source = Path(p['baseline'])/job['family']
                checkpoint = folder/'latest_state.npz'
                if (folder/'result.json').exists():
                    result = json.loads((folder/'result.json').read_text())
                    if not result['passed'] or result['protocol_sha256'] != ph:
                        raise ValueError('Invalid saved result')
                    done += 1
                    continue
                reference = json.loads((source/name(job)/'history.json').read_text())
                if checkpoint.exists():
                    with np.load(checkpoint) as data:
                        meta = json.loads(str(data['long_experiment']))
                    if meta['protocol_hash'] != ph or meta['job'] != job:
                        raise ValueError('Restart mismatch')
                    sim = GpuSimulation.restore(checkpoint, device)
                    audit, history = meta['audit'], meta['history']
                else:
                    host = AttributeSimulation.restore(source/'source.npz')
                    if job['target'] is not None:
                        host.activator, host.inhibitor = pulse(np.array([host.activator, host.inhibitor]),
                                                             list(host.ids).index(job['target']), job['factor'])
                    sim = GpuSimulation(host, device)
                    audit = dict.fromkeys(('chemical_log_max', 'polarity_abs_max', 'relative_axis_max',
                        'relative_volume_max', 'relative_transport_max', 'max_volume_error', 'max_clipping',
                        'max_boundary', 'dilution_amount_error', 'wall_seconds'), 0.)
                    audit['min_radius'] = 1e100
                    history = []
                start_clock, previous_wall = time.perf_counter(), audit['wall_seconds']
                step0 = _steps(p['start'], p['dt'])
                stop = step0+_steps(p['duration'], p['dt'])
                every = _steps(p['interval'], p['dt'])
                save_every = _steps(p['checkpoint_interval'], p['dt'])
                last = history[-1]['time'] if history else -1.
                write_json(root/'status.json', dict(state='running', completed=done, total=len(p['jobs']), current=name(job)))
                while sim.step_number <= stop:
                    quality = sim.audit()
                    for key, value in quality.items():
                        audit[key] = min(audit[key], value) if key == 'min_radius' else max(audit[key], value)
                    if audit['dilution_amount_error'] > p['criteria']['dilution_amount_error']:
                        raise RuntimeError('Concentration dilution changed amount')
                    if (sim.step_number-step0)%every == 0 and sim.time > last+1e-9:
                        row = sim.observe(sim.time-p['start'])
                        errors = discrepancies(row, reference[(sim.step_number-step0)//every])
                        for key, value in errors.items():
                            audit[key] = max(audit[key], value)
                        audit['max_boundary'] = max(audit['max_boundary'], row['boundary_occupancy'])
                        if any(audit[k] > p['criteria'][k] for k in errors) or audit['max_boundary'] >= .01:
                            raise RuntimeError('GPU reference agreement gate exceeded: '+str(audit))
                        history.append(row)
                        last = row['time']
                        audit['wall_seconds'] = previous_wall+time.perf_counter()-start_clock
                        write_json(folder/'history.json', history)
                        write_json(folder/'status.json', dict(state='running', elapsed=row['elapsed'], audit=audit))
                        if (sim.step_number-step0)%save_every == 0 or sim.step_number == stop:
                            save_checkpoint(sim, checkpoint, audit, history, ph, job)
                            print(name(job), row['elapsed'], audit['chemical_log_max'], flush=True)
                    if sim.step_number == stop:
                        break
                    sim.step()
                with np.load(source/name(job)/'latest_state.npz') as data:
                    phi_error = float(abs(sim.phi.cpu().numpy()-data['phi']).max())
                if phi_error > p['criteria']['final_phi_abs_max']:
                    raise RuntimeError('Final GPU phase fields disagree: '+str(phi_error))
                write_json(folder/'result.json', dict(job=job, passed=True, protocol_sha256=ph,
                                                      audit=audit, final_phi_abs_max=phi_error))
                write_json(folder/'status.json', dict(state='completed', elapsed=p['duration'], audit=audit))
                done += 1
            result = assess(root)
            write_json(root/'status.json', dict(state='completed', passed=result['passed'], completed=done, total=len(p['jobs'])))
        except Exception as exc:
            write_json(root/'status.json', dict(state='failed', completed=done, total=len(p['jobs']), error=str(exc)))
            raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('prepare', 'run', 'assess'))
    parser.add_argument('--output', type=Path, default=Path('outputs/gpu-backend-validation'))
    parser.add_argument('--device', default='cuda:0')
    args = parser.parse_args()
    if args.command == 'prepare':
        prepare(args.output)
    elif args.command == 'run':
        run(args.output, args.device)
    else:
        assess(args.output)
