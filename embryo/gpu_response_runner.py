"""Validation-gated GPU execution of mature moving-response job protocols.

Prepares a separate output directory; accepted CPU study files are never edited.
Only the validated physical parameters and supported mature regime are admitted.
"""
import argparse
import fcntl
import json
from pathlib import Path
import shutil
import tempfile
import time

import numpy as np
import torch

from .attribute_development import AttributeSimulation
from .cell_exchange_response_moving import verify
from .cell_response_moving import initialize, name
from .feedback_long import digest, save_checkpoint
from .gpu_backend import GpuSimulation, library
from .resolution import write_json, _steps
from .validate_gpu_backend import assess

PARAMETERS = ('grid', 'extent', 'dt', 'max_cells', 'interface_width', 'surface_tension',
              'volume_stiffness', 'adhesion', 'repulsion', 'fate_tension', 'fate_adhesion',
              'signal_transport', 'signaling', 'signal_beta', 'signal_da', 'signal_dh',
              'graph_contact_cutoff', 'polarity_enabled', 'polarity_rate', 'polarity_alignment',
              'polarity_decay', 'polarity_tension', 'polarity_contact_cutoff')


def require_validation(root):
    root = Path(root).resolve()
    status = json.loads((root/'status.json').read_text())
    if status.get('state') != 'completed' or not status.get('passed') or status.get('completed') != 4:
        raise ValueError('Full four-run GPU validation must complete and pass first')
    p = json.loads((root/'protocol.json').read_text())
    verify(p)
    saved = json.loads((root/'comparison.json').read_text())
    if not saved.get('scientific_ready') or saved.get('protocol_sha256') != digest(root/'protocol.json'):
        raise ValueError('GPU acceptance evidence is missing or invalid')
    # Reassess without rewriting the accepted comparison artifact.
    with tempfile.TemporaryDirectory(prefix='embryo-gpu-gate-') as temporary:
        check = Path(temporary)
        shutil.copy2(root/'protocol.json', check/'protocol.json')
        for job in p['jobs']:
            (check/name(job)).symlink_to(root/name(job), target_is_directory=True)
        result = assess(check)
    if not result['passed'] or not result['scientific_ready']:
        raise ValueError('GPU agreement no longer passes')
    if result != saved:
        raise ValueError('GPU acceptance report differs from current evidence')
    device = json.loads((root/'device.json').read_text())
    if (device['name'] != torch.cuda.get_device_name(0) or device['torch'] != torch.__version__ or
            device['torch_cuda'] != torch.version.cuda or device['binary_sha256'] != digest(library()._name)):
        raise ValueError('GPU software/hardware changed; repeat validation')
    checkpoint = Path(p['baseline'])/'unexchanged/source.npz'
    accepted = AttributeSimulation.restore(checkpoint).config
    evidence = [root/f for f in ('protocol.json', 'comparison.json', 'status.json', 'device.json')]
    for job in p['jobs']:
        evidence.extend(root/name(job)/f for f in ('result.json', 'history.json', 'latest_state.npz'))
    return accepted, evidence


def compatible(config, accepted):
    changed = [key for key in PARAMETERS if getattr(config, key) != getattr(accepted, key)]
    if changed:
        raise ValueError('Parameters outside validated GPU regime: '+', '.join(changed))


def prepare(output, source, validation):
    output, source = Path(output).resolve(), Path(source).resolve()
    if output.exists():
        raise FileExistsError(output)
    accepted, evidence = require_validation(validation)
    original = json.loads((source/'protocol.json').read_text())
    verify(original)
    checkpoint = Path(original['checkpoint'])
    with np.load(checkpoint) as data:
        if 'fertilization_cue' in data.files:
            raise ValueError('Cue forcing requires a separately validated GPU implementation')
    host = AttributeSimulation.restore(checkpoint)
    compatible(host.config, accepted)
    if abs(host.time-original['start']) > 1e-9 or host.config.dt != original['dt']:
        raise ValueError('Source clock mismatch')
    steps = _steps(original['duration'], original['dt'])
    every = _steps(original['interval'], original['dt'])
    if steps % every or _steps(3., original['dt']) % every:
        raise ValueError('Observation, checkpoint and endpoint clocks must align')
    if not original['jobs'] or len({name(j) for j in original['jobs']}) != len(original['jobs']):
        raise ValueError('Expected nonempty distinct jobs')
    # Reject unsupported regimes before writing any output.
    GpuSimulation(host).audit()
    with np.load(source/'initial_states.npz') as data:
        if not np.array_equal(data['ids'], host.ids):
            raise ValueError('Initial cell IDs mismatch')
        if not np.allclose(data['masses'], host.volumes(), rtol=1e-12, atol=1e-14):
            raise ValueError('Initial compartment masses mismatch')
        for job in original['jobs']:
            initialize(checkpoint, data[job['family']], job['target'], job['factor'])
    output.mkdir(parents=True)
    shutil.copy2(checkpoint, output/'source.npz')
    shutil.copy2(source/'initial_states.npz', output/'initial_states.npz')
    code = [Path(__file__)]
    inputs = [source/'protocol.json', checkpoint, source/'initial_states.npz',
              output/'source.npz', output/'initial_states.npz', *evidence]
    p = {**original, 'checkpoint': str(output/'source.npz'), 'backend': 'resident-gpu',
         'gpu_validation': str(Path(validation).resolve()), 'original_protocol': str(source/'protocol.json'),
         'original_protocol_sha256': digest(source/'protocol.json'),
         'source_sha256': {**original.get('source_sha256', {}),
                          **{str(f.resolve()): digest(f) for f in code}},
         'input_sha256': {**original.get('input_sha256', {}), **{str(f): digest(f) for f in inputs}}}
    write_json(output/'protocol.json', p)
    write_json(output/'status.json', dict(state='prepared', completed=0, total=len(p['jobs']), backend=p['backend']))
    return p


def run(output):
    output = Path(output).resolve()
    with (output/'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        p = json.loads((output/'protocol.json').read_text())
        verify(p)
        accepted, _ = require_validation(p['gpu_validation'])
        ph = digest(output/'protocol.json')
        torch.set_num_threads(1)
        done = 0
        try:
            for job in p['jobs']:
                folder = output/name(job)
                folder.mkdir(exist_ok=True)
                checkpoint = folder/'latest_state.npz'
                if (folder/'result.json').exists():
                    result = json.loads((folder/'result.json').read_text())
                    if not result['quality_pass'] or result['protocol_sha256'] != ph:
                        raise ValueError('Invalid saved GPU result')
                    done += 1
                    continue
                if checkpoint.exists():
                    with np.load(checkpoint) as data:
                        meta = json.loads(str(data['long_experiment']))
                    if meta['protocol_hash'] != ph or meta['job'] != job:
                        raise ValueError('GPU restart mismatch')
                    sim = GpuSimulation.restore(checkpoint)
                    audit, history = meta['audit'], meta['history']
                else:
                    with np.load(output/'initial_states.npz') as data:
                        host = initialize(p['checkpoint'], data[job['family']], job['target'], job['factor'])
                    sim = GpuSimulation(host)
                    audit = dict(max_volume_error=0., min_radius=1e100, max_clipping=0.,
                                 max_sampled_boundary=0., dilution_amount_error=0., elapsed_seconds=0.)
                    history = []
                compatible(sim.config, accepted)
                step0 = 0 if p['start'] == 0 else _steps(p['start'], p['dt'])
                stop = step0+_steps(p['duration'], p['dt'])
                every, checkpoint_every = _steps(p['interval'], p['dt']), _steps(3., p['dt'])
                last = history[-1]['time'] if history else -1.
                started, previous_wall = time.perf_counter(), audit['elapsed_seconds']
                write_json(output/'status.json', dict(state='running', completed=done, total=len(p['jobs']), current=name(job), backend='resident-gpu'))
                while sim.step_number <= stop:
                    for key, value in sim.audit().items():
                        audit[key] = min(audit[key], value) if key == 'min_radius' else max(audit[key], value)
                    if audit['dilution_amount_error'] > 2e-14:
                        raise RuntimeError('GPU dilution conservation failure')
                    if (sim.step_number-step0)%every == 0 and sim.time > last+1e-9:
                        row = sim.observe(sim.time-p['start'])
                        audit['max_sampled_boundary'] = max(audit['max_sampled_boundary'], row['boundary_occupancy'])
                        if audit['max_sampled_boundary'] >= .01:
                            raise RuntimeError('GPU boundary screen failed')
                        history.append(row); last = row['time']
                        audit['elapsed_seconds'] = previous_wall+time.perf_counter()-started
                        write_json(folder/'history.json', history)
                        write_json(folder/'status.json', dict(state='running', elapsed=row['elapsed'], audit=audit))
                        if (sim.step_number-step0)%checkpoint_every == 0 or sim.step_number == stop:
                            save_checkpoint(sim, checkpoint, audit, history, ph, job)
                            print(name(job), row['elapsed'], flush=True)
                    if sim.step_number == stop:
                        break
                    sim.step()
                write_json(folder/'result.json', dict(job=job, quality_pass=True, protocol_sha256=ph,
                    gpu_validation_protocol_sha256=digest(Path(p['gpu_validation'])/'protocol.json'), audit=audit))
                write_json(folder/'status.json', dict(state='completed', elapsed=p['duration'], audit=audit))
                done += 1
            write_json(output/'status.json', dict(state='completed', completed=done, total=len(p['jobs']), backend='resident-gpu'))
        except Exception as exc:
            write_json(output/'status.json', dict(state='failed', completed=done, total=len(p['jobs']), error=str(exc)))
            raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('prepare', 'run'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--source', type=Path)
    parser.add_argument('--validation', type=Path, default=Path('outputs/gpu-backend-validation'))
    args = parser.parse_args()
    if args.command == 'prepare':
        if args.source is None:
            parser.error('prepare requires --source child-study directory')
        prepare(args.output, args.source, args.validation)
    else:
        run(args.output)
