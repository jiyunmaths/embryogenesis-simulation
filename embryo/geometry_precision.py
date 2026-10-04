"""Short factorial moving-geometry precision diagnosis, not backend acceptance."""
import argparse
from dataclasses import asdict
import fcntl
from pathlib import Path
import time

import numpy as np
import torch

from .attribute_development import AttributeSimulation
from .feedback_long import digest, save_checkpoint, restore_checkpoint
from .feedback_survival_validation import retime
from .gpu_backend import GpuSimulation, library as accepted_library
from .gpu_precision_control import PrecisionSimulation, ARMS, library
from .gpu_response_runner import require_validation
from .neighbor_context import read, validate_graph
from .parameter_robustness import verify
from .parameter_robustness_moving import initialize, retain
from .polarity_robustness import assert_same_physical_start
from .polarity_robustness_refinement import check_halving
from .resolution import write_json, _steps

DTS = (.001875, .0009375)


def exact_state(a, b):
    for key in ('phi', 'activator', 'inhibitor', 'polarity', 'geometry'):
        if not torch.equal(getattr(a, key), getattr(b, key)):
            raise ValueError('New baseline differs from accepted kernel: '+key)
    if a.time != b.time or a.step_number != b.step_number or not torch.equal(a.matrices()[1], b.matrices()[1]):
        raise ValueError('New baseline differs in clock or transport')


def implementation_gate(root, source, chemical, point, config):
    """Exact experimental context, not an extrapolated synthetic-only check."""
    make = lambda: initialize(source, chemical, 'uniform', point, config)
    original = GpuSimulation(make()); baseline = PrecisionSimulation(make())
    carried = PrecisionSimulation(make(), 'phase_carry')
    original.step(); baseline.step(); carried.step()
    exact_state(original, baseline); exact_state(original, carried)
    for _ in range(15):
        original.step(); baseline.step(); exact_state(original, baseline)
    both = PrecisionSimulation(make(), 'both')
    for _ in range(8): both.step()
    path = root/'gate-checkpoint.npz'; both.checkpoint(path)
    restart = PrecisionSimulation.restore(path)
    for _ in range(4):
        both.step(); restart.step(); exact_state(both, restart)
        if not torch.equal(both.phase_carry, restart.phase_carry): raise ValueError('Carry lost on restart')
    quality = both.audit(); del original, baseline, carried, both, restart
    torch.cuda.empty_cache()
    return dict(passed=True, baseline_bit_exact_steps=16, zero_carry_first_step_bit_exact=True,
        carry_restart_bit_exact_steps=4, quality=quality, checkpoint_sha256=digest(path))


def prepare(root, parent=Path('outputs/polarity-robustness-refined'), chemistry=Path('outputs/chemistry-accuracy')):
    root, parent, chemistry = [Path(x).resolve() for x in (root, parent, chemistry)]
    if root.exists(): raise FileExistsError(root)
    old = read(parent/'protocol.json'); verify(old)
    cp = read(chemistry/'protocol.json'); verify(cp); summary = read(chemistry/'summary.json')
    if (read(chemistry/'status.json')['state'] != 'completed' or not summary['passed'] or
            summary['protocol_sha256'] != digest(chemistry/'protocol.json')):
        raise ValueError('Complete and pass the frozen-chemistry diagnostic first')
    inputs = [parent/'protocol.json', chemistry/'protocol.json', chemistry/'summary.json', chemistry/'status.json',
        chemistry/'references/results.json']
    for r in summary['results']:
        folder = chemistry/r['job']['key']; path = folder/'paths.npz'
        if r['protocol_sha256'] != digest(chemistry/'protocol.json') or digest(path) != r['paths_sha256']:
            raise ValueError('Changed frozen-chemistry result')
        inputs += [path, folder/'result.json']
    for file, expected in read(chemistry/'references/results.json')['paths_sha256'].items():
        if digest(file) != expected: raise ValueError('Changed independent reference')
        inputs.append(Path(file))
    accepted, evidence = require_validation(old['validation_root']); inputs += evidence
    if asdict(accepted) != old['validated_fine_config']: raise ValueError('Changed accepted fine configuration')
    for folder in (Path(old['reference_folder']), parent/old['jobs'][0]['key']):
        file = folder/'prefix/comparison.json'; prefix = read(file)
        if not prefix['passed']: raise ValueError('Original exact-context native/GPU gate failed')
        inputs.append(file)
        for name, expected in prefix['evidence_sha256'].items():
            if digest(name) != expected: raise ValueError('Changed original context gate')
            inputs.append(Path(name))
    root.mkdir(parents=True); jobs = []; configs = {}; gates = []
    reference = old['reference_job']; chemical = Path(reference['chemical_file'])
    inputs += [Path(reference['source']), chemical]
    for level, dt in zip(('fine', 'finer'), DTS):
        folder = root/level; folder.mkdir(); source = folder/'source.npz'
        raw = retime(reference['source'], source, dt)
        if level == 'fine':
            assert_same_physical_start(AttributeSimulation.restore(reference['source']), raw)
        else: check_halving(AttributeSimulation.restore(reference['source']), raw)
        configs[level] = asdict(raw.config)
        gates.append(dict(level=level, **implementation_gate(folder, source, chemical, reference['point'], configs[level])))
        inputs += [source, folder/'gate-checkpoint.npz']
        for arm in ARMS:
            jobs.append(dict(key=f'{arm}_{level}', arm=arm, level=level, dt=dt, seed=9,
                source=str(source), chemical_file=str(chemical), family='uniform', point=reference['point']))
    binary = Path(library()._name); original_binary = Path(accepted_library()._name)
    device = dict(name=torch.cuda.get_device_name(0), torch=torch.__version__, torch_cuda=torch.version.cuda,
        experimental_binary=str(binary), experimental_binary_sha256=digest(binary),
        accepted_binary=str(original_binary), accepted_binary_sha256=digest(original_binary))
    write_json(root/'implementation-gate.json', dict(passed=all(g['passed'] for g in gates), gates=gates, device=device))
    inputs += [root/'implementation-gate.json', binary, original_binary]
    sources = [Path(__file__), *[Path(__file__).with_name(f) for f in ('gpu_precision_control.py',
        'gpu_precision_kernels.cu', 'gpu_backend.py', 'gpu_kernels.cu', 'cuda_spatial_bench.cu',
        'cuda_mechanics.cu', 'attribute_development.py', 'model.py', 'native_mechanics.py',
        'native_mechanics.cpp', 'fast_polarity.c', 'polarity.py', 'transport.py', 'signaling.py',
        'parameter_robustness_moving.py', 'parameter_robustness_reference.py', 'parameter_robustness.py',
        'feedback_long.py', 'feedback_survival_validation.py', 'polarity_robustness.py',
        'polarity_robustness_refinement.py', 'gpu_response_runner.py', 'resolution.py', 'neighbor_context.py')]]
    p = dict(jobs=jobs, accepted_configs=configs, start=150., duration=6., interval=.15,
        checkpoint_interval=3., dts=list(DTS), independent_histories=1, existing_history=9, new_histories=0,
        criteria=dict(boundary_max=.01, dilution_error_max=2e-14),
        source_sha256={str(f.resolve()):digest(f) for f in sources},
        input_sha256={str(f.resolve()):digest(f) for f in inputs}, device=device,
        design='Same mature history-9 t=150 geometry, polarity, lineage/RNG and near-uniform chemical start. D_a=.02, D_b=.55, beta=2, directional tension chi=0; activity tension/adhesion .25/.35 retained. Four numerical arms crossed with two timesteps, all for six time units. Sequential GPU jobs; no new developmental histories.',
        interventions='baseline reproduces accepted arithmetic. phase_carry retains unrounded phase increments in float64 while visible fields stay float32. contact64 accumulates contacts in float64 from the same float32 shell and therefore changes both conservative chemistry and polarity-alignment weights. both combines these interventions. Other geometry/cue/reduction arithmetic and equations stay unchanged.',
        interpretation='Diagnostic early moving window, not the nonlinear formation interval near elapsed 118. Report quality separately from numerical sensitivity. Reduced early timestep discrepancy cannot establish the cause of the long failure, backend scientific acceptance, spatial convergence, or identity differentiation. No threshold is relaxed and the original long failures remain unresolved.')
    write_json(root/'protocol.json', p); write_json(root/'status.json', dict(state='prepared', completed=0, total=8))
    return p


def checked_history(history, p, job, horizon):
    times = np.arange(_steps(horizon, p['interval'])+1)*p['interval']
    if len(history) != len(times) or not np.allclose([r['elapsed'] for r in history], times, rtol=0, atol=1e-9):
        raise ValueError('Incomplete or misaligned precision observation clocks')
    with np.load(job['chemical_file']) as z:
        if not np.array_equal(history[0]['chemistry'], z['uniform']): raise ValueError('Changed chemical start')
        ids = z['ids']
    for r in history:
        if not np.array_equal(ids, r['ids']) or abs(r['time']-p['start']-r['elapsed']) > 1e-9:
            raise ValueError('Changed IDs or physical clock')
        if np.any(np.array(r['chemistry']) <= 0) or not np.isfinite(r['chemistry']).all():
            raise ValueError('Invalid chemistry')
        validate_graph(r['delta'], r['volumes'])
    return history


def worker(root, job, p):
    root = Path(root); folder = root/job['key']; folder.mkdir(exist_ok=True); ph = digest(root/'protocol.json')
    checkpoint = folder/'latest_state.npz'
    if (folder/'result.json').exists():
        r = read(folder/'result.json')
        if (r['protocol_sha256'] != ph or r['job'] != job or r['history_sha256'] != digest(folder/'history.json')
                or r['checkpoint_sha256'] != digest(checkpoint) or not r['quality_pass']):
            raise ValueError('Changed completed precision evidence')
        checked_history(read(folder/'history.json'), p, job, p['duration']); return r
    if checkpoint.exists():
        host, audit, history = restore_checkpoint(checkpoint, ph, job)
        sim = PrecisionSimulation.restore(checkpoint)
        if sim.precision_arm != job['arm']: raise ValueError('Changed checkpoint precision arm')
        elapsed = sim.time-p['start']; checked_history(history, p, job, elapsed)
        if not np.array_equal(history[-1]['chemistry'], [host.activator, host.inhibitor]):
            raise ValueError('Checkpoint chemistry and history disagree')
    else:
        sim = PrecisionSimulation(initialize(job['source'], job['chemical_file'], 'uniform', job['point'],
            p['accepted_configs'][job['level']]), job['arm'])
        history = []; audit = dict(max_volume_error=0., min_radius=1e100, max_clipping=0.,
            boundary_max=0., dilution_error_max=0., wall_seconds=0.)
    start = _steps(p['start'], job['dt']); stop = _steps(p['start']+p['duration'], job['dt'])
    every = _steps(p['interval'], job['dt']); save_every = _steps(p['checkpoint_interval'], job['dt'])
    if sim.step_number < start or sim.step_number > stop or abs(sim.time-sim.step_number*job['dt']) > 1e-9:
        raise ValueError('Invalid checkpoint clock')
    clock = time.perf_counter(); base_wall = audit['wall_seconds']
    while sim.step_number <= stop:
        quality = sim.audit()
        for key in ('max_volume_error', 'max_clipping'): audit[key] = max(audit[key], quality[key])
        audit['min_radius'] = min(audit['min_radius'], quality['min_radius'])
        audit['dilution_error_max'] = max(audit['dilution_error_max'], quality['dilution_amount_error'])
        if audit['dilution_error_max'] > p['criteria']['dilution_error_max']: raise RuntimeError('Dilution failure')
        offset = sim.step_number-start
        if offset % every == 0:
            elapsed = offset*job['dt']
            if not history or abs(history[-1]['elapsed']-elapsed) > 1e-10:
                row = retain(sim.observe(elapsed), job); row['precision'] = sim.precision_diagnostics(); history.append(row)
            audit['boundary_max'] = max(audit['boundary_max'], history[-1]['boundary_occupancy'])
            if audit['boundary_max'] >= p['criteria']['boundary_max']: raise RuntimeError('Boundary failure')
            audit['wall_seconds'] = base_wall+time.perf_counter()-clock
            write_json(folder/'history.json', history)
            write_json(folder/'status.json', dict(state='running', elapsed=elapsed, until=p['duration'], audit=audit))
            if offset % save_every == 0 or sim.step_number == stop:
                save_checkpoint(sim, checkpoint, audit, history, ph, job)
                print(f'{job["key"]}: elapsed={elapsed:g}, wall={audit["wall_seconds"]:.1f}s', flush=True)
        if sim.step_number == stop: break
        sim.step()
    checked_history(history, p, job, p['duration'])
    result = dict(job=job, protocol_sha256=ph, history_sha256=digest(folder/'history.json'),
        checkpoint_sha256=digest(checkpoint), quality_pass=True, audit=audit,
        final_log_sd=history[-1]['log_activator_sd'], precision=history[-1]['precision'])
    write_json(folder/'result.json', result); write_json(folder/'status.json', dict(state='completed', elapsed=p['duration']))
    del sim; torch.cuda.empty_cache(); return result


def errors(a, b):
    if not np.array_equal([r['elapsed'] for r in a], [r['elapsed'] for r in b]): raise ValueError('Unaligned comparison')
    result = {}
    for key, name, transform in (('chemistry', 'chemical_log_max', np.log),
            ('volumes', 'relative_volume_max', lambda x:x), ('delta', 'relative_transport_max', lambda x:x),
            ('polarity', 'polarity_abs_max', lambda x:x), ('centers', 'center_abs_max', lambda x:x)):
        x, y = [np.array([r[key] for r in rows]) for rows in (a, b)]
        difference = abs(transform(x)-transform(y))
        result[name] = float(np.max(difference/np.maximum(abs(y), 1e-30))) if key == 'volumes' else float(
            np.max(np.linalg.norm(x-y, axis=(1, 2))/np.maximum(np.linalg.norm(y, axis=(1, 2)), 1e-30))) if key == 'delta' else float(difference.max())
    result['growth_abs_max'] = float(max(abs(x['uniform_growth_max']-y['uniform_growth_max']) for x, y in zip(a, b)))
    return result


def assess(root):
    root = Path(root); p = read(root/'protocol.json'); verify(p); rows = []; histories = {}
    for job in p['jobs']:
        folder = root/job['key']; r = read(folder/'result.json')
        if (r['job'] != job or r['protocol_sha256'] != digest(root/'protocol.json') or not r['quality_pass']
                or r['history_sha256'] != digest(folder/'history.json') or r['checkpoint_sha256'] != digest(folder/'latest_state.npz')):
            raise ValueError('Changed precision result')
        host, audit, h = restore_checkpoint(folder/'latest_state.npz', digest(root/'protocol.json'), job)
        checked_history(h, p, job, p['duration'])
        if (h != read(folder/'history.json') or audit != r['audit'] or host.time != p['start']+p['duration']
                or not np.array_equal(h[-1]['chemistry'], [host.activator, host.inhibitor])):
            raise ValueError('Checkpoint and history disagree')
        with np.load(folder/'latest_state.npz') as z:
            if str(z['precision_arm']) != job['arm'] or z['phase_carry'].dtype != np.float64 or not np.isfinite(z['phase_carry']).all():
                raise ValueError('Invalid saved precision state')
        histories[job['arm'], job['level']] = h; rows.append(r)
    comparisons = [dict(arm=arm, **errors(histories[arm, 'fine'], histories[arm, 'finer'])) for arm in ARMS]
    interventions = [dict(arm=arm, level=level, **errors(histories[arm, level], histories['baseline', level]))
        for arm in ARMS[1:] for level in ('fine', 'finer')]
    # Endpoint fields complement graph diagnostics; not a spatial-convergence test.
    for comparison in comparisons:
        arm = comparison['arm']
        with np.load(root/f'{arm}_fine/latest_state.npz') as a, np.load(root/f'{arm}_finer/latest_state.npz') as b:
            x, y = a['phi'].astype(np.float64), b['phi'].astype(np.float64)
            comparison['endpoint_phase_abs_max'] = float(abs(x-y).max())
            comparison['endpoint_phase_relative_l2'] = float(np.linalg.norm(x-y)/np.linalg.norm(y))
    summary = dict(protocol_sha256=digest(root/'protocol.json'), completed=len(rows), total=8,
        quality_pass=all(r['quality_pass'] for r in rows), independent_histories=1, new_histories=0,
        timestep_comparisons=comparisons, intervention_comparisons=interventions, results=rows, scope=p['interpretation'])
    write_json(root/'summary.json', summary); return summary


def run(root):
    root = Path(root).resolve()
    with (root/'coordinator.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        p = read(root/'protocol.json'); verify(p)
        d = p['device']
        if (d['name'] != torch.cuda.get_device_name(0) or d['torch'] != torch.__version__ or
                d['torch_cuda'] != torch.version.cuda or digest(library()._name) != d['experimental_binary_sha256'] or
                digest(accepted_library()._name) != d['accepted_binary_sha256']):
            raise ValueError('Changed GPU hardware/software')
        done = 0
        try:
            for job in p['jobs']:
                write_json(root/'status.json', dict(state='running', current_job=job['key'], completed=done, total=8))
                worker(root, job, p); done += 1
            summary = assess(root)
            write_json(root/'status.json', dict(state='completed', completed=done, total=8, quality_pass=summary['quality_pass']))
            return summary
        except Exception as error:
            write_json(root/'status.json', dict(state='failed', completed=done, total=8, error=str(error))); raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'run', 'assess'))
    parser.add_argument('--output', type=Path, default=Path('outputs/geometry-precision')); args = parser.parse_args()
    report = prepare(args.output) if args.action == 'prepare' else run(args.output) if args.action == 'run' else assess(args.output)
    print({k:v for k,v in report.items() if k in ('completed', 'total', 'quality_pass', 'timestep_comparisons')})
