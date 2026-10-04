"""Sparse two-parameter GPU follow-up to the frozen robustness map.

New physical parameters receive their own native/GPU prefix gate. The existing
full-horizon baseline gate is retained and is not extended by assertion.
"""
import argparse
from dataclasses import asdict
import fcntl
from pathlib import Path
import time

import numpy as np

from .attribute_development import AttributeSimulation
from .attribute_persistence import perturb
from .benchmark_gpu_backend import cpu_audit
from .cell_response_moving import observe
from .feedback_long import digest, save_checkpoint, restore_checkpoint
from .feedback_survival_validation import retime
from .gpu_response_runner import require_validation, compatible, PARAMETERS
from .neighbor_context import read, validate_graph
from .parameter_robustness import verify, frozen_spectrum, worker as frozen_worker, CRITERIA
from .resolution import write_json, _steps
from .validate_gpu_backend import CRITERIA as GPU_CRITERIA, discrepancies


def validate_point(point):
    ratio, chi = point['ratio'], point['chi']
    if not (np.isfinite(ratio) and ratio > 0 and np.isfinite(chi) and 0 <= chi < 1):
        raise ValueError('Positive ratio and 0 <= polarity-tension contrast < 1 required')


def initialize(source, chemical_file, family, point, accepted_config):
    """Change only D_b and the directional-tension coefficient, after baseline checking."""
    validate_point(point)
    host = AttributeSimulation.restore(source)
    from .model import Config
    accepted = Config(**accepted_config)
    compatible(host.config, accepted)
    if (host.divisions or host.attribute_mode != 'direct' or not host.config.feedback or
            len(host.ids) != host.config.max_cells or abs(host.time-150.) > 1e-9):
        raise ValueError('Requires a direct-feedback mature t=150 starting geometry')
    with np.load(chemical_file) as z:
        if family not in ('pattern', 'uniform') or not np.array_equal(z['ids'], host.ids):
            raise ValueError('Chemical family or cell order mismatch')
        state = z[family].copy()
    if state.shape != (2, len(host.ids)) or not np.isfinite(state).all() or np.any(state <= 0):
        raise ValueError('Invalid positive chemical start')
    host.activator, host.inhibitor = state
    host.config.signal_dh = host.config.signal_da*point['ratio']
    host.config.polarity_tension = point['chi']
    host.config.validate()
    unexpected = [k for k in PARAMETERS if k not in ('signal_dh', 'polarity_tension') and
                  getattr(host.config, k) != getattr(accepted, k)]
    if unexpected:
        raise ValueError('Unexpected parameter changes: '+', '.join(unexpected))
    return host


def prepare(root, frozen=Path('outputs/parameter-robustness'), validation=Path('outputs/gpu-backend-validation')):
    from .gpu_backend import library
    import torch
    root, frozen, validation = [Path(x).resolve() for x in (root, frozen, validation)]
    if root.exists():
        raise FileExistsError(root)
    fp = read(frozen/'protocol.json'); verify(fp)
    status = read(frozen/'status.json'); summary = read(frozen/'summary.json')
    if status['state'] != 'completed' or not summary['numerical_pass']:
        raise ValueError('Complete the numerically accepted frozen map first')
    selection = read(frozen/'moving-selection.json')
    if (selection['frozen_summary_sha256'] != digest(frozen/'summary.json') or
            selection['protocol_sha256'] != digest(frozen/'protocol.json')):
        raise ValueError('Moving selection no longer matches frozen evidence')
    accepted, evidence = require_validation(validation)
    points = selection['points']
    for point in points:
        validate_point(point)
    inputs = {str(f):digest(f) for f in evidence}
    for f in ('protocol.json', 'status.json', 'summary.json', 'moving-selection.json'):
        inputs[str(frozen/f)] = digest(frozen/f)
    inputs.update(fp['input_sha256'])
    root.mkdir(parents=True)
    jobs = []
    for g in fp['graphs']:
        if g['branch'] != 'switch_on':
            continue
        folder = root/f"seed-{g['seed']}"; folder.mkdir()
        source, chemical = folder/'source.npz', folder/'initial_states.npz'
        host = retime(g['checkpoint'], source, accepted.dt)
        compatible(host.config, accepted)
        with np.load(g['source']) as z:
            pattern = z['trajectories'][0, -1].copy()
            uniform = perturb(np.ones_like(pattern), z['volumes'], fp['uniform_log_noise'],
                              np.random.default_rng(g['seed']))
            np.savez_compressed(chemical, pattern=pattern, uniform=uniform, ids=z['ids'], masses=z['volumes'])
        for point in points:
            for family in ('uniform', 'pattern'):
                jobs.append(dict(key=f"seed-{g['seed']}_{point['key']}_{family}", seed=g['seed'],
                    point=point, family=family, source=str(source), chemical_file=str(chemical)))
        for f in (source, chemical):
            inputs[str(f)] = digest(f)
    source_names = ('gpu_backend.py', 'gpu_kernels.cu', 'cuda_mechanics.cu', 'cuda_spatial_bench.cu',
        'native_mechanics.py', 'native_mechanics.cpp', 'fast_mechanics.py', 'fast_mechanics.c',
        'attribute_development.py', 'attribute_persistence.py', 'model.py', 'polarity.py', 'signaling.py',
        'transport.py', 'benchmark_gpu_backend.py', 'validate_gpu_backend.py', 'gpu_response_runner.py',
        'cell_response_moving.py', 'feedback_long.py', 'feedback_survival_validation.py',
        'neighbor_context.py', 'parameter_robustness.py', 'cell_response_exchange.py',
        'feedback_endpoint_bistability.py', 'resolution.py')
    files = [Path(__file__), *[Path(__file__).with_name(f) for f in source_names]]
    p = dict(points=points, selection_criterion_satisfied=selection['criterion_satisfied'],
        frozen_study=str(frozen), gpu_validation=str(validation), accepted_config=asdict(accepted),
        jobs=jobs, total_jobs=len(jobs), independent_histories=3, start=150., duration=60.,
        dt=accepted.dt, interval=.15, checkpoint_interval=3., late_window=12.,
        prefix_duration=.6, prefix_native_threads=4, prefix_criteria=GPU_CRITERIA,
        criteria=dict(boundary_max=.01, dilution_error_max=GPU_CRITERIA['dilution_amount_error'],
                      late_log_sd_min=.1), endpoint_criteria=CRITERIA,
        endpoint_horizons=fp['horizons'], endpoint_interval=fp['interval'],
        design='At each of three parameter points and three existing histories, start matched near-uniform and previously developed baseline-ratio chemistry on identical t=150 switch-on geometry/polarity. All tension/adhesion response parameters stay fixed. Only D_b and polarity-tension contrast change. Chemistry and mechanics co-evolve to t=210; analyze instantaneous spectra and frozen endpoints separately.',
        gate='Retain full-horizon four-run baseline GPU acceptance and verify hardware/software/binary hashes. Before EVERY new job, require .6 units of native/GPU agreement at its exact new parameters, geometry and chemical start using existing strict field/chemistry/transport thresholds. No bypass of the original compatible() gate; this module supplies a separately scoped experimental acceptance layer.',
        limits='New-parameter GPU acceptance is short-horizon only, not full-horizon backend equivalence or timestep/spatial convergence. Three histories, eighteen nested moving interventions. Mature near-uniform restart tests formation opportunity, not initiation along fresh zygote trajectories. Frozen endpoints assess chemical local attractors, not full moving-system stability. No interpolation over untested points.',
        endpoint_design='At each moving endpoint freeze its OWN measured graph. Start from its current chemistry plus two 1% log perturbations and two 0.1% near-uniform log perturbations; use the same dual-solver stationarity and local-return tests as the original frozen map. Endpoint chemistry may be uniform: no patterned basin is then inferred from its loss.',
        device=dict(name=torch.cuda.get_device_name(0), torch=torch.__version__, torch_cuda=torch.version.cuda,
                    binary_sha256=digest(library()._name)),
        source_sha256={str(f.resolve()):digest(f) for f in files}, input_sha256=inputs)
    write_json(root/'protocol.json', p)
    write_json(root/'status.json', dict(state='prepared', completed=0, total=len(jobs)))
    return p


def prefix(root, job, p):
    from .gpu_backend import GpuSimulation
    from .native_mechanics import NativeSimulation
    folder = Path(root)/job['key']/'prefix'; folder.mkdir(parents=True, exist_ok=True)
    ph = digest(Path(root)/'protocol.json')
    if (folder/'comparison.json').exists():
        result = read(folder/'comparison.json')
        if result['protocol_sha256'] != ph or not result['passed']:
            raise ValueError('Invalid saved new-parameter gate')
        for f, h in result['evidence_sha256'].items():
            if digest(f) != h:
                raise ValueError('Changed prefix evidence')
        return result
    cpu = initialize(job['source'], job['chemical_file'], job['family'], job['point'], p['accepted_config'])
    cpu.__class__ = NativeSimulation; cpu.native_threads = p['prefix_native_threads']
    gpu = GpuSimulation(initialize(job['source'], job['chemical_file'], job['family'], job['point'], p['accepted_config']))
    steps = _steps(p['prefix_duration'], p['dt']); every = _steps(p['interval'], p['dt'])
    histories = [[], []]
    errors = dict.fromkeys(('chemical_log_max', 'polarity_abs_max', 'relative_axis_max',
                           'relative_volume_max', 'relative_transport_max'), 0.)
    started = time.perf_counter()
    for step in range(steps+1):
        cpu_audit(cpu); audit = gpu.audit()
        if audit['dilution_amount_error'] > p['criteria']['dilution_error_max']:
            raise RuntimeError('Prefix dilution conservation failure')
        if step % every == 0:
            a, b = gpu.observe(step*p['dt']), observe(cpu, step*p['dt'])
            if max(a['boundary_occupancy'], b['boundary_occupancy']) >= p['criteria']['boundary_max']:
                raise RuntimeError('New-parameter prefix boundary failure')
            for key, value in discrepancies(a, b).items():
                errors[key] = max(errors[key], value)
            histories[0].append(a); histories[1].append(b)
            write_json(folder/'status.json', dict(state='running', elapsed=step*p['dt'], errors=errors))
        if step < steps:
            cpu.step(); gpu.step()
    phi_error = float(abs(gpu.phi.cpu().numpy()-cpu.phi).max())
    passed = all(v <= p['prefix_criteria'][k] for k, v in errors.items()) and phi_error <= p['prefix_criteria']['final_phi_abs_max']
    files = [folder/'gpu-history.json', folder/'cpu-history.json', folder/'gpu-end.npz', folder/'cpu-end.npz']
    write_json(files[0], histories[0]); write_json(files[1], histories[1]); gpu.checkpoint(files[2]); cpu.checkpoint(files[3])
    result = dict(passed=bool(passed), protocol_sha256=ph, job=job, errors=errors, phi_abs_max=phi_error,
        wall_seconds=time.perf_counter()-started, evidence_sha256={str(f.resolve()):digest(f) for f in files},
        scope='Same-state, new-parameter .6-unit native/GPU check; not full-horizon acceptance.')
    write_json(folder/'comparison.json', result)
    write_json(folder/'status.json', dict(state='completed', passed=bool(passed)))
    if not passed:
        raise RuntimeError('New-parameter native/GPU disagreement; scientific continuation blocked')
    return result


def retain(row, job):
    spectrum = frozen_spectrum(np.array(row['delta']), np.array(row['volumes']), 2., .02, job['point']['ratio'])
    return {**row, 'uniform_growth_max':spectrum['uniform_jacobian_max_real'],
            'unstable_modes':spectrum['unstable_modes'], 'laplacian_lambdas':spectrum['lambdas']}


def endpoint_assay(root, job, p, final):
    """Reuse the hashed frozen chemical worker on each actual moving endpoint."""
    folder = Path(root)/job['key']/'endpoint'; folder.mkdir(exist_ok=True)
    if not (folder/'protocol.json').exists():
        state = np.array(final['chemistry']); delta = np.array(final['delta']); masses = np.array(final['volumes'])
        validate_graph(delta, masses)
        source = folder/'source.npz'
        np.savez_compressed(source, trajectories=state[None, None], delta=delta, volumes=masses, ids=final['ids'])
        fp = dict(noise_seeds=[0, 1], pattern_log_noise=.01, uniform_log_noise=.001,
            horizons=p['endpoint_horizons'], interval=p['endpoint_interval'], criteria=p['endpoint_criteria'],
            source_sha256=p['source_sha256'], input_sha256={str(source.resolve()):digest(source)},
            moving_protocol_sha256=digest(Path(root)/'protocol.json'),
            moving_history_sha256=digest(Path(root)/job['key']/'history.json'))
        write_json(folder/'protocol.json', fp)
    fp = read(folder/'protocol.json'); verify(fp)
    if (fp['moving_protocol_sha256'] != digest(Path(root)/'protocol.json') or
            fp['moving_history_sha256'] != digest(Path(root)/job['key']/'history.json')):
        raise ValueError('Endpoint assay source changed')
    frozen_job = dict(job_key='assay', seed=job['seed'], branch=job['family'], ratio=job['point']['ratio'],
                      db=.02*job['point']['ratio'], beta=2., da=.02, source=str(folder/'source.npz'))
    return frozen_worker((folder, frozen_job, fp))


def continuation(root, job, p):
    from .gpu_backend import GpuSimulation
    folder = Path(root)/job['key']; folder.mkdir(exist_ok=True)
    ph = digest(Path(root)/'protocol.json')
    checkpoint = folder/'latest_state.npz'
    prefix(root, job, p)
    if (folder/'result.json').exists():
        r = read(folder/'result.json')
        if (r['protocol_sha256'] != ph or digest(checkpoint) != r['checkpoint_sha256'] or
                digest(folder/'history.json') != r['history_sha256']):
            raise ValueError('Changed completed moving evidence')
        endpoint_assay(root, job, p, read(folder/'history.json')[-1])
        return r
    if checkpoint.exists():
        host, audit, history = restore_checkpoint(checkpoint, ph, job)
        sim = GpuSimulation(host)
    else:
        sim = GpuSimulation(initialize(job['source'], job['chemical_file'], job['family'], job['point'], p['accepted_config']))
        history = []
        audit = dict(max_volume_error=0., min_radius=1e100, max_clipping=0., boundary_max=0.,
                     dilution_error_max=0., wall_seconds=0.)
    start_step = _steps(p['start'], p['dt']); stop = _steps(p['start']+p['duration'], p['dt'])
    every = _steps(p['interval'], p['dt']); checkpoint_every = _steps(p['checkpoint_interval'], p['dt'])
    if not start_step <= sim.step_number <= stop:
        raise ValueError('Invalid resumed simulation clock')
    started = time.perf_counter(); base_wall = audit['wall_seconds']
    try:
        while sim.step_number <= stop:
            current = sim.audit()
            for key in ('max_volume_error', 'max_clipping'):
                audit[key] = max(audit[key], current[key])
            audit['min_radius'] = min(audit['min_radius'], current['min_radius'])
            audit['dilution_error_max'] = max(audit['dilution_error_max'], current['dilution_amount_error'])
            if audit['dilution_error_max'] > p['criteria']['dilution_error_max']:
                raise RuntimeError('Moving dilution conservation failure')
            offset = sim.step_number-start_step
            if offset % every == 0:
                elapsed = offset*p['dt']
                if not history or abs(history[-1]['elapsed']-elapsed) > 1e-10:
                    row = retain(sim.observe(elapsed), job); history.append(row)
                else:
                    row = history[-1]
                audit['boundary_max'] = max(audit['boundary_max'], row['boundary_occupancy'])
                if audit['boundary_max'] >= p['criteria']['boundary_max']:
                    raise RuntimeError('Moving boundary screen failed')
                audit['wall_seconds'] = base_wall+time.perf_counter()-started
                write_json(folder/'status.json', dict(state='running', elapsed=elapsed, audit=audit,
                    log_activator_sd=row['log_activator_sd'], uniform_growth_max=row['uniform_growth_max']))
                if offset % checkpoint_every == 0 or sim.step_number == stop:
                    save_checkpoint(sim.to_cpu(), checkpoint, audit, history, ph, job)
                    write_json(folder/'history.json', history)
                    print(f"{job['key']} t={row['time']:.2f} SD={row['log_activator_sd']:.4f} growth={row['uniform_growth_max']:.4f}", flush=True)
            if sim.step_number == stop:
                break
            sim.step()
        late = [r for r in history if r['elapsed'] >= p['duration']-p['late_window']-1e-9]
        spread = min(r['log_activator_sd'] for r in late)
        result = dict(job=job, quality_pass=True, protocol_sha256=ph, audit=audit,
            history_sha256=digest(folder/'history.json'), checkpoint_sha256=digest(checkpoint),
            late_min_log_activator_sd=spread,
            persistent_contrast=bool(spread > p['criteria']['late_log_sd_min']),
            interpretation='Late chemical contrast in a moving continuation; not chemical stability or irreversible identity.')
        write_json(folder/'result.json', result)
        endpoint = endpoint_assay(root, job, p, history[-1])
        write_json(folder/'status.json', dict(state='completed', elapsed=p['duration'], quality_pass=True,
                                            endpoint_numerical_pass=endpoint['numerical_pass']))
        return result
    except Exception as error:
        write_json(folder/'status.json', dict(state='failed', time=sim.time, error=str(error), audit=audit))
        raise


def assess(root):
    root = Path(root); p = read(root/'protocol.json'); verify(p)
    rows, failures = [], []
    for job in p['jobs']:
        folder = root/job['key']
        sf = folder/'status.json'
        if sf.exists() and read(sf)['state'] == 'failed':
            failures.append(dict(job=job, error=read(sf)['error']))
        if not (folder/'result.json').exists():
            continue
        r = read(folder/'result.json')
        if (r['protocol_sha256'] != digest(root/'protocol.json') or
                digest(folder/'history.json') != r['history_sha256'] or
                digest(folder/'latest_state.npz') != r['checkpoint_sha256']):
            raise ValueError('Changed moving evidence')
        endpoint_file = folder/'endpoint'/'assay'/'result.json'
        endpoint = read(endpoint_file) if endpoint_file.exists() else None
        if endpoint is not None:
            ep = read(folder/'endpoint'/'protocol.json'); verify(ep)
            if (endpoint['protocol_sha256'] != digest(folder/'endpoint'/'protocol.json') or
                    endpoint['paths_sha256'] != digest(endpoint_file.with_name('paths.npz'))):
                raise ValueError('Changed endpoint assay evidence')
        rows.append(dict(job=job, persistent_contrast=r['persistent_contrast'],
            late_min_log_activator_sd=r['late_min_log_activator_sd'], quality_pass=r['quality_pass'], endpoint=endpoint))
    endpoints = [r['endpoint'] for r in rows if r['endpoint'] is not None]
    summary = dict(independent_histories=3, completed=len(rows), endpoint_completed=len(endpoints),
                   endpoint_numerical_pass=bool(len(endpoints) == len(p['jobs']) and all(e['numerical_pass'] for e in endpoints)),
                   endpoint_unresolved=sum(e['phase'] == 'unresolved' for e in endpoints),
                   total=len(p['jobs']), results=rows, failures=failures, protocol_sha256=digest(root/'protocol.json'))
    write_json(root/'summary.json', summary)
    plot(root, p, rows)
    return summary


def plot(root, p, rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharex=True, sharey=True)
    for ax, family in zip(axes, ('uniform', 'pattern')):
        for point in p['points']:
            completed = [r for r in rows if r['job']['point'] == point and r['job']['family'] == family]
            counts = sum(r['persistent_contrast'] for r in completed)
            color = '#d2d5da' if not completed else '#3879bd' if counts == len(completed) else '#db8a25' if counts else '#52565d'
            ax.scatter(point['ratio'], point['chi'], c=color, s=170, edgecolors='black', zorder=3)
            ax.annotate(f"{point['key'].replace('_',' ')}\ncontrast {counts}/{len(completed)} histories", (point['ratio'], point['chi']),
                        xytext=(8, 5), textcoords='offset points', fontsize=8)
        ax.set(title='Near-uniform start' if family == 'uniform' else 'Developed-pattern start',
               xlabel='D_b / D_a (D_a = 0.02)', ylabel='Polarity-tension contrast chi',
               xlim=(min(x['ratio'] for x in p['points'])-3, max(x['ratio'] for x in p['points'])+9), ylim=(-.1, .85))
        ax.grid(alpha=.2)
    fig.suptitle('Moving parameter spot checks: sampled points only; grey = pending', fontsize=11)
    fig.tight_layout(); fig.savefig(root/'moving-parameter-points.png', dpi=180); plt.close(fig)


def run(root):
    import torch
    root = Path(root); p = read(root/'protocol.json'); verify(p)
    accepted, _ = require_validation(p['gpu_validation'])
    if asdict(accepted) != p['accepted_config']:
        raise ValueError('Baseline GPU gate changed')
    torch.set_num_threads(1)
    with (root/'coordinator.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        for job in p['jobs']:
            completed = sum((root/j['key']/'result.json').exists() for j in p['jobs'])
            write_json(root/'status.json', dict(state='running', completed=completed, total=len(p['jobs']), current_job=job['key']))
            try:
                continuation(root, job, p)
            except (RuntimeError, FloatingPointError) as error:
                folder = root/job['key']; folder.mkdir(exist_ok=True)
                write_json(folder/'failure.json', dict(job=job, error=str(error), protocol_sha256=digest(root/'protocol.json')))
                if not (folder/'status.json').exists() or read(folder/'status.json')['state'] != 'failed':
                    write_json(folder/'status.json', dict(state='failed', error=str(error)))
                print(f"FAILED {job['key']}: {error}", flush=True)
            assess(root)
        summary = assess(root)
        accepted_all = summary['completed'] == p['total_jobs'] and summary['endpoint_completed'] == p['total_jobs'] and not summary['failures']
        write_json(root/'status.json', dict(state='completed' if accepted_all else 'completed_with_failures',
            completed=summary['completed'], total=len(p['jobs']), endpoint_completed=summary['endpoint_completed']))
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'run', 'assess'))
    parser.add_argument('--output', type=Path, default=Path('outputs/parameter-robustness-moving'))
    parser.add_argument('--frozen', type=Path, default=Path('outputs/parameter-robustness'))
    parser.add_argument('--validation', type=Path, default=Path('outputs/gpu-backend-validation'))
    args = parser.parse_args()
    result = prepare(args.output, args.frozen, args.validation) if args.action == 'prepare' else run(args.output) if args.action == 'run' else assess(args.output)
    print({k:v for k,v in result.items() if k not in ('results','jobs','accepted_config','source_sha256','input_sha256')})
