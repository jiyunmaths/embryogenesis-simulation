"""Geometry-conditional initiation/maintenance map on accepted endpoint graphs.

Polarity-tension contrast is NOT an argument of a frozen chemical operator.
The second parameter axis is tested separately by parameter-gated moving runs.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import fcntl
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq

from .attribute_development import AttributeSimulation
from .attribute_persistence import rhs, perturb, distance
from .cell_response_exchange import checked_solve
from .feedback_endpoint_bistability import chemical_jacobian
from .feedback_long import digest
from .neighbor_context import read, validate_graph
from .resolution import write_json


RATIOS = (8., 10., 12., 15., 17.5, 20., 22.5, 25., 27.5, 30., 40.)
CRITERIA = dict(solver_log_max=1e-5, stationarity_rhs_max=1e-6,
                stable_growth_margin=1e-8, pattern_log_sd_min=.1,
                uniform_log_max=1e-4, pattern_return_log_rms_max=1e-4)
SEEDS = (0, 1)


def verify(p):
    for field in ('source_sha256', 'input_sha256'):
        for filename, expected in p[field].items():
            if digest(filename) != expected:
                raise ValueError('Changed study dependency: '+filename)


def validate_parameters(beta, da, ratio):
    if not all(np.isfinite(v) and v > 0 for v in (beta, da, ratio)):
        raise ValueError('Positive finite beta, D_a and D_b/D_a required')
    if beta <= 1:
        raise ValueError('This initiation assay requires stable isolated-cell uniform chemistry (beta > 1)')
    return da*ratio


def frozen_spectrum(delta, masses, beta, da, ratio):
    """Exact finite-graph dispersion; positive lambda are eigenvalues of -Delta_V."""
    delta, masses = validate_graph(delta, masses)
    db = validate_parameters(beta, da, ratio)
    symmetric = -np.sqrt(masses)[:, None]*delta/np.sqrt(masses)[None, :]
    lambdas = np.linalg.eigvalsh(symmetric)
    if lambdas[0] < -1e-10:
        raise ValueError('Negative conservative Laplacian eigenvalue')
    lambdas = np.maximum(lambdas, 0.)
    growth = np.array([np.linalg.eigvals([[1-da*l, -1], [2*beta, -beta-db*l]]).real.max()
                       for l in lambdas])
    nonconstant = lambdas > max(lambdas[-1], 1.)*1e-12
    maximum = float(growth[nonconstant].max()) if nonconstant.any() else None
    # det(J_lambda) = beta + (beta D_a-D_b)lambda + D_a D_b lambda^2.
    discriminant = (beta*da-db)**2-4*da*db*beta
    band = None
    if db > beta*da and discriminant > 0:
        band = sorted(((db-beta*da-np.sqrt(discriminant))/(2*da*db),
                       (db-beta*da+np.sqrt(discriminant))/(2*da*db)))
    return dict(lambdas=lambdas.tolist(), modal_growth=growth.tolist(),
                maximum_spatial_growth=maximum,
                unstable_modes=int(np.count_nonzero(nonconstant & (growth > CRITERIA['stable_growth_margin']))),
                uniform_jacobian_max_real=float(growth.max()),
                continuous_instability_band=band)


def uniform_threshold(delta, masses, beta, da, lower=8., upper=40.):
    fun = lambda r: frozen_spectrum(delta, masses, beta, da, r)['uniform_jacobian_max_real']
    if fun(lower)*fun(upper) >= 0:
        return None
    return float(brentq(fun, lower, upper, xtol=1e-10))


def prepare(root, ratios=RATIOS):
    root = Path(root).resolve()
    if root.exists():
        raise FileExistsError(root)
    ratios = sorted(set(float(r) for r in ratios))
    if len(ratios) < 3 or 20. not in ratios:
        raise ValueError('Use at least three distinct ratios including baseline 20')
    for ratio in ratios:
        validate_parameters(2., .02, ratio)
    graphs, inputs = [], {}
    for seed in (7, 8, 9):
        source = Path('outputs/feedback-endpoint-bistability'+('' if seed == 7 else f'-seed-{seed}')).resolve()
        prior = read(source/'protocol.json')
        for f, h in prior['sha256'].items():
            if digest(f) != h:
                raise ValueError('Changed accepted endpoint dependency: '+f)
        results = read(source/'results.json')
        cp_root = Path('outputs/feedback-survival' if seed == 7 else
                       f'outputs/feedback-survival-validation/seed-{seed}/survival').resolve()
        for branch in ('switch_on', 'keep_off'):
            if not results[branch]['chemical_bistability_supported']:
                raise ValueError('Requires accepted endpoint bistability')
            file = source/(branch+'.npz')
            checkpoint = cp_root/branch/'latest_state.npz'
            host = AttributeSimulation.restore(checkpoint)
            with np.load(file) as z:
                delta, masses = validate_graph(z['delta'], z['volumes'])
                graph = host.signaling_graph()
                if (abs(host.time-150.) > 1e-9 or not np.array_equal(z['ids'], host.ids) or
                        not np.allclose(delta, graph.delta, rtol=1e-12, atol=1e-12) or
                        not np.allclose(masses, graph.masses, rtol=1e-12, atol=1e-14)):
                    raise ValueError('Endpoint geometry/checkpoint mismatch')
            c = host.config
            if (c.signal_beta, c.signal_da, c.signal_dh) != (2., .02, .4):
                raise ValueError('Unexpected accepted chemical baseline')
            graphs.append(dict(key=f'seed-{seed}_{branch}', seed=seed, branch=branch,
                source=str(file), checkpoint=str(checkpoint), beta=c.signal_beta, da=c.signal_da,
                baseline_ratio=20., geometry_time=150., polarity_tension=c.polarity_tension,
                mechanical_feedback=host.attribute_mode,
                uniform_instability_ratio=uniform_threshold(delta, masses, c.signal_beta, c.signal_da)))
            for f in (file, checkpoint):
                inputs[str(f)] = digest(f)
        for f in (source/'protocol.json', source/'results.json'):
            inputs[str(f)] = digest(f)
    sources = [Path(__file__), *[Path(__file__).with_name(f) for f in
        ('attribute_persistence.py', 'feedback_endpoint_bistability.py', 'cell_response_exchange.py',
         'neighbor_context.py', 'feedback_long.py', 'resolution.py')]]
    root.mkdir(parents=True)
    jobs = [dict(**g, ratio=r, db=g['da']*r, job_key=f"{g['key']}_ratio-{r:g}") for g in graphs for r in ratios]
    p = dict(graphs=graphs, jobs=jobs, ratios=ratios, criteria=CRITERIA, noise_seeds=list(SEEDS),
        uniform_log_noise=.001, pattern_log_noise=.01, horizons=[240., 960.], interval=2.,
        independent_histories=3, graph_count=6, total_jobs=len(jobs), trajectories_per_job=5,
        design='Fix beta=2 and D_a=.02; vary D_b=.02 ratio. On each accepted t=150 graph, integrate the baseline-ratio developed equilibrium and two 1% log perturbations; integrate two 0.1% near-uniform log perturbations. Species amounts are preserved at perturbation. Exact uniform stability is analytical. Repeat the same initial states at each ratio, without parameter continuation or replacing lost branches.',
        extension='Extend the identical initial condition from 240 to 960 only when stationarity/classification is unresolved; no favorable-outcome selection. Keep cross-solver failure and slow/oscillatory endpoints unresolved.',
        classification='Local frozen coexistence requires stable uniform chemistry AND a stationary linearly stable patterned endpoint AND both developed-pattern perturbations returning to it. Initiation additionally requires positive homogeneous modal growth and both small-noise trials settling to stable patterned endpoints. Failure of sampled pattern starts is not proof that no patterned attractor exists.',
        moving_selection='Three parameter points: baseline (20,.35), polarity-tension ablation (20,0), and (r,.7) where r is the smallest sampled ratio >20 with homogeneous growth >.02 and supported maintenance on all six graphs; if none qualifies, use largest sampled ratio and flag selection as inconclusive. Three source histories and two chemical starts per point, identical t=150 switch-on geometry within each history.',
        limits='Frozen chemistry has no polarity-tension parameter. Six geometries are nested in three histories, not six independent replicates. Full/no-feedback geometries differ in several couplings: their spectral difference does not isolate polarity. Moving spot checks are mature-state pilot robustness, not a new zygote initiation study, a spatial-convergence proof, or an exhaustive bifurcation/attractor census.',
        source_sha256={str(f.resolve()): digest(f) for f in sources}, input_sha256=inputs)
    write_json(root/'protocol.json', p)
    write_json(root/'status.json', dict(state='prepared', completed=0, total=len(jobs)))
    return p


def metrics(state, fun, jac, criteria):
    residual = float(abs(fun(0., state.ravel())).max())
    growth = float(np.linalg.eigvals(jac(0., state.ravel())).real.max())
    spread = float(np.std(np.log(state[0])))
    deviation = float(abs(np.log(state)).max())
    stable = growth < -criteria['stable_growth_margin']
    stationary = residual < criteria['stationarity_rhs_max']
    kind = 'unresolved'
    if stationary and stable and spread > criteria['pattern_log_sd_min']:
        kind = 'patterned'
    elif stationary and stable and deviation < criteria['uniform_log_max']:
        kind = 'uniform'
    return dict(rhs_max=residual, jacobian_max_real=growth, log_activator_sd=spread,
                uniform_log_max=deviation, endpoint_kind=kind)


def classify(trials, uniform_growth, criteria):
    numerical = len(trials) == 5 and all(t.get('solver_pass', False) for t in trials)
    settled = numerical and all(t['endpoint_kind'] != 'unresolved' for t in trials)
    pattern = trials[:3]
    maintenance = numerical and len(pattern) == 3 and all(
        t['endpoint_kind'] == 'patterned' and t['pattern_return_log_rms'] < criteria['pattern_return_log_rms_max']
        for t in pattern)
    uniform_stable = uniform_growth < -criteria['stable_growth_margin']
    coexistence = bool(maintenance and uniform_stable)
    initiation = bool(numerical and uniform_growth > criteria['stable_growth_margin'] and
                      all(t['endpoint_kind'] == 'patterned' for t in trials[3:]))
    phase = 'unresolved'
    if settled:
        if coexistence:
            phase = 'coexistence'
        elif initiation and maintenance:
            phase = 'initiation_and_maintenance'
        elif all(t['endpoint_kind'] == 'uniform' for t in trials):
            phase = 'uniform_from_sampled_starts'
        else:
            phase = 'mixed_sampled_basins'
    return dict(numerical_pass=bool(numerical), all_trials_settled=bool(settled),
                uniform_stable=bool(uniform_stable), maintenance_supported=bool(maintenance),
                local_bistability_supported=coexistence, initiation_supported=initiation, phase=phase)


def worker(args):
    root, job, p = args; root = Path(root)
    folder = root/job['job_key']; folder.mkdir(exist_ok=True)
    ph = digest(root/'protocol.json')
    if (folder/'result.json').exists():
        result = read(folder/'result.json')
        if result['protocol_sha256'] != ph or digest(folder/'paths.npz') != result['paths_sha256']:
            raise ValueError('Changed completed frozen evidence')
        return result
    with np.load(job['source']) as z:
        delta, masses = validate_graph(z['delta'], z['volumes'])
        pattern = z['trajectories'][0, -1].copy()
        ids = z['ids'].copy()
    fun = rhs(delta, job['beta'], job['da'], job['db'])
    jac = lambda t, y: chemical_jacobian(y.reshape(2, -1), delta, job['beta'], job['da'], job['db'])
    starts = [pattern]
    labels = ['pattern_control']
    for family, base, scale in (('pattern', pattern, p['pattern_log_noise']),
                                ('uniform', np.ones_like(pattern), p['uniform_log_noise'])):
        for seed in p['noise_seeds']:
            starts.append(perturb(base, masses, scale, np.random.default_rng(seed)))
            labels.append(f'{family}_seed-{seed}')
    rows, payload = [], dict(delta=delta, masses=masses, ids=ids, initial_states=np.array(starts))
    started = time.perf_counter()
    for i, (label, initial) in enumerate(zip(labels, starts)):
        try:
            for horizon in p['horizons']:
                times = np.arange(0., horizon+p['interval']/2, p['interval'])
                path, error = checked_solve(fun, jac, initial, times)
                row = dict(label=label, horizon=horizon, solver_log_error=error,
                           solver_pass=bool(error < p['criteria']['solver_log_max']),
                           **metrics(path[-1], fun, jac, p['criteria']))
                if row['endpoint_kind'] != 'unresolved':
                    break
            row['pattern_return_log_rms'] = (float(distance(path[-1], payload['path_0'][-1], masses))
                if i and 'path_0' in payload else 0. if not i else None)
            payload[f'path_{i}'] = path; payload[f'times_{i}'] = times
        except (RuntimeError, FloatingPointError, ValueError) as error:
            row = dict(label=label, solver_pass=False, endpoint_kind='unresolved', error=str(error))
        rows.append(row)
        write_json(folder/'status.json', dict(state='running', completed=len(rows), total=5))
    spectrum = frozen_spectrum(delta, masses, job['beta'], job['da'], job['ratio'])
    result = dict(seed=job['seed'], branch=job['branch'], ratio=job['ratio'], db=job['db'],
        spectrum=spectrum, trials=rows, wall_seconds=time.perf_counter()-started,
        protocol_sha256=ph, **classify(rows, spectrum['uniform_jacobian_max_real'], p['criteria']))
    temporary = folder/'paths.tmp.npz'
    np.savez_compressed(temporary, **payload); temporary.replace(folder/'paths.npz')
    result['paths_sha256'] = digest(folder/'paths.npz')
    write_json(folder/'result.json', result)
    write_json(folder/'status.json', dict(state='completed', numerical_pass=result['numerical_pass']))
    return result


def assess(root):
    root = Path(root); p = read(root/'protocol.json'); verify(p)
    results = []
    for job in p['jobs']:
        file = root/job['job_key']/'result.json'
        if file.exists():
            r = read(file)
            if r['protocol_sha256'] != digest(root/'protocol.json') or digest(file.with_name('paths.npz')) != r['paths_sha256']:
                raise ValueError('Changed frozen result evidence')
            results.append(r)
    summary = dict(independent_histories=3, completed=len(results), total=len(p['jobs']),
        numerical_pass=bool(len(results) == len(p['jobs']) and all(r['numerical_pass'] for r in results)),
        unresolved_jobs=[dict(seed=r['seed'], branch=r['branch'], ratio=r['ratio']) for r in results if r['phase'] == 'unresolved'],
        uniform_thresholds=[{k:g[k] for k in ('seed', 'branch', 'uniform_instability_ratio')} for g in p['graphs']],
        results=results, protocol_sha256=digest(root/'protocol.json'))
    write_json(root/'summary.json', summary)
    if len(results) == len(p['jobs']):
        candidates = [r for r in p['ratios'] if r > 20. and all(
            x['spectrum']['uniform_jacobian_max_real'] > .02 and x['maintenance_supported']
            for x in results if x['ratio'] == r)]
        chosen = min(candidates) if candidates else max(p['ratios'])
        points = [dict(key='baseline', ratio=20., chi=.35),
                  dict(key='polarity_ablation', ratio=20., chi=0.),
                  dict(key='initiation_challenge', ratio=chosen, chi=.7)]
        write_json(root/'moving-selection.json', dict(points=points, criterion_satisfied=bool(candidates),
            protocol_sha256=digest(root/'protocol.json'), frozen_summary_sha256=digest(root/'summary.json'),
            interpretation='Three points in the actual two-parameter plane; no interpolation or frozen-chi surrogate. Each point gets three histories and two matched chemical starts.'))
        plot(root, p, results)
    return summary


def plot(root, p, results):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    phases = ['uniform_from_sampled_starts', 'coexistence', 'initiation_and_maintenance', 'mixed_sampled_basins', 'unresolved']
    names = ['Uniform from sampled starts', 'Local coexistence', 'Initiation + maintenance', 'Mixed sampled basins', 'Unresolved']
    cmap = ListedColormap(['#d7dbe3', '#3879bd', '#db8a25', '#884ca8', '#333333'])
    norm = BoundaryNorm(np.arange(-.5, 5.5), cmap.N)
    matrix = np.array([[phases.index(next(r['phase'] for r in results if r['seed'] == g['seed'] and
        r['branch'] == g['branch'] and r['ratio'] == ratio)) for ratio in p['ratios']] for g in p['graphs']])
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), gridspec_kw=dict(height_ratios=[1, 1.4]))
    for g in p['graphs']:
        rows = sorted((r for r in results if r['seed'] == g['seed'] and r['branch'] == g['branch']), key=lambda r:r['ratio'])
        axes[0].plot([r['ratio'] for r in rows], [r['spectrum']['uniform_jacobian_max_real'] for r in rows],
                     marker='.', label=f"History {g['seed']} / {'feedback on' if g['branch']=='switch_on' else 'feedback off'}")
    axes[0].axhline(0, color='black', linewidth=.8)
    axes[0].set(xlabel='D_b / D_a (D_a = 0.02)', ylabel='Largest uniform growth rate',
                title='Frozen t=150 endpoint graphs: three developmental histories')
    axes[0].legend(fontsize=8, ncol=2)
    im = axes[1].imshow(matrix, cmap=cmap, norm=norm, aspect='auto', interpolation='nearest')
    axes[1].set_xticks(range(len(p['ratios'])), [f'{r:g}' for r in p['ratios']])
    axes[1].set_yticks(range(len(p['graphs'])), [f"History {g['seed']} / {g['branch']}" for g in p['graphs']])
    axes[1].set(xlabel='D_b / D_a (sampled values)', title='Finite-horizon local chemical basin tests; polarity contrast is absent on a frozen graph')
    bar = fig.colorbar(im, ax=axes[1], ticks=range(5), fraction=.045, pad=.02)
    bar.ax.set_yticklabels(names, fontsize=8)
    fig.tight_layout(); fig.savefig(root/'frozen-phase-map.png', dpi=180); plt.close(fig)


def run(root, workers=4):
    root = Path(root); p = read(root/'protocol.json'); verify(p)
    if not 1 <= workers <= 6:
        raise ValueError('Use 1–6 CPU workers with one BLAS thread each')
    with (root/'coordinator.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        completed = sum((root/j['job_key']/'result.json').exists() for j in p['jobs'])
        write_json(root/'status.json', dict(state='running', completed=completed, total=len(p['jobs'])))
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(worker, (root, j, p)) for j in p['jobs']]
            for future in as_completed(futures):
                r = future.result()
                completed = sum((root/j['job_key']/'result.json').exists() for j in p['jobs'])
                write_json(root/'status.json', dict(state='running', completed=completed, total=len(p['jobs'])))
                print(f"{r['seed']} {r['branch']} ratio={r['ratio']:g}: {r['phase']} ({completed}/{len(p['jobs'])})", flush=True)
        summary = assess(root)
        write_json(root/'status.json', dict(state='completed', completed=len(p['jobs']), total=len(p['jobs']),
                                          numerical_pass=summary['numerical_pass']))
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'run', 'assess'))
    parser.add_argument('--output', type=Path, default=Path('outputs/parameter-robustness'))
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    result = prepare(args.output) if args.action == 'prepare' else run(args.output, args.workers) if args.action == 'run' else assess(args.output)
    print({k:v for k,v in result.items() if k not in ('jobs', 'graphs', 'results', 'source_sha256', 'input_sha256')})
