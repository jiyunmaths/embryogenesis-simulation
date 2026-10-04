"""Matched finite-horizon frozen references for moving parameter spot checks.

Near-critical linear growth may need longer than the moving 60-unit window.
Compare like-for-like horizons before attributing a missing pattern to motion.
"""
import argparse
from pathlib import Path

import numpy as np

from .attribute_development import AttributeSimulation
from .attribute_persistence import rhs
from .cell_response_exchange import checked_solve
from .feedback_endpoint_bistability import chemical_jacobian
from .feedback_long import digest
from .neighbor_context import read, validate_graph
from .parameter_robustness import verify, frozen_spectrum
from .resolution import write_json, _steps


def prepare(root, moving):
    root, moving = Path(root).resolve(), Path(moving).resolve()
    if root.exists(): raise FileExistsError(root)
    p = read(moving/'protocol.json'); verify(p)
    inputs = {str(moving/'protocol.json'):digest(moving/'protocol.json')}
    root.mkdir(parents=True)
    graphs = {}
    for job in p['jobs']:
        if job['seed'] in graphs: continue
        host = AttributeSimulation.restore(job['source']); graph = host.signaling_graph()
        file = root/f"seed-{job['seed']}.npz"
        with np.load(job['chemical_file']) as z:
            if not np.array_equal(z['ids'], host.ids): raise ValueError('Cell order mismatch')
            np.savez_compressed(file, delta=graph.delta, masses=graph.masses, ids=host.ids,
                                pattern=z['pattern'], uniform=z['uniform'])
        graphs[job['seed']] = str(file)
        for f in (Path(job['source']), Path(job['chemical_file']), file): inputs[str(f)] = digest(f)
    protocol = dict(moving=str(moving), parent_protocol_sha256=digest(moving/'protocol.json'),
        jobs=[dict(**job, graph=graphs[job['seed']]) for job in p['jobs']], duration=p['duration'],
        interval=p['interval'], late_window=p['late_window'], late_log_sd_min=p['criteria']['late_log_sd_min'],
        independent_histories=3,
        interpretation='Same initial geometry, chemistry, parameters and observation horizon as each moving run. Polarity contrast is absent in the frozen reference. Trajectory contrast is descriptive, not endpoint equilibrium. Missing contrast in both references cannot establish suppression by moving geometry.',
        source_sha256={str(Path(__file__).resolve()):digest(__file__), **p['source_sha256']}, input_sha256=inputs)
    write_json(root/'protocol.json', protocol)
    return protocol


def run(root):
    root = Path(root); p = read(root/'protocol.json'); verify(p)
    rows = []
    times = np.arange(_steps(p['duration'], p['interval'])+1)*p['interval']
    for job in p['jobs']:
        with np.load(job['graph']) as z:
            delta, masses = validate_graph(z['delta'], z['masses']); initial = z[job['family']].copy()
        fun = rhs(delta, 2., .02, .02*job['point']['ratio'])
        jac = lambda t,y:chemical_jacobian(y.reshape(2,-1), delta, 2., .02, .02*job['point']['ratio'])
        path, error = checked_solve(fun, jac, initial, times)
        spread = np.std(np.log(path[:,0]), axis=1)
        late = spread[times >= p['duration']-p['late_window']-1e-9]
        file = root/(job['key']+'.npz')
        np.savez_compressed(file, times=times, trajectory=path, initial=initial, delta=delta, masses=masses)
        rows.append(dict(key=job['key'], seed=job['seed'], family=job['family'], point=job['point'],
            solver_log_error=error, final_log_sd=float(spread[-1]), late_min_log_sd=float(late.min()),
            persistent_contrast=bool(late.min() > p['late_log_sd_min']), paths_sha256=digest(file),
            spectrum=frozen_spectrum(delta, masses, 2., .02, job['point']['ratio'])))
    result = dict(completed=len(rows), total=len(p['jobs']), independent_histories=3,
        protocol_sha256=digest(root/'protocol.json'), max_solver_log_error=max(r['solver_log_error'] for r in rows),
        results=rows,
        horizon_limited_initiation=[r['key'] for r in rows if r['family']=='uniform' and
            r['spectrum']['uniform_jacobian_max_real'] > 0 and not r['persistent_contrast']])
    write_json(root/'results.json', result)
    return result


def assess(root):
    root = Path(root); p = read(root/'protocol.json'); verify(p)
    r = read(root/'results.json')
    if r['protocol_sha256'] != digest(root/'protocol.json'): raise ValueError('Changed reference protocol')
    pairs = []
    for row in r['results']:
        if digest(root/(row['key']+'.npz')) != row['paths_sha256']: raise ValueError('Changed frozen reference')
        folder = Path(p['moving'])/row['key']
        if not (folder/'result.json').exists(): continue
        mr = read(folder/'result.json')
        if mr['protocol_sha256'] != p['parent_protocol_sha256'] or digest(folder/'history.json') != mr['history_sha256']:
            raise ValueError('Changed matched moving evidence')
        history = read(folder/'history.json')
        with np.load(root/(row['key']+'.npz')) as z:
            if not np.allclose([h['elapsed'] for h in history], z['times'], rtol=0, atol=1e-9):
                raise ValueError('Moving/frozen observation horizon mismatch')
            if not np.allclose(history[0]['chemistry'], z['initial'], rtol=1e-14, atol=1e-14):
                raise ValueError('Moving/frozen initial chemistry mismatch')
        pairs.append(dict(key=row['key'], seed=row['seed'], point=row['point'], family=row['family'],
            frozen_late_min_log_sd=row['late_min_log_sd'], moving_late_min_log_sd=mr['late_min_log_activator_sd'],
            frozen_contrast=row['persistent_contrast'], moving_contrast=mr['persistent_contrast'],
            frozen_horizon_too_short_for_observed_initiation=row['key'] in r['horizon_limited_initiation']))
    result = dict(completed_pairs=len(pairs), total=len(r['results']), independent_histories=3,
                  horizon_limited_initiation=r['horizon_limited_initiation'], pairs=pairs,
                  reference_results_sha256=digest(root/'results.json'))
    write_json(root/'comparison.json', result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'run', 'assess'))
    parser.add_argument('--output', type=Path, default=Path('outputs/parameter-robustness-reference'))
    parser.add_argument('--moving', type=Path, default=Path('outputs/parameter-robustness-moving'))
    args = parser.parse_args()
    result = prepare(args.output,args.moving) if args.action=='prepare' else run(args.output) if args.action=='run' else assess(args.output)
    print({k:v for k,v in result.items() if k not in ('jobs','results','pairs','source_sha256','input_sha256')})
