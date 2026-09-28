"""Timestep refinement and independent developmental histories for switch survival."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
from pathlib import Path
import time
import numpy as np
from .attribute_development import AttributeSimulation, attributes
from .model import Config
from .feedback_long import digest, save_checkpoint, restore_checkpoint
from . import feedback_survival as survival
from .resolution import write_json, _steps


def retime(source, destination, dt):
    """Preserve the physical state and random streams, including the clock."""
    sim = AttributeSimulation.restore(source)
    old_dt = sim.config.dt
    sim.step_number = _steps(sim.time, dt)
    sim.config.steps = _steps(sim.config.steps * old_dt, dt)
    sim.config.save_every = _steps(sim.config.save_every * old_dt, dt)
    sim.config.dt = dt
    sim.config.validate()
    sim.checkpoint(destination)
    return sim


def prepare(root, baseline, seeds=(8, 9)):
    root, baseline = Path(root).resolve(), Path(baseline).resolve()
    if root.exists():
        raise FileExistsError(root)
    if not seeds or len(set(seeds)) != len(seeds):
        raise ValueError('Use distinct developmental seeds')
    old = json.loads((baseline / 'protocol.json').read_text())
    if not (baseline / 'comparison.json').exists():
        raise ValueError('Requires completed coarse comparison')
    source = Path(old['checkpoint'])
    if digest(source) != old['checkpoint_sha256']:
        raise ValueError('Original developmental checkpoint changed')
    if old['source_config']['seed'] in seeds:
        raise ValueError('New histories must differ from the original seed')
    root.mkdir(parents=True)
    retime(source, root / 'refined-source.npz', old['dt'] / 2)
    survival.prepare(root / 'refined-source.npz', root / 'refined', old['duration'], old['interval'])
    inputs = [baseline / 'protocol.json', baseline / 'comparison.json', source,
              root / 'refined-source.npz']
    inputs += [baseline / arm / name for arm in survival.ARMS for name in ('history.json', 'result.json')]
    p = dict(baseline=str(baseline), seeds=list(seeds), source_config=old['source_config'],
             start=90., duration=old['duration'], interval=old['interval'], fine_dt=old['dt']/2,
             criteria=dict(max_chemical_log_rms=.02, max_relative_axis_ratio_error=.01,
                           require_same_survival_classification=True),
             scope='Two new unselected zygote histories at coarse dt; same-state timestep halving for original seed. Pilot robustness study, not population inference or spatial convergence.',
             source_sha256={str(f.resolve()): digest(f) for f in Path(__file__).parent.glob('*.py')},
             input_sha256={str(f): digest(f) for f in inputs})
    write_json(root / 'protocol.json', p)
    write_json(root / 'status.json', dict(state='prepared', protocol_sha256=digest(root/'protocol.json')))
    return p


def development(root, seed, p):
    folder = Path(root) / f'seed-{seed}' / 'development'
    folder.mkdir(parents=True, exist_ok=True)
    ph = digest(Path(root)/'protocol.json')
    job = dict(stage='development', seed=seed)
    checkpoint = folder/'latest_state.npz'
    if (folder/'result.json').exists():
        return json.loads((folder/'result.json').read_text())
    if checkpoint.exists():
        sim, audit, history = restore_checkpoint(checkpoint, ph, job)
    else:
        config = dict(p['source_config'], seed=seed, feedback=False,
                      steps=_steps(p['start'], p['source_config']['dt']))
        sim = AttributeSimulation(Config(**config), 'no_feedback')
        history = []
        audit = dict(max_volume_error=0., min_radius=1e100, max_clipping=0., max_sampled_boundary=0.)
    stop = _steps(p['start'], sim.config.dt)
    every = _steps(p['interval'], sim.config.dt)
    last = history[-1]['metrics']['time'] if history else -1.
    try:
        while sim.step_number <= stop:
            v = sim.volumes()
            audit['max_volume_error'] = max(audit['max_volume_error'], float(np.max(abs(v/sim.target-1))))
            audit['min_radius'] = min(audit['min_radius'], float(np.min((3*v/(4*np.pi))**(1/3))/sim.dx))
            audit['max_clipping'] = max(audit['max_clipping'], sim.clipped_fraction)
            if audit['max_volume_error'] >= .05 or audit['min_radius'] < 4 or audit['max_clipping'] > 0:
                raise RuntimeError('Development numerical quality limit exceeded')
            if sim.step_number % every == 0 and sim.time > last + 1e-9:
                row = attributes(sim); history.append(row); last = sim.time
                audit['max_sampled_boundary'] = max(audit['max_sampled_boundary'], row['metrics']['boundary_occupancy'])
                if audit['max_sampled_boundary'] >= .01:
                    raise RuntimeError('Development boundary quality limit exceeded')
                write_json(folder/'history.json', history)
                write_json(folder/'status.json', dict(state='running', time=sim.time, audit=audit))
                if sim.step_number % (10*every) == 0 or sim.step_number == stop:
                    save_checkpoint(sim, checkpoint, audit, history, ph, job)
                    print(f'development seed {seed}: t={sim.time:g}', flush=True)
            if sim.step_number == stop:
                break
            sim.step()
        mature = len(sim.phi) == 16 and not sim.divisions
        spread = float(np.std(np.log(sim.activator)))
        result = dict(seed=seed, quality_pass=True, mature=bool(mature), initial_log_activator_sd=spread,
                      eligible_for_survival=bool(mature and spread > .1), audit=audit,
                      checkpoint=str(checkpoint.resolve()))
        write_json(folder/'result.json', result)
        write_json(folder/'status.json', dict(state='completed', time=sim.time))
        return result
    except Exception as error:
        write_json(folder/'status.json', dict(state='failed', time=sim.time, error=str(error), audit=audit))
        raise


def compare_refinement(root, p):
    baseline, refined = Path(p['baseline']), Path(root)/'refined'
    report = {}
    for arm in survival.ARMS:
        coarse = json.loads((baseline/arm/'history.json').read_text())
        fine = json.loads((refined/arm/'history.json').read_text())
        if len(coarse) != len(fine):
            raise ValueError('Unmatched observation counts')
        errors, shapes = [], []
        for a, b in zip(coarse, fine):
            if abs(a['metrics']['time']-b['metrics']['time']) > 1e-9 or a['cell_order'] != b['cell_order']:
                raise ValueError('Unmatched times or cell identities')
            def chemical(row):
                return np.array([[c['attributes'][k] for c in row['cells']] for k in ('log_activator','log_inhibitor')])
            weights = np.array([c['context']['volume'] for c in coarse[0]['cells']])
            weights /= weights.sum()
            errors.append(float(np.sqrt(np.sum(weights*(chemical(a)-chemical(b))**2)/2)))
            shapes.append(abs(a['metrics']['axis_ratio']/b['metrics']['axis_ratio']-1))
        a = json.loads((baseline/arm/'result.json').read_text())
        b = json.loads((refined/arm/'result.json').read_text())
        same = all(a[k] == b[k] for k in ('contrast_survives','initial_ordering_retained'))
        gates = dict(chemical=max(errors)<=p['criteria']['max_chemical_log_rms'],
                     shape=max(shapes)<=p['criteria']['max_relative_axis_ratio_error'],
                     classification=same, quality=a['quality_pass'] and b['quality_pass'])
        report[arm] = dict(max_chemical_log_rms=max(errors), max_relative_axis_ratio_error=max(shapes),
                           gates=gates, passed=all(gates.values()), coarse=a, fine=b)
    write_json(Path(root)/'refinement.json', report)
    return report


def refinement(root, p):
    survival.run(Path(root)/'refined')
    return compare_refinement(root, p)


def cohort(root, p):
    rows = []
    for seed in p['seeds']:
        try:
            row = development(root, seed, p)
            if row['eligible_for_survival']:
                output = Path(root)/f'seed-{seed}'/'survival'
                if not output.exists():
                    survival.prepare(row['checkpoint'], output, p['duration'], p['interval'])
                survival.run(output)
                if json.loads((output/'status.json').read_text())['state'] != 'completed':
                    raise RuntimeError('Paired survival failed; see arm status')
                row['survival'] = json.loads((output/'comparison.json').read_text())
            else:
                row['outcome'] = 'no_pattern_at_switch' if row['mature'] else 'development_incomplete'
        except Exception as error:
            row = dict(seed=seed, outcome='failed', error=str(error))
        rows.append(row)
        write_json(Path(root)/'cohort.json', dict(planned_seeds=p['seeds'], histories=rows))
    return rows


def run(root):
    root = Path(root); p = json.loads((root/'protocol.json').read_text())
    ph = digest(root/'protocol.json')
    if ph != json.loads((root/'status.json').read_text())['protocol_sha256']:
        raise ValueError('Protocol changed')
    for name, expected in {**p['source_sha256'], **p['input_sha256']}.items():
        if digest(name) != expected:
            raise ValueError('Frozen input changed: '+name)
    write_json(root/'status.json', dict(state='running', protocol_sha256=ph, started_unix=time.time()))
    results, failures = {}, []
    # At most four new mechanics processes: two refinement arms and a cohort
    # pipeline using one development worker or two paired survival workers.
    with ProcessPoolExecutor(max_workers=2) as pool:
        jobs = {pool.submit(fn, root, p): name for name, fn in [('refinement',refinement),('cohort',cohort)]}
        for future in as_completed(jobs):
            name = jobs[future]
            try:
                results[name] = future.result()
                if name == 'cohort' and any(r.get('outcome') == 'failed' for r in results[name]):
                    failures.append(dict(stage=name, error='One or more planned histories failed'))
            except Exception as error:
                failures.append(dict(stage=name, error=str(error)))
            write_json(root/'status.json', dict(state='running', protocol_sha256=ph, finished_stages=list(results), failures=failures))
    write_json(root/'summary.json', dict(results=results, failures=failures,
        original_history=json.loads((Path(p['baseline'])/'comparison.json').read_text()),
        planned_total_histories=1+len(p['seeds']),
        scope=p['scope']))
    write_json(root/'status.json', dict(state='failed' if failures else 'completed', protocol_sha256=ph,
                                      finished_stages=list(results), failures=failures))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['prepare','run'])
    parser.add_argument('--output', type=Path, default=Path('outputs/feedback-survival-validation'))
    parser.add_argument('--baseline', type=Path, default=Path('outputs/feedback-survival'))
    args = parser.parse_args()
    if args.command == 'prepare': prepare(args.output, args.baseline)
    else: run(args.output)
