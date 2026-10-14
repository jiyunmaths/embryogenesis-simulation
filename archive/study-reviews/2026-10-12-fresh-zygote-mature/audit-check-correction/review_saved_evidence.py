"""Independent saved-array audit of mature clocks; no embryo/GPU imports.

Recompute conservative generalized eigenproblems, weighted empirical CDFs,
homologous errors, and every native density field. Earlier independent review
algorithms are reused; production comparison routines are never called.
"""
from datetime import datetime, timezone
from pathlib import Path
import importlib.util
import json
import time

import numpy as np

ROOT = Path('output/studies/fresh-zygote-mature-clock-fix')
OLD = Path('output/studies/fresh-zygote-mature-qualification')
HERE = Path(__file__).resolve().parent
EXPECTED = '057f18713e37c562cc70743dd0dbff37be6e76b53040455a4535452752c19b04'
HELPER = Path('archive/study-reviews/2026-10-10-development-time-diagnosis/review_saved_arrays.py')
BRANCHES = ('off', 'tension_adhesion', 'polarity', 'full')
COUPLINGS = {'off': (False, 0., 0., 0.), 'tension_adhesion': (True, .25, .35, 0.),
             'polarity': (True, 0., 0., .35), 'full': (True, .25, .35, .35)}
LEVELS = {'coarse': .00375, 'fine': .001875}


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


t = module('independent_dense_review', HELPER)
r = t.r


def compare(a, b):
    """Generalize the independently written CDF comparison to 41 observations."""
    maxima, homologous, peaks = {}, {}, {}
    assert len(a) == len(b) == 41
    for x, y in zip(a, b):
        assert abs(x['time'] - y['time']) < 1e-9
        assert x['ids'] == y['ids'] and x['parents'] == y['parents']
        u, w = np.asarray(x['chemistry']), np.asarray(y['chemistry'])
        v, z = np.asarray(x['volumes']), np.asarray(y['volumes'])
        sx, sy = r.eig(x)[0], r.eig(y)[0]
        values = dict(
            concentration_wasserstein_max=max(r.w1(i, j, v, z) for i, j in zip(u, w)),
            log_wasserstein_max=max(r.w1(np.log(i), np.log(j), v, z) for i, j in zip(u, w)),
            contrast_abs_max=float(abs(u.std(1) - w.std(1)).max()),
            relative_axis_max=abs(x['axis_ratio'] / y['axis_ratio'] - 1),
            relative_total_volume_max=abs(v.sum() / z.sum() - 1),
            whole_amount_relative_max=float(abs((u @ v) / (w @ z) - 1).max()),
            mean_polarity_abs_max=float(abs(v @ np.array(x['polarity']) / v.sum()
                - z @ np.array(y['polarity']) / z.sum()).max()),
            relative_spectrum_max=float(np.linalg.norm(sx - sy) / np.linalg.norm(sy)),
            cell_count_difference_max=0., unmatched_observations_max=0)
        for name, value in values.items():
            if name not in maxima or value > maxima[name]:
                maxima[name] = float(value)
                peaks[name] = x['time']
        values = dict(chemical_log_max=float(abs(np.log(u / w)).max()),
            polarity_abs_max=float(abs(np.array(x['polarity']) - np.array(y['polarity'])).max()),
            relative_volume_max=float(abs(v / z - 1).max()),
            relative_transport_max=float(np.linalg.norm(np.array(x['delta']) - np.array(y['delta']))
                / np.linalg.norm(np.array(y['delta']))), unmatched_observations_max=0)
        for name, value in values.items():
            homologous[name] = max(homologous.get(name, 0.), float(value))
    return maxima, homologous, peaks


def invariant(z):
    metadata = json.loads(str(z['metadata']))
    for key in ('time', 'step_number', 'volume_projection_max', 'clipped_fraction'):
        metadata.pop(key, None)
    return dict(metadata=metadata, arrays={name: (str(z[name].dtype), z[name].shape, z[name].tobytes())
        for name in ('target', 'fate', 'ids', 'parents', 'due', 'rounding', 'precision_arm',
                     'chemical_method', 'attribute_mode', 'attribute_schema',
                     'development_schema', 'event_audits')})


def endpoint(row, graph, cfg):
    lam = r.eig(graph)[0]
    beta = cfg['signal_beta']
    jac = np.array([[1., -1.], [2 * beta, -beta]])
    growth = [np.linalg.eigvals(jac - ell * np.diag([cfg['signal_da'], cfg['signal_dh']])).real.max()
              for ell in lam]
    return dict(time=row['time'], log_activator_sd=float(np.log(np.array(row['chemistry'])[0]).std()),
        axis_ratio=row['axis_ratio'], summed_undirected_conductance=float(np.triu(graph['conductance'], 1).sum()),
        frozen_uniform_unstable_modes=int(np.sum((lam > 1e-10) & (np.array(growth) > 1e-10))))


def main():
    started = time.perf_counter()
    output = HERE / 'independent_review.json'
    assert not output.exists(), 'Never overwrite a completed review'
    for file in (Path(__file__), HELPER, t.HELPER):
        r.pin(file)
    p = r.read(ROOT / 'protocol.json')
    ph = r.pin(ROOT / 'protocol.json', EXPECTED)
    status, summary = r.read(ROOT / 'status.json'), r.read(ROOT / 'summary.json')
    assert status['state'] == 'completed_provisional'
    assert status['completed'] == status['total'] == 24
    assert status['pairs_completed'] == status['pairs_total'] == 12
    assert summary['execution_completed'] and summary['mature_screens_passed']
    assert not summary['qualification_accepted'] and summary['independent_review_pending']
    assert summary['protocol_sha256'] == status['protocol_sha256'] == ph
    assert summary['new_scientific_histories'] == 0 and summary['new_numerical_paths'] == 24
    assert not summary['long_scientific_continuation_started']
    print('Checking 58 source and 494 input pins', flush=True)
    for name, expected in {**p['source_sha256'], **p['input_sha256']}.items():
        r.pin(name, expected)
    parent = r.read('output/studies/fresh-zygote-control-completion/protocol.json')
    r.pin('output/studies/fresh-zygote-control-completion/protocol.json', p['parent_protocol_sha256'])
    for key in ('physical_limits', 'time_limits', 'dense_time_spectrum_limit', 'prefix_limits',
                'prefix_steps', 'method', 'gpu_uuid', 'runtime_binary_sha256'):
        assert p[key] == parent[key], key
    accepted = r.read('docs/fresh_zygote_control_completion_assessment.json')
    assert accepted['scoped_combined_developmental_preparations_accepted']
    prior = r.read(accepted['independent_review']['path'])
    r.pin(accepted['independent_review']['path'], accepted['independent_review']['sha256'])
    assert prior['scoped_combined_developmental_preparations_accepted']
    original = r.read(OLD / 'protocol.json')
    r.pin(OLD / 'protocol.json', 'f6c099718c3d368646912660fb91fc2ce1280ce11130fb105cf07790f9811131')
    assert r.read(OLD / 'status.json')['state'] == 'failed'
    assert r.read('output/studies/fresh-zygote-pilot/status.json')['state'] == 'failed'
    assert not r.read('output/studies/fresh-zygote-pilot/time-h1103-off.json')['passed']
    for key, value in original.items():
        if key not in ('jobs', 'study_root', 'source_sha256', 'input_sha256', 'stop_policy'):
            assert p[key] == value, key
    assert p['histories'] == [1101, 1102, 1103] and p['branches'] == list(BRANCHES)
    assert p['start_time'] == 18. and p['duration'] == 6. and p['snapshots'] == [0, 3, 6]
    assert p['interval'] == .15 and p['diagnostic_interval'] == .015
    assert p['restart_elapsed'] == 3. and p['restart_steps'] == 8
    runtime = r.read(ROOT / 'runtime.json')
    r.pin(runtime['binary'], p['runtime_binary_sha256'])
    assert runtime['binary_sha256'] == p['runtime_binary_sha256']
    expected_jobs = [f'h{h}_{b}_{l}' for h in p['histories'] for b in BRANCHES for l in LEVELS]
    assert [j['key'] for j in p['jobs']] == expected_jobs
    histories, networks, reports, paths, prefixes = {}, {}, {}, [], []
    limits = p['physical_limits']
    for job, old_job in zip(p['jobs'], original['jobs']):
        key, cfg = job['key'], job['config']
        folder = ROOT / key
        assert job == dict(old_job, source=str(Path(p['study_root']) / key / 'source.npz'))
        assert Path(job['source']).is_symlink()
        assert Path(job['source']).resolve() == Path(old_job['source']).resolve()
        assert cfg['dt'] == LEVELS[job['level']] and cfg['grid'] == 120 and cfg['extent'] == 2.24
        assert cfg['steps'] * cfg['dt'] == 24. and cfg['save_every'] * cfg['dt'] == .15
        assert cfg['seed'] == job['seed'] and cfg['signal_transport'] == 'conservative'
        assert tuple(cfg[k] for k in ('feedback', 'fate_tension', 'fate_adhesion', 'polarity_tension')) == COUPLINGS[job['branch']]
        report = r.read(folder / 'gpu-result.json')
        assert report['passed'] and report['failure'] is None and report['job'] == job
        assert report['protocol_sha256'] == ph
        assert report['final_time'] == 24. and report['final_cells'] == 16 and len(report['events']) == 15
        assert report['exact_restart_steps'] == 8 and report['lineage_cycle_rng_unchanged']
        audit = report['audit']
        assert audit['max_volume_error'] < limits['volume_error_max'] and audit['min_radius'] >= limits['radius_min']
        assert audit['max_clipping'] == 0 and audit['dilution_amount_error'] <= limits['dilution_amount_error_max']
        for name, expected in report['evidence_sha256'].items():
            r.pin(name, expected)
        history = r.read(folder / 'gpu-history.json')
        trace = t.jsonl(folder / 'network.jsonl')
        assert len(history) == 41 and len(trace) == report['network_trace_rows'] == 401
        closure, boundary = 0., 0.
        for i, row in enumerate(history):
            r.close(row['time'], 18. + i * .15, 'regular time')
            r.close(row['elapsed'], i * .15, 'regular elapsed')
            assert len(row['ids']) == 16 and not row['dividing_ids']
            assert row['nondividing_components'] == [1] * 16
            assert row['boundary_occupancy'] < limits['boundary_max']
            u, v = np.asarray(row['chemistry']), np.asarray(row['volumes'])
            assert u.shape == (2, 16) and np.isfinite(u).all() and u.min() > 0
            assert np.isfinite(v).all() and v.min() > 0
            r.close(row['log_activator_sd'], np.log(u[0]).std(), 'observed contrast')
            closure = max(closure, r.eig(row)[1])
            boundary = max(boundary, row['boundary_occupancy'])
        for i, row in enumerate(trace):
            r.close(row['time'], 18. + i * .015, 'dense time')
            r.close(row['time'], row['step'] * cfg['dt'], 'dense clock')
            assert row['ids'] == history[0]['ids'] and row['parents'] == history[0]['parents']
            assert not row['active_divisions'] and row['completed_events'] == 15
            assert row['last_step_projection_max'] == 0
            closure = max(closure, t.graph(row, cfg)[1])
        prefix = r.read(folder / 'actual-start-prefix.json')
        assert prefix['passed'] and prefix['failure'] is None and prefix['steps'] == 8 and prefix['exact_gpu_restart']
        assert prefix['protocol_sha256'] == ph
        assert prefix == next(x for x in summary['prefixes'] if x['case'] == key)
        for name, value in prefix['errors'].items():
            assert value <= p['prefix_limits'][name]
        for name, expected in prefix['evidence_sha256'].items():
            r.pin(name, expected)
        if key == p['reused_prefix_case']:
            assert prefix['reused'] and prefix['original_protocol_sha256'] == r.pin(OLD / 'protocol.json')
            r.pin(prefix['original_report'], prefix['original_report_sha256'])
        prefixes.append(dict(case=key, recorded_passed=True, errors=prefix['errors'],
                             saved_evidence_hashed=True, independently_rerun=False))
        histories[key], networks[key], reports[key] = history, trace, report
        paths.append(dict(case=key, audit=audit, boundary_max=boundary, graph_closure_max=closure,
                          regular_rows=41, dense_rows=401, recorded_exact_restart_steps=8))
        print('Independently checked histories, graph closure, provenance:', key, flush=True)
    pairs, samples, endpoints = [], [], []
    jobs = {j['key']: j for j in p['jobs']}
    for pair in p['time_pairs']:
        a, b = pair['cases']
        saved = r.read(ROOT / f"time-h{pair['seed']}-{pair['branch']}.json")
        assert saved == next(x for x in summary['pairs'] if x['cases'] == [a, b])
        assert saved['passed'] and saved['protocol_sha256'] == ph
        for name, expected in saved['evidence_sha256'].items():
            r.pin(name, expected)
        errors, homologous, peaks = compare(histories[a], histories[b])
        for name, value in errors.items():
            r.close(value, saved['errors'][name], (a, name))
            assert value <= p['time_limits'][name]
        for name, value in homologous.items():
            r.close(value, saved['homologous_errors'][name], (a, name))
            assert value <= p['homologous_limits'][name]
        dense, topology = [], []
        for x, y in zip(networks[a], networks[b]):
            assert abs(x['time'] - y['time']) < 1e-9
            sx, sy = r.eig(x)[0], r.eig(y)[0]
            dense.append(dict(time=x['time'], relative_spectrum=float(np.linalg.norm(sx - sy) / np.linalg.norm(sy))))
            if not np.array_equal(np.array(x['conductance']) > 0, np.array(y['conductance']) > 0):
                topology.append(x['time'])
        peak = max(dense, key=lambda x: x['relative_spectrum'])
        r.close(peak['relative_spectrum'], saved['dense_peak']['relative_spectrum'], 'dense maximum')
        assert peak['time'] == saved['dense_peak']['time'] and topology == saved['topology_mismatch_times']
        assert peak['relative_spectrum'] <= p['dense_time_spectrum_limit']
        fields, constants = [], {}
        for elapsed in (0, 3, 6):
            densities = {}
            for key in (a, b):
                cfg = jobs[key]['config']
                z = r.load(ROOT / key / f'gpu-e{elapsed:02d}.npz')
                observed = histories[key][round(elapsed / .15)]
                check = r.checkpoint(z, cfg, observed, False)
                metadata = json.loads(str(z['metadata']))
                assert check['cells'] == 16 and check['active'] == 0
                assert metadata['time'] == 18. + elapsed and metadata['step_number'] * cfg['dt'] == metadata['time']
                assert metadata['next_id'] == 31 and len(metadata['lineage']) == 31
                assert json.loads(str(z['event_audits'])) == reports[key]['events']
                assert z['rounding'].dtype == np.int32 and z['rounding'].shape == (2,)
                assert check['volume_error'] < limits['volume_error_max'] and check['minimum_radius'] >= limits['radius_min']
                if elapsed == 0:
                    source = r.load(jobs[key]['source'])
                    assert set(source) == set(z)
                    for name in z:
                        assert z[name].dtype == source[name].dtype and z[name].shape == source[name].shape
                        assert z[name].tobytes() == source[name].tobytes(), (key, 'initial inheritance', name)
                    del source
                    constants[key] = invariant(z)
                else:
                    assert invariant(z) == constants[key], (key, 'lineage/target/RNG drift')
                if elapsed == 3:
                    # Identical complete archive bytes establish all replay-start arrays match.
                    assert r.pin(ROOT / key / 'restart.npz') == r.pin(ROOT / key / 'gpu-e03.npz')
                densities[key] = r.density(z)
                samples.append(dict(case=key, **check, inherited_constants_exact=True))
                del z
            values = r.field_error(densities[a], densities[b])
            expected = next(x for x in saved['fields'] if x['elapsed'] == elapsed)
            for name, value in values.items():
                r.close(value, expected['errors'][name], ('density', a, elapsed, name))
                assert value <= p['time_limits']['field_max']
            fields.append(dict(elapsed=elapsed, time=18. + elapsed, errors=values))
            del densities
            print('Independently reconstructed both checkpoint states and densities:', a, elapsed, flush=True)
        field_max = max(v for x in fields for v in x['errors'].values())
        r.close(field_max, saved['field_max'], 'maximum native field')
        pairs.append(dict(history=pair['seed'], branch=pair['branch'], cases=[a, b], independently_passed=True,
                          errors=errors, homologous_errors=homologous, peak_times=peaks,
                          dense_peak=peak, topology_mismatch_times=topology, fields=fields, field_max=field_max))
        first, last = [endpoint(histories[b][i], networks[b][i * 10], jobs[b]['config']) for i in (0, 40)]
        endpoints.append(dict(history=pair['seed'], branch=pair['branch'], dt=LEVELS['fine'],
            start=first, end=last, log_contrast_growth=last['log_activator_sd'] / first['log_activator_sd'],
            conductance_relative_change=last['summed_undirected_conductance'] / first['summed_undirected_conductance'] - 1))
    assert len(samples) == 72 and len(pairs) == len(endpoints) == 12
    print('Checking that every pinned input and production record is unchanged', flush=True)
    for name, expected in list(r.WATCHED.items()):
        assert r.sha(name) == expected, name
    result = dict(review_completed=True, scoped_mature_timestep_accepted=True,
        reviewed_utc=datetime.now(timezone.utc).isoformat(), seconds=time.perf_counter() - started,
        protocol_sha256=ph, source_count=len(p['source_sha256']), input_count=len(p['input_sha256']),
        histories=p['histories'], independent_histories=3, new_scientific_histories=0,
        new_numerical_paths=24, paths=paths, prefixes=prefixes, pairs=pairs,
        checkpoints=samples, fine_endpoint_observations=endpoints,
        original_scheduling_failure_retained=True, original_developmental_failure_retained=True,
        full_cpu_horizon_rerun=False, exact_restart_independently_rerun=False,
        long_window_qualification_accepted=False, biological_identity_qualified=False,
        general_developmental_convergence_accepted=False, long_scientific_continuation_started=False,
        scope='Independent NumPy/SciPy saved-array audit: 24 histories (984 regular observations), '
              '9624 dense graphs, all 12 timestep comparisons, 72 complete checkpoint states/native density '
              'reconstructions, exact initial-state inheritance and unchanged lineage/targets/RNG. '
              'All 24 prefixes and eight-step restarts checked through immutable saved evidence; no fresh '
              'CPU/GPU stepping. Restart-start archives are byte-identical to their t21 checkpoints. '
              'Acceptance applies only to dt0.00375/0.001875, grid120, physical t18–24, these twelve contexts. '
              'No long-window formation, persistence, biological identity or universal convergence claim.',
        original_files_sha256=dict(r.WATCHED))
    output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ('review_completed', 'scoped_mature_timestep_accepted', 'seconds')}, indent=2), flush=True)


if __name__ == '__main__':
    alias = module('verified_study_alias', Path('scripts/with_study_output_path_alias.py'))
    with alias.base.validated_alias(Path('archive/output-relocation/2026-10-10-output/studies-manifest.json'),
            '8fb2545d8bae35a1a673ec3a43d3cc191148f03cac58892f29bd368182a8aff6',
            source_root=alias.SOURCE_ROOT, target_root=alias.TARGET_ROOT):
        main()
