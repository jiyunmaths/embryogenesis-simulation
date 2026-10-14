"""Independent NumPy/SciPy formation audit; no solver/GPU imports or stepping.

Reuse earlier independent generalized-eigenproblem, CDF and checkpoint audits.
Keep production files and all prior failures immutable; write a separate review.
"""
from datetime import datetime, timezone
from pathlib import Path
import importlib.util
import json
import time

import numpy as np

ROOT = Path('output/studies/fresh-zygote-formation')
PARENT = Path('output/studies/fresh-zygote-mature-clock-fix')
HERE = Path(__file__).resolve().parent
EXPECTED = 'e193c41f6d018cf4cfc7ed88777e00f7e62e99aa89fb98e783baab0c8e199d06'
HELPER = Path('archive/study-reviews/2026-10-12-fresh-zygote-mature/review_saved_evidence.py')


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


prior = module('independent_mature_audit', HELPER)
t, r = prior.t, prior.r


def growth(row, cfg):
    lam = r.eig(row)[0]
    beta = cfg['signal_beta']
    jac = np.array([[1., -1.], [2 * beta, -beta]])
    rates = np.array([np.linalg.eigvals(jac - ell * np.diag([cfg['signal_da'], cfg['signal_dh']])).real.max()
                      for ell in lam[lam > 1e-10]])
    return dict(max_growth=float(rates.max()), unstable_modes=int(np.sum(rates > 1e-10)))


def outcomes(history, p):
    times = np.array([row['elapsed'] for row in history])
    contrast = np.array([np.log(np.array(row['chemistry'])[0]).std() for row in history])
    results = []
    for threshold in p['thresholds']:
        above = contrast > threshold
        crossing = np.flatnonzero(above)
        onset = None
        for index in crossing:
            end = times[index] + p['sustained_window']
            if end > times[-1] + 1e-9:
                break
            interval = (times >= times[index] - 1e-9) & (times <= end + 1e-9)
            if above[interval].all():
                onset = float(times[index])
                break
        late = times >= p['duration'] - p['late_window'] - 1e-9
        changes = np.flatnonzero(above[1:] != above[:-1]) + 1
        results.append(dict(threshold=threshold,
            first_crossing_elapsed=float(times[crossing[0]]) if len(crossing) else None,
            first_sustained_onset_elapsed=onset, sustained_formation=onset is not None,
            late_persistent=bool(above[late].all()), late_min_log_sd=float(contrast[late].min()),
            final_log_sd=float(contrast[-1]), sampled_transitions=[dict(elapsed=float(times[i]),
                direction='up' if above[i] else 'down') for i in changes]))
    return results


def same_outcomes(independent, recorded):
    assert len(independent) == len(recorded) == 3
    for x, y in zip(independent, recorded):
        assert x.keys() == y.keys()
        for name, value in x.items():
            if isinstance(value, float):
                r.close(value, y[name], ('formation', name))
            else:
                assert value == y[name], (name, value, y[name])


def compare(a, b, cfg):
    assert len(a) == len(b) == 401
    maxima, homologous, peaks = {}, {}, {}
    growth_difference = 0.
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
                maxima[name], peaks[name] = float(value), x['time']
        values = dict(chemical_log_max=float(abs(np.log(u / w)).max()),
            polarity_abs_max=float(abs(np.array(x['polarity']) - np.array(y['polarity'])).max()),
            relative_volume_max=float(abs(v / z - 1).max()),
            relative_transport_max=float(np.linalg.norm(np.array(x['delta']) - np.array(y['delta']))
                / np.linalg.norm(np.array(y['delta']))), unmatched_observations_max=0)
        for name, value in values.items():
            homologous[name] = max(homologous.get(name, 0.), float(value))
        growth_difference = max(growth_difference, abs(growth(x, cfg)['max_growth'] - growth(y, cfg)['max_growth']))
    return maxima, homologous, peaks, growth_difference


def main():
    started = time.perf_counter()
    output = HERE / 'independent_review.json'
    assert not output.exists(), 'Never overwrite a completed review'
    for file in (Path(__file__), HELPER, prior.HELPER, t.HELPER):
        r.pin(file)
    p = r.read(ROOT / 'protocol.json')
    ph = r.pin(ROOT / 'protocol.json', EXPECTED)
    status, summary = r.read(ROOT / 'status.json'), r.read(ROOT / 'summary.json')
    assert status['state'] == 'completed_provisional'
    assert status['completed'] == status['total'] == 24
    assert status['pairs_completed'] == status['pairs_total'] == 12
    assert summary['execution_completed'] and summary['formation_screens_passed']
    assert not summary['qualification_accepted'] and not summary['scientific_accepted']
    assert summary['independent_review_pending'] and not summary['maintenance_extension_started']
    assert summary['protocol_sha256'] == status['protocol_sha256'] == ph
    assert summary['new_scientific_histories'] == 0 and summary['new_numerical_paths'] == 24
    print('Checking 63 source and 523 input pins', flush=True)
    for name, expected in {**p['source_sha256'], **p['input_sha256']}.items():
        r.pin(name, expected)
    parent = r.read(PARENT / 'protocol.json')
    r.pin(PARENT / 'protocol.json', p['parent_protocol_sha256'])
    for key in ('histories', 'branches', 'levels', 'time_pairs', 'physical_limits', 'time_limits',
                'homologous_limits', 'dense_time_spectrum_limit', 'prefix_limits', 'prefix_steps',
                'method', 'gpu_uuid', 'runtime_binary_sha256', 'interval', 'diagnostic_interval'):
        assert p[key] == parent[key], key
    for mapping in ('source_sha256', 'input_sha256'):
        assert all(p[mapping][name] == value for name, value in parent[mapping].items())
    accepted = r.read('docs/fresh_zygote_mature_clock_fix_assessment.json')
    assert accepted['scoped_mature_timestep_accepted']
    preceding = r.read(accepted['independent_review']['path'])
    r.pin(accepted['independent_review']['path'], accepted['independent_review']['sha256'])
    assert preceding['scoped_mature_timestep_accepted']
    assert r.read('output/studies/fresh-zygote-mature-qualification/status.json')['state'] == 'failed'
    assert r.read('output/studies/fresh-zygote-pilot/status.json')['state'] == 'failed'
    assert not r.read('output/studies/fresh-zygote-pilot/time-h1103-off.json')['passed']
    assert p['histories'] == [1101, 1102, 1103] and p['branches'] == list(prior.BRANCHES)
    assert p['start_time'] == 18. and p['duration'] == 60.
    assert p['snapshots'] == [0, 12, 24, 36, 48, 60] and p['restart_elapsed'] == [12, 36]
    assert p['interval'] == .15 and p['diagnostic_interval'] == .015 and p['restart_steps'] == 8
    assert p['thresholds'] == [.05, .1, .2] and p['sustained_window'] == p['late_window'] == 24.
    assert p['onset_time_limit'] == .3 and p['growth_abs_limit'] == .001
    runtime = r.read(ROOT / 'runtime.json')
    r.pin(runtime['binary'], p['runtime_binary_sha256'])
    assert runtime['binary_sha256'] == p['runtime_binary_sha256']
    assert len(p['jobs']) == len(parent['jobs']) == 24
    jobs = {j['key']: j for j in p['jobs']}
    for job, old in zip(p['jobs'], parent['jobs']):
        assert job['config'] == dict(old['config'], steps=round(78 / old['config']['dt']))
        assert job['source_case'] == old['key']
        assert Path(job['source_checkpoint']).resolve() == Path(old['source']).resolve()
        assert Path(job['source']).resolve() == (ROOT / job['key'] / 'source.npz').resolve()
        assert tuple(job['config'][k] for k in ('feedback', 'fate_tension', 'fate_adhesion', 'polarity_tension')) == prior.COUPLINGS[job['branch']]
    paths, prefixes, pairs, samples, endpoints = [], [], [], [], []
    checked_start_arrays = 0
    for pair in p['time_pairs']:
        histories, networks, reports = {}, {}, {}
        a, b = pair['cases']
        for key in (a, b):
            job, cfg, folder = jobs[key], jobs[key]['config'], ROOT / key
            assert cfg['dt'] == prior.LEVELS[job['level']] and cfg['grid'] == 120 and cfg['extent'] == 2.24
            report = r.read(folder / 'gpu-result.json')
            assert report['passed'] and report['failure'] is None and report['job'] == job
            assert report['protocol_sha256'] == ph and report['final_time'] == 78. and report['final_cells'] == 16
            assert len(report['events']) == 15 and report['lineage_cycle_rng_unchanged']
            assert report['exact_restart_steps'] == {'12': 8, '36': 8}
            audit, limits = report['audit'], p['physical_limits']
            assert audit['max_volume_error'] < limits['volume_error_max'] and audit['min_radius'] >= limits['radius_min']
            assert audit['max_clipping'] == 0 and audit['dilution_amount_error'] <= limits['dilution_amount_error_max']
            for name, expected in report['evidence_sha256'].items():
                r.pin(name, expected)
            history = r.read(folder / 'gpu-history.json')
            assert len(history) == 401
            closure, boundary = 0., 0.
            for index, row in enumerate(history):
                r.close(row['time'], 18. + index * .15, 'regular time')
                r.close(row['elapsed'], index * .15, 'regular elapsed')
                assert row['ids'] == history[0]['ids'] and row['parents'] == history[0]['parents']
                assert len(row['ids']) == 16 and not row['dividing_ids']
                assert row['nondividing_components'] == [1] * 16 and row['boundary_occupancy'] < limits['boundary_max']
                u, v = np.asarray(row['chemistry']), np.asarray(row['volumes'])
                assert u.shape == (2, 16) and np.isfinite(u).all() and u.min() > 0
                assert np.isfinite(v).all() and v.min() > 0 and np.isfinite(row['polarity']).all()
                r.close(row['log_activator_sd'], np.log(u[0]).std(), 'observed contrast')
                calculated = growth(row, cfg)
                r.close(calculated['max_growth'], row['frozen_uniform_max_growth'], 'uniform growth')
                assert calculated['unstable_modes'] == row['frozen_uniform_unstable_modes']
                closure = max(closure, r.eig(row)[1])
                boundary = max(boundary, row['boundary_occupancy'])
            calculated_outcomes = outcomes(history, p)
            same_outcomes(calculated_outcomes, report['formation'])
            trace = t.jsonl(folder / 'network.jsonl')
            assert len(trace) == report['network_trace_rows'] == 4001
            for index, row in enumerate(trace):
                r.close(row['time'], 18. + index * .015, 'dense time')
                r.close(row['time'], row['step'] * cfg['dt'], 'dense clock')
                assert row['ids'] == history[0]['ids'] and row['parents'] == history[0]['parents']
                assert not row['active_divisions'] and row['completed_events'] == 15
                assert row['last_step_projection_max'] == 0
                closure = max(closure, t.graph(row, cfg)[1])
                if index % 10 == 0:
                    regular = history[index // 10]
                    np.testing.assert_allclose(row['delta'], regular['delta'], rtol=1e-13, atol=1e-14)
                    np.testing.assert_allclose(row['volumes'], regular['volumes'], rtol=1e-13, atol=1e-14)
            prefix = r.read(folder / 'actual-start-prefix.json')
            assert prefix['passed'] and prefix['failure'] is None and prefix['steps'] == 8 and prefix['exact_gpu_restart']
            assert prefix['protocol_sha256'] == ph
            assert prefix == next(x for x in summary['prefixes'] if x['case'] == key)
            for name, value in prefix['errors'].items():
                assert value <= p['prefix_limits'][name]
            for name, expected in prefix['evidence_sha256'].items():
                r.pin(name, expected)
            prefixes.append(dict(case=key, recorded_passed=True, errors=prefix['errors'],
                                 saved_evidence_hashed=True, independently_rerun=False))
            histories[key], networks[key], reports[key] = history, trace, report
            paths.append(dict(case=key, audit=audit, boundary_max=boundary, graph_closure_max=closure,
                              regular_rows=401, dense_rows=4001, recorded_exact_restart_steps={'12': 8, '36': 8},
                              independently_recomputed_formation=calculated_outcomes))
            print('Independently checked regular chemistry, growth and all dense graphs:', key, flush=True)
        saved = r.read(ROOT / f"time-h{pair['seed']}-{pair['branch']}.json")
        assert saved == next(x for x in summary['pairs'] if x['cases'] == [a, b])
        assert saved['passed'] and not saved['failed_gates'] and saved['protocol_sha256'] == ph
        for name, expected in saved['evidence_sha256'].items():
            r.pin(name, expected)
        errors, homologous, peaks, grow = compare(histories[a], histories[b], jobs[b]['config'])
        for name, value in errors.items():
            r.close(value, saved['errors'][name], (a, name))
            assert value <= p['time_limits'][name]
        for name, value in homologous.items():
            r.close(value, saved['homologous_errors'][name], (a, name))
            assert value <= p['homologous_limits'][name]
        r.close(grow, saved['growth_abs_max'], 'paired growth maximum')
        assert grow <= p['growth_abs_limit']
        formation = [outcomes(histories[key], p) for key in (a, b)]
        for x, y in zip(formation, saved['formation']):
            same_outcomes(x, y)
        formation_checks = []
        for x, y in zip(*formation):
            agree = x['sustained_formation'] == y['sustained_formation'] and x['late_persistent'] == y['late_persistent']
            onset_error = abs(x['first_sustained_onset_elapsed'] - y['first_sustained_onset_elapsed']) if x['sustained_formation'] and y['sustained_formation'] else None
            check = dict(threshold=x['threshold'], classifications_agree=agree, sustained_onset_time_error=onset_error,
                         passed=agree and (onset_error is None or onset_error <= p['onset_time_limit'] + 1e-9))
            assert check['passed']
            formation_checks.append(check)
        assert formation_checks == saved['formation_checks']
        dense, topology = [], []
        for x, y in zip(networks[a], networks[b]):
            assert abs(x['time'] - y['time']) < 1e-9
            sx, sy = r.eig(x)[0], r.eig(y)[0]
            dense.append(dict(time=x['time'], relative_spectrum=float(np.linalg.norm(sx - sy) / np.linalg.norm(sy)),
                              edges=[x['edges'], y['edges']]))
            if not np.array_equal(np.array(x['conductance']) > 0, np.array(y['conductance']) > 0):
                topology.append(x['time'])
        peak = max(dense, key=lambda x: x['relative_spectrum'])
        r.close(peak['relative_spectrum'], saved['dense_peak']['relative_spectrum'], 'dense maximum')
        assert peak['time'] == saved['dense_peak']['time'] and peak['edges'] == saved['dense_peak']['edges']
        assert topology == saved['topology_mismatch_times'] and peak['relative_spectrum'] <= p['dense_time_spectrum_limit']
        fields, constants = [], {}
        for elapsed in p['snapshots']:
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
                assert check['volume_error'] < p['physical_limits']['volume_error_max']
                assert check['minimum_radius'] >= p['physical_limits']['radius_min']
                if elapsed == 0:
                    source = r.load(jobs[key]['source'])
                    assert set(source) == set(z)
                    for name in z:
                        assert z[name].dtype == source[name].dtype and z[name].shape == source[name].shape
                        assert z[name].tobytes() == source[name].tobytes(), (key, 'initial inheritance', name)
                    old = r.load(jobs[key]['source_checkpoint'])
                    assert set(old) == set(source)
                    ma, mb = [json.loads(str(q['metadata'])) for q in (old, source)]
                    assert mb['config'] == dict(ma['config'], steps=round(78 / cfg['dt']))
                    ma.pop('config'); mb.pop('config'); assert ma == mb
                    for name in set(source) - {'metadata'}:
                        assert source[name].dtype == old[name].dtype and source[name].shape == old[name].shape
                        assert source[name].tobytes() == old[name].tobytes(), (key, 'retimed inheritance', name)
                        checked_start_arrays += 1
                    del source, old
                    constants[key] = prior.invariant(z)
                else:
                    assert prior.invariant(z) == constants[key], (key, 'lineage/target/RNG drift')
                if elapsed in p['restart_elapsed']:
                    assert r.pin(ROOT / key / f'restart-e{elapsed:02d}.npz') == r.pin(ROOT / key / f'gpu-e{elapsed:02d}.npz')
                densities[key] = r.density(z)
                edge = densities[key][0]
                boundary = min(1., max(float(face.max()) for face in
                    (edge[0], edge[-1], edge[:, 0], edge[:, -1], edge[:, :, 0], edge[:, :, -1])))
                r.close(boundary, observed['boundary_occupancy'], ('boundary', key, elapsed))
                assert boundary < p['physical_limits']['boundary_max']
                samples.append(dict(case=key, **check, boundary_occupancy=boundary, inherited_constants_exact=True))
                del z
            values = r.field_error(densities[a], densities[b])
            expected = next(x for x in saved['fields'] if x['elapsed'] == elapsed)
            for name, value in values.items():
                r.close(value, expected['errors'][name], ('density', a, elapsed, name))
                assert value <= p['time_limits']['field_max']
            fields.append(dict(elapsed=elapsed, time=18. + elapsed, errors=values))
            del densities
            print('Independently reconstructed paired checkpoint states and native fields:', a, elapsed, flush=True)
        field_max = max(v for x in fields for v in x['errors'].values())
        r.close(field_max, saved['field_max'], 'maximum native field')
        pairs.append(dict(history=pair['seed'], branch=pair['branch'], cases=[a, b], independently_passed=True,
            errors=errors, homologous_errors=homologous, peak_times=peaks, dense_peak=peak,
            topology_mismatch_times=topology, growth_abs_max=grow, fields=fields, field_max=field_max,
            formation=formation, formation_checks=formation_checks))
        first, last = [prior.endpoint(histories[b][i], networks[b][i * 10], jobs[b]['config']) for i in (0, 400)]
        endpoints.append(dict(history=pair['seed'], branch=pair['branch'], dt=prior.LEVELS['fine'],
            start=first, end=last, formation=formation[1],
            log_contrast_growth=last['log_activator_sd'] / first['log_activator_sd'],
            conductance_relative_change=last['summed_undirected_conductance'] / first['summed_undirected_conductance'] - 1,
            frozen_uniform_max_growth_start=growth(histories[b][0], jobs[b]['config'])['max_growth'],
            frozen_uniform_max_growth_end=growth(histories[b][-1], jobs[b]['config'])['max_growth']))
        del histories, networks, reports
    assert len(samples) == 144 and len(pairs) == len(endpoints) == 12 and checked_start_arrays == 432
    transfer = r.read(ROOT / 'transfer-verification.json')
    assert transfer['passed'] and transfer['checked_arrays'] == 432 and transfer['all24_starts_lossless']
    print('Rechecking every pinned input and production record for changes', flush=True)
    for name, expected in list(r.WATCHED.items()):
        assert r.sha(name) == expected, name
    result = dict(review_completed=True, scoped_formation_timestep_accepted=True,
        reviewed_utc=datetime.now(timezone.utc).isoformat(), seconds=time.perf_counter() - started,
        protocol_sha256=ph, source_count=len(p['source_sha256']), input_count=len(p['input_sha256']),
        histories=p['histories'], independent_histories=3, new_scientific_histories=0, new_numerical_paths=24,
        paths=paths, prefixes=prefixes, pairs=pairs, checkpoints=samples, fine_endpoint_observations=endpoints,
        initial_arrays_losslessly_checked=checked_start_arrays,
        original_scheduling_failure_retained=True, original_developmental_failure_retained=True,
        full_cpu_horizon_rerun=False, exact_restart_independently_rerun=False,
        maintenance_240_qualified=False, biological_identity_qualified=False,
        general_developmental_convergence_accepted=False, maintenance_extension_started=False,
        scope='Independent NumPy/SciPy saved-array audit: 24 numerical paths within three histories, '
              '9624 regular observations, 96024 dense conservative graphs, all twelve complete timestep '
              'comparisons, 144 complete checkpoint states/native density reconstructions, formation '
              'at thresholds0.05/0.1/0.2 with full24-unit windows, frozen chemical growth, exact inherited '
              'starts and unchanged lineage/targets/RNG. All24 prefixes and eight-step restarts checked '
              'through immutable saved evidence, not independently rerun. Both restart-start archives '
              'are byte-identical to their corresponding full checkpoints. Acceptance applies only '
              'to dt0.00375/0.001875, grid120, physical t18–78, twelve inherited contexts. '
              'No240-unit maintenance, biological identity or universal convergence claim.',
        original_files_sha256=dict(r.WATCHED))
    output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ('review_completed', 'scoped_formation_timestep_accepted', 'seconds')}, indent=2), flush=True)


if __name__ == '__main__':
    alias = module('verified_study_alias', Path('scripts/with_study_output_path_alias.py'))
    with alias.base.validated_alias(Path('archive/output-relocation/2026-10-10-output/studies-manifest.json'),
            '8fb2545d8bae35a1a673ec3a43d3cc191148f03cac58892f29bd368182a8aff6',
            source_root=alias.SOURCE_ROOT, target_root=alias.TARGET_ROOT):
        main()
