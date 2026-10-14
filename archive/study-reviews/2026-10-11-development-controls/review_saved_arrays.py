"""Independent NumPy/SciPy review of completed matched developmental controls.

Reuse the separately written, accepted saved-array review algorithms, never
production comparisons or GPU stepping. Preserve original reports and limits.
"""
from datetime import datetime, timezone
import importlib.util
import json
from pathlib import Path
import time

import numpy as np

ROOT = Path('output/studies/development-control-qualification')
HERE = Path(__file__).parent
EXPECTED_PROTOCOL_SHA256 = 'c57f1e21732680319147c1c1bff11f2f8aea2a48753d32fdab4f37bdc0a66f88'
HELPER = Path('archive/study-reviews/2026-10-10-development-space-refinement/review_saved_arrays.py')
spec = importlib.util.spec_from_file_location('independent_space_review', HELPER)
s = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s)
t, r = s.t, s.r
read, load, pin, close = r.read, r.load, r.pin, r.close


def review_path(job, protocol, protocol_hash):
    key = job['key']; folder = ROOT / key; cfg = job['config']
    report = read(folder / 'gpu-result.json')
    assert report['passed'] and report['failure'] is None
    if job['reused']:
        old = protocol['reuse_provenance'][key]
        assert folder.is_symlink() and folder.resolve() == Path(old['folder'])
        assert pin(folder / 'gpu-result.json') == old['report_sha256']
        assert report['job'] == old['job'] and report['protocol_sha256'] == old['protocol_sha256']
        assert report['job']['config'] == cfg
    else:
        assert not folder.is_symlink() and report['job'] == job
        assert report['protocol_sha256'] == protocol_hash
    for file, expected in report['evidence_sha256'].items(): pin(file, expected)
    history = read(folder / 'gpu-history.json')
    events = read(folder / 'gpu-events.json'); lineage = read(folder / 'gpu-lineage.json')
    assert events == report['events']
    assert len(history) == 121 and len(events) == 15 and len(lineage) == 31
    assert report['final_time'] == 18 and report['final_cells'] == 16
    limits = protocol['physical_limits']; audit = report['audit']
    assert audit['max_volume_error'] < limits['volume_error_max']
    assert audit['min_radius'] >= limits['radius_min'] and audit['max_clipping'] == 0
    assert audit['dilution_amount_error'] <= limits['dilution_amount_error_max']
    assert report['exact_restart_steps'] == {'2': 8, '4': 8, '8': 8}
    assert report['independent_chemical_steps'] == 8 and report['independent_chemical_log_max'] <= 1e-11
    assert 18 - max(e['time'] for e in events) >= protocol['minimum_mature_window']
    assert len({e['parent'] for e in events}) == 15
    for e in events:
        assert e['neck'] <= .15 and e['overlap'] <= .002 and .1 < e['fraction'] < .9
        assert e['daughter_components'] == [1, 1]
        for k, limit in [('amount_relative_jump', 'event_amount_jump_max'),
                         ('target_relative_jump', 'event_target_jump_max'),
                         ('volume_relative_jump', 'event_volume_jump_max'),
                         ('effective_occupancy_absolute_jump', 'event_occupancy_jump_max')]:
            assert e[k] <= limits[limit]
        for child in e['children']: assert lineage[child]['parent'] == e['parent']
        np.testing.assert_allclose(np.array(e['chemistry_before']) @ np.array(e['volumes_before']),
                                   np.array(e['chemistry_after']) @ np.array(e['volumes_after']),
                                   rtol=2e-14, atol=1e-15)
    closure = 0.; boundary = 0.
    for i, row in enumerate(history):
        assert abs(row['time'] - i * .15) < 1e-9
        n = len(row['ids']); u = np.array(row['chemistry']); v = np.array(row['volumes'])
        assert u.shape == (2, n) and np.isfinite(u).all() and u.min() > 0
        assert np.isfinite(v).all() and v.min() > 0
        assert all(c == 1 for c in row['nondividing_components'])
        assert len(row['nondividing_components']) == n - len(row['dividing_ids'])
        assert row['boundary_occupancy'] < limits['boundary_max']
        _, error = r.eig(row); closure = max(closure, error)
        boundary = max(boundary, row['boundary_occupancy'])
    trace = t.jsonl(folder / 'network.jsonl'); assert len(trace) == report['network_trace_rows']
    previous = -1.; regular = {}
    for row in trace:
        assert row['time'] > previous; previous = row['time']
        close(row['time'], row['step'] * cfg['dt'], 'network clock')
        _, error = t.graph(row, cfg); closure = max(closure, error)
        if abs(row['time'] / .015 - round(row['time'] / .015)) < 1e-7:
            regular[round(row['time'], 9)] = row
    assert len(regular) == 1201
    pressure = t.jsonl(folder / 'projection.jsonl'); seen = set()
    parents = {e['parent'] for e in events}
    for row in pressure:
        assert row['parent'] in parents and (row['step'], row['parent']) not in seen
        seen.add((row['step'], row['parent']))
        close(row['time_start'], row['step'] * cfg['dt'], 'projection clock')
        close(row['time_end'] - row['time_start'], cfg['dt'], 'projection dt')
        assert row['capacity'] > 0 and np.isfinite(list(row.values())).all()
        assert abs(row['after_relative_volume_error']) <= 1e-7 + 1e-14
        assert row['correction_magnitude'] >= 0
    extrema = dict(calls=len(pressure),
                   correction_magnitude_max=max(x['correction_magnitude'] for x in pressure),
                   before_relative_volume_error_max=max(abs(x['before_relative_volume_error']) for x in pressure),
                   after_relative_volume_error_max=max(abs(x['after_relative_volume_error']) for x in pressure))
    for k, value in extrema.items(): close(value, report['projection'][k], ('projection', key, k))
    result = dict(case=key, reused=job['reused'], seconds=report['seconds'], audit=audit,
                  boundary_max=boundary, graph_closure_max=closure, projection=extrema,
                  network_rows=len(trace), recorded_runtime_event_gates_passed=True,
                  log_activator_sd_at_18=history[-1]['log_activator_sd'])
    print('Reviewed history, dense graph, events and pressure:', key, flush=True)
    return result, report, history, regular


def review_fields(branch, original, protocol, histories, reports, networks):
    spatial = [s.compare_histories(row, histories, reports, networks) for row in original['space_rows']]
    temporal = s.compare_histories(original['time_row'], histories, reports, networks)
    for actual, saved in zip([row['dense'] for row in spatial] + [temporal['dense']],
                             original['dense_space_rows'] + [original['dense_time']]):
        close(actual['peak']['relative_spectrum'], saved['peak']['relative_spectrum'], 'dense spectrum')
        assert actual['peak']['time'] == saved['peak']['time']
        assert actual['unmatched_times'] == saved['unmatched_times']
        assert actual['topology_mismatch_times'] == saved['topology_mismatch_times']
    jobs = [j for j in protocol['jobs'] if j.get('branch') == branch]
    samples = []; live = []
    for stamp in (0, 3, 6, 9, 12, 15, 18):
        primary = {}; native = {}; ranges = {}
        for job in jobs:
            key = job['key']; z = load(ROOT / key / f'gpu-t{stamp:02d}.npz')
            check = r.checkpoint(z, job['config'], histories[key][int(round(stamp / .15))], stamp == 18)
            assert check['volume_error'] < protocol['physical_limits']['volume_error_max']
            assert check['minimum_radius'] >= protocol['physical_limits']['radius_min']
            assert json.loads(str(z['event_audits'])) == [e for e in reports[key]['events'] if e['time'] <= stamp + 1e-12]
            samples.append(dict(case=key, **check))
            if job['level'] != 'time120_half': primary[key], ranges[key] = s.density(z, 112, 3)
            if key in protocol['time_pairs'][branch]: native[key], _ = s.density(z, None)
            del z
        for i, row in enumerate(spatial):
            a, b = row['cases']; values = r.field_error(primary[a], primary[b])
            expected = next(f for f in original['space_rows'][i]['fields']
                            if f['time'] == stamp and f['probe'] == 112 and f['order'] == 3)
            for k, value in values.items(): close(value, expected['errors'][k], ('space primary', branch, stamp, k))
            np.testing.assert_allclose([ranges[a], ranges[b]], expected['phase_ranges'], rtol=1e-6, atol=1e-9)
            row['field_checks'].append(dict(time=stamp, probe=112, order=3, errors=values,
                                           phase_ranges=[ranges[a], ranges[b]]))
        a, b = temporal['cases']; values = r.field_error(native[a], native[b])
        expected = next(f for f in original['time_row']['fields'] if f['time'] == stamp)
        for k, value in values.items(): close(value, expected['errors'][k], ('native time', branch, stamp, k))
        temporal['field_checks'].append(dict(time=stamp, probe=None, order=3, errors=values))
        del native, primary
        if stamp == 18:
            for probe, order in ((112, 5), (128, 3), (128, 5)):
                fields = {}; ranges = {}
                for job in jobs:
                    if job['level'] == 'time120_half': continue
                    key = job['key']; z = load(ROOT / key / 'gpu-t18.npz')
                    fields[key], ranges[key] = s.density(z, probe, order); del z
                for i, row in enumerate(spatial):
                    a, b = row['cases']; values = r.field_error(fields[a], fields[b])
                    expected = next(f for f in original['space_rows'][i]['fields']
                                    if f['time'] == 18 and f['probe'] == probe and f['order'] == order)
                    for k, value in values.items(): close(value, expected['errors'][k], ('space variants', branch, probe, order, k))
                    np.testing.assert_allclose([ranges[a], ranges[b]], expected['phase_ranges'], rtol=1e-6, atol=1e-9)
                    row['field_checks'].append(dict(time=18, probe=probe, order=order, errors=values,
                                                   phase_ranges=[ranges[a], ranges[b]]))
                del fields
        print('Rebuilt independent spatial/native fields:', branch, stamp, flush=True)
    for job in jobs:
        for stage in (2, 4, 8):
            z = load(ROOT / job['key'] / f'restart-stage-{stage}.npz')
            check = r.checkpoint(z, job['config'])
            assert stage < check['cells'] < 2 * stage and check['active'] > 0
            live.append(dict(case=job['key'], stage=stage, **check)); del z
    for row, saved in zip(spatial, original['space_rows']):
        s.decision(row, saved, protocol['space_limits'], protocol['diagnostic_spread_limit'])
    s.decision(temporal, original['time_row'], protocol['time_limits'])
    sequences = {k: [row['errors'][k] for row in spatial]
                 for k in ('concentration_wasserstein_max', 'log_wasserstein_max', 'relative_axis_max')}
    sequences.update({k: [max(f['errors'][k] for f in row['field_checks']) for row in spatial] for k in r.FIELDS})
    trend = all(all(b <= a + protocol['trend_floor'] for a, b in zip(values, values[1:]))
                for values in sequences.values())
    assert trend == original['spatial_error_decay_passed']
    dense_pass = temporal['dense']['peak']['relative_spectrum'] <= protocol['dense_time_spectrum_limit']
    passed = bool(spatial[-1]['independently_passed'] and temporal['independently_passed'] and dense_pass and trend)
    assert passed == original['scoped_screens_passed']
    return dict(branch=branch, independently_passed=passed, space_rows=spatial, time_row=temporal,
                dense_time_passed=dense_pass, spatial_error_decay_passed=trend,
                spatial_error_sequences=sequences, regular_checkpoints=samples, live_checkpoints=live)


def review_replay(key, stage, protocol, protocol_hash, histories, configs):
    folder = ROOT / 'replays' / key / f'stage-{stage}'
    saved = read(folder / 'result.json')
    assert saved['case'] == key and saved['stage'] == stage and saved['protocol_sha256'] == protocol_hash
    assert saved['passed'] and saved['failure'] is None and saved['complete_state_bitwise_equal']
    assert saved['end'] == 18 and saved['final_cells'] == 16 and saved['final_events'] == 15
    for name, expected in {**saved['input_sha256'], **saved['output_sha256']}.items(): pin(name, expected)
    source = load(ROOT / key / f'restart-stage-{stage}.npz')
    check = r.checkpoint(source, configs[key]); meta = json.loads(str(source['metadata']))
    assert stage < check['cells'] < 2 * stage and check['active'] > 0
    close(saved['start'], meta['time'], 'replay start')
    assert saved['steps'] == configs[key]['steps'] - meta['step_number']
    expected_times = [row['time'] for row in histories[key] if row['time'] >= saved['start'] - 1e-12]
    assert saved['exact_observation_times'] == expected_times
    expected_checkpoints = [float(stamp) for stamp in (0, 3, 6, 9, 12, 15, 18) if stamp >= saved['start'] - 1e-12]
    assert saved['checkpoint_times'] == expected_checkpoints
    del source
    actual = load(folder / 'final.npz'); expected = load(ROOT / key / 'gpu-t18.npz')
    assert actual.keys() == expected.keys()
    for name in actual:
        assert actual[name].dtype == expected[name].dtype and actual[name].shape == expected[name].shape
        assert actual[name].tobytes() == expected[name].tobytes(), (key, stage, name)
    final = r.checkpoint(actual, configs[key], histories[key][-1], True)
    limits = protocol['physical_limits']; audit = saved['audit']
    assert audit['max_volume_error'] < limits['volume_error_max'] and audit['min_radius'] >= limits['radius_min']
    assert audit['max_clipping'] == 0 and audit['dilution_amount_error'] <= limits['dilution_amount_error_max']
    print('Verified replay evidence and bitwise final arrays:', key, stage, flush=True)
    return dict(case=key, stage=stage, independently_checked_final_arrays=True,
                recorded_exact_observations=len(expected_times), recorded_exact_checkpoints=expected_checkpoints,
                recorded_full_replay_passed=True, seconds=saved['seconds'], final_checkpoint=final,
                note='Intermediate exact comparisons are recorded production checks, not a new independent mechanics replay.')


def main():
    target = HERE / 'independent_review.json'
    if target.exists(): raise FileExistsError(target)
    started = time.perf_counter()
    for file in (__file__, HELPER, s.HELPER, t.HELPER): pin(file)
    p = read(ROOT / 'protocol.json'); ph = pin(ROOT / 'protocol.json')
    assert ph == EXPECTED_PROTOCOL_SHA256
    summary = read(ROOT / 'summary.json'); status = read(ROOT / 'status.json'); timing = read(ROOT / 'timing.json')
    assert status['state'] in ('completed_provisional', 'completed_with_unresolved_checks')
    assert summary['execution_completed'] and not summary['qualification_accepted']
    assert summary['protocol_sha256'] == ph and summary['new_scientific_histories'] == 0
    assert len(summary['paths']) == 13 and len(summary['branch_rows']) == 3 and len(summary['replay_rows']) == 12
    for file, expected in {**p['source_sha256'], **p['input_sha256'], **summary['evidence_sha256']}.items(): pin(file, expected)
    runtime = read(ROOT / 'runtime.json'); pin(runtime['binary'], runtime['binary_sha256'])
    assert runtime['binary_sha256'] == p['runtime_binary_sha256']
    calibration = read(p['calibration_reference']); cm = read(p['calibration_manifest'])
    assert calibration['passed'] and len(calibration['rows']) == 48
    assert calibration['protocol_sha256'] == cm['protocol_sha256'] == p['parent_protocol_sha256']
    assert cm['calibration_sha256'] == pin(p['calibration_reference'])
    assert max(row['relative_l2'] for row in calibration['rows']) <= p['calibration_limit']
    gate = read(p['instrumentation_gate'])
    assert gate['passed'] and gate['protocol_sha256'] == p['instrumentation_protocol_sha256']
    for file, expected in gate['evidence_sha256'].items(): pin(file, expected)
    histories = {}; reports = {}; networks = {}; paths = []
    configs = {job['key']: job['config'] for job in p['jobs']}
    for job in p['jobs']:
        checked, report, history, network = review_path(job, p, ph)
        key = job['key']; paths.append(checked); reports[key] = report; histories[key] = history; networks[key] = network
    branches = []
    for branch in p['branches']:
        original = read(ROOT / (branch + '-refinement.json'))
        assert original == next(row for row in summary['branch_rows'] if row['branch'] == branch)
        branches.append(review_fields(branch, original, p, histories, reports, networks))
    reused_samples = []; reused_live = []
    job = next(job for job in p['jobs'] if job['key'] == 'off120')
    for stamp in (0, 3, 6, 9, 12, 15, 18):
        z = load(ROOT / 'off120' / f'gpu-t{stamp:02d}.npz')
        reused_samples.append(r.checkpoint(z, job['config'], histories['off120'][int(round(stamp / .15))], stamp == 18)); del z
    for stage in (2, 4, 8):
        z = load(ROOT / 'off120' / f'restart-stage-{stage}.npz'); checked = r.checkpoint(z, job['config'])
        assert stage < checked['cells'] < 2 * stage and checked['active'] > 0
        reused_live.append(dict(stage=stage, **checked)); del z
    required = {(key, stage) for key in p['replay_cases'] for stage in p['restart_stages']}
    assert required == {(row['case'], row['stage']) for row in summary['replay_rows']}
    replays = []
    for key in p['replay_cases']:
        for stage in p['restart_stages']:
            assert read(ROOT / 'replays' / key / f'stage-{stage}' / 'result.json') == next(
                row for row in summary['replay_rows'] if row['case'] == key and row['stage'] == stage)
            replays.append(review_replay(key, stage, p, ph, histories, configs))
    preserved = read(p['preserved_manifest'])
    for file, expected in preserved['preserved_sha256'].items(): pin(file, expected)
    passed = all(row['independently_passed'] for row in branches) and len(replays) == 12
    assert passed == summary['scoped_screens_passed'] == status['scoped_screens_passed']
    # Re-read hashes without the pin cache to detect any concurrent mutation.
    for file, expected in r.WATCHED.items(): assert r.sha(file) == expected, file
    result = dict(review_completed=True, scoped_qualification_accepted=bool(passed),
                  qualification_accepted=bool(passed), biological_identity_qualified=False,
                  general_continuum_qualified=False, fine_cpu_ladder_performed=False,
                  independently_repeated_mechanics=False, protocol_sha256=ph,
                  reviewed_at_utc=datetime.now(timezone.utc).isoformat(), source_count=len(p['source_sha256']),
                  input_count=len(p['input_sha256']), preserved_count=len(preserved['preserved_sha256']),
                  paths=paths, branch_rows=branches, replay_rows=replays,
                  reused_regular_checkpoints=reused_samples, reused_live_checkpoints=reused_live,
                  recorded_checks_not_independently_replayed=[
                      'Every-step extrema, eight-step prefixes and intermediate full-replay exact checks; logged evidence and saved arrays inspected.',
                      'Inherited reconstruction calibration and instrumentation; original records and hashes verified.'],
                  consumed_sha256=dict(r.WATCHED), new_scientific_histories=0,
                  new_numerical_paths=12, reused_numerical_paths=1,
                  resumed_coordinator_seconds=timing['seconds'], review_seconds=time.perf_counter()-started,
                  scope=p['scope'], original_reports_unchanged=True)
    with target.open('x') as stream: stream.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ('review_completed', 'scoped_qualification_accepted',
                                          'preserved_count', 'review_seconds')}), flush=True)


if __name__ == '__main__':
    alias_path = Path('scripts/with_study_output_path_alias.py')
    alias_spec = importlib.util.spec_from_file_location('study_review_alias', alias_path)
    alias = importlib.util.module_from_spec(alias_spec); alias_spec.loader.exec_module(alias)
    record = Path('archive/output-relocation/2026-10-10-output/studies-manifest.json')
    with alias.base.validated_alias(record,
            '8fb2545d8bae35a1a673ec3a43d3cc191148f03cac58892f29bd368182a8aff6',
            source_root=alias.SOURCE_ROOT, target_root=alias.TARGET_ROOT):
        for file in (alias_path, 'scripts/with_output_path_alias.py', record): pin(file)
        main()
