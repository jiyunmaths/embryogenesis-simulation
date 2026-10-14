"""Saved-evidence review; NumPy/SciPy only, no embryo imports or GPU stepping.

Use the separately written preceding review's independent CDF, generalized
eigenproblem, checkpoint and density reconstructions. Do not call production
comparison helpers. Preserve all production reports and earlier failures.
"""
from pathlib import Path
from datetime import datetime, timezone
import importlib.util
import json
import time
import numpy as np

ROOT = Path('output/studies/development-time-diagnosis')
HERE = Path(__file__).parent
HELPER = Path('archive/study-reviews/2026-10-10-repeated-development-gpu/review_saved_arrays.py')
spec = importlib.util.spec_from_file_location('independent_prior_review', HELPER)
r = importlib.util.module_from_spec(spec)
spec.loader.exec_module(r)
read, load, pin, close = r.read, r.load, r.pin, r.close


def jsonl(path):
    pin(path)
    return [json.loads(line) for line in Path(path).open()]


def graph(row, config):
    values, error = r.eig(row)
    v = np.array(row['volumes']); d = np.array(row['delta'])
    g = np.array(row['conductance']); w = np.array(row['contact'])
    n = len(v); upper = np.triu_indices(n, 1)
    assert np.isfinite(g).all() and np.isfinite(w).all()
    np.testing.assert_allclose(g, g.T, rtol=1e-13, atol=1e-14)
    np.testing.assert_allclose(w, w.T, rtol=1e-13, atol=1e-14)
    assert g.min() >= 0 and w.min() >= 0
    assert not np.diag(g).any() and not np.diag(w).any()
    expected = g/v[:, None]
    np.fill_diagonal(expected, -expected.sum(1))
    np.testing.assert_allclose(d, expected, rtol=1e-12, atol=1e-13)
    close(row['cutoff'], config['graph_contact_cutoff']*w.max(), 'cutoff')
    np.testing.assert_allclose(row['contact_cutoff_margins'],
        w[upper]-row['cutoff'], rtol=1e-13, atol=1e-14)
    assert row['edges'] == np.count_nonzero(g[upper])
    np.testing.assert_array_equal(g[upper] > 0, w[upper] > row['cutoff'])
    for a in row['active_divisions']:
        assert a['parent'] in row['ids']
        p = np.clip((row['time']-a['start'])/config['cytokinesis_duration'], 0., 1.)
        close(a['progress'], p, 'furrow progress')
        close(a['smooth_ramp'], p*p*(3-2*p), 'ring ramp')
        close(a['neck_margin'], config['neck_threshold']-a['neck'], 'neck margin')
        close(a['overlap_margin'], config['division_overlap_tolerance']-a['overlap'], 'overlap margin')
        assert 0 < a['fraction'] < 1
    return values, error


def main():
    started = time.perf_counter(); pin(__file__); pin(HELPER)
    p = read(ROOT/'protocol.json'); ph = pin(ROOT/'protocol.json')
    summary = read(ROOT/'summary.json'); status = read(ROOT/'status.json')
    timing = read(ROOT/'timing.json'); runtime = read(ROOT/'runtime.json')
    assert status['state'] == 'completed_diagnostic' and not summary['qualification_accepted']
    assert summary['protocol_sha256'] == ph and summary['new_scientific_histories'] == 0
    for file, expected in {**p['source_sha256'], **p['input_sha256'], **summary['evidence_sha256']}.items():
        pin(file, expected)
    pin(runtime['binary'], runtime['binary_sha256'])
    gate = read(ROOT/'instrumentation-gate.json')
    assert gate['passed'] and gate['protocol_sha256'] == ph
    assert len(gate['rows']) == 3 and sum(i['steps'] for i in gate['rows']) == 24
    for file, expected in gate['evidence_sha256'].items(): pin(file, expected)
    histories = {}; reports = {}; networks = {}; lineages = {}; paths = []
    for job in p['jobs']:
        key = job['key']; folder = ROOT/key; cfg = job['config']
        report = read(folder/'gpu-result.json'); reports[key] = report
        assert report['passed'] and report['failure'] is None and report['job'] == job
        assert report['protocol_sha256'] == ph
        for file, expected in report['evidence_sha256'].items(): pin(file, expected)
        history = read(folder/'gpu-history.json'); histories[key] = history
        events = read(folder/'gpu-events.json'); assert events == report['events']
        lineage = read(folder/'gpu-lineage.json'); lineages[key] = lineage
        assert len(history) == 121 and len(events) == 15 and len(lineage) == 31
        assert report['final_time'] == 18 and report['final_cells'] == 16
        audit = report['audit']; limits = p['physical_limits']
        assert audit['max_volume_error'] < limits['volume_error_max']
        assert audit['min_radius'] >= limits['radius_min'] and audit['max_clipping'] == 0
        assert audit['dilution_amount_error'] <= limits['dilution_amount_error_max']
        assert report['exact_restart_steps'] == {'2': 8, '4': 8, '8': 8}
        assert report['independent_chemical_steps'] == 8 and report['independent_chemical_log_max'] <= 1e-11
        assert 18-max(e['time'] for e in events) >= p['minimum_mature_window']
        for e in events:
            assert e['neck'] <= .15 and e['overlap'] <= .002 and .1 < e['fraction'] < .9
            assert e['daughter_components'] == [1, 1]
            for k, l in [('amount_relative_jump','event_amount_jump_max'),
                ('target_relative_jump','event_target_jump_max'),('volume_relative_jump','event_volume_jump_max'),
                ('effective_occupancy_absolute_jump','event_occupancy_jump_max')]:
                assert e[k] <= limits[l]
            for child in e['children']: assert lineage[child]['parent'] == e['parent']
            np.testing.assert_allclose(np.array(e['chemistry_before'])@np.array(e['volumes_before']),
                np.array(e['chemistry_after'])@np.array(e['volumes_after']), rtol=2e-14, atol=1e-15)
        closure = 0.; boundary = 0.
        for i, row in enumerate(history):
            assert abs(row['time']-i*.15) < 1e-9
            n = len(row['ids']); u = np.array(row['chemistry']); v = np.array(row['volumes'])
            assert u.shape == (2, n) and np.isfinite(u).all() and u.min() > 0
            assert np.isfinite(v).all() and v.min() > 0
            assert all(c == 1 for c in row['nondividing_components'])
            assert len(row['nondividing_components']) == n-len(row['dividing_ids'])
            assert row['boundary_occupancy'] < limits['boundary_max']
            _, er = r.eig(row); closure = max(closure, er); boundary = max(boundary, row['boundary_occupancy'])
        trace = jsonl(folder/'network.jsonl'); assert len(trace) == report['network_trace_rows']
        regular = {}; previous = -1.
        for row in trace:
            assert row['time'] > previous; previous = row['time']
            close(row['time'], row['step']*cfg['dt'], 'network clock')
            _, er = graph(row, cfg); closure = max(closure, er)
            if abs(row['time']/.015-round(row['time']/.015)) < 1e-7:
                regular[round(row['time'], 9)] = row
        assert len(regular) == 1201
        networks[key] = regular
        pressure = jsonl(folder/'projection.jsonl')
        extrema = dict(calls=len(pressure), correction_magnitude_max=max(i['correction_magnitude'] for i in pressure),
            before_relative_volume_error_max=max(abs(i['before_relative_volume_error']) for i in pressure),
            after_relative_volume_error_max=max(abs(i['after_relative_volume_error']) for i in pressure))
        for name, value in extrema.items(): close(value, report['projection'][name], ('projection', name))
        seen = set(); parent_ids = {e['parent'] for e in events}
        for row in pressure:
            assert row['parent'] in parent_ids and (row['step'],row['parent']) not in seen
            seen.add((row['step'],row['parent']))
            close(row['time_start'], row['step']*cfg['dt'], 'projection clock')
            close(row['time_end']-row['time_start'], cfg['dt'], 'projection dt')
            assert row['capacity'] > 0 and np.isfinite(list(row.values())).all()
            assert abs(row['after_relative_volume_error']) <= 1e-7+1e-14
            assert row['correction_magnitude'] >= 0
        old_reproduction = None
        if job['level'] in ('coarse', 'fine'):
            old = Path(p['parent_study'])/key
            for name in ('gpu-history.json','gpu-events.json','gpu-lineage.json'):
                assert read(folder/name) == read(old/name)
            for stamp in (0,3,6,9,12,15,18):
                a = load(folder/f'gpu-t{stamp:02d}.npz'); b = load(old/f'gpu-t{stamp:02d}.npz')
                assert a.keys() == b.keys()
                for name in a: np.testing.assert_array_equal(a[name], b[name])
            old_reproduction = dict(exact=True, checkpoints=7)
        paths.append(dict(case=key, seconds=report['seconds'], recorded_runtime_gates_pass=True,
            boundary_max=boundary, graph_closure_max=closure, audit=audit, projection=extrema,
            network_rows=len(trace), regular_network_rows=len(regular), old_reproduction=old_reproduction,
            log_activator_sd_at_18=history[-1]['log_activator_sd']))
        print('Reviewed histories, every-call logs, dense spectra and old reproduction:', key, flush=True)
    pairs = []; names = [j['key'] for j in p['jobs']]
    for original in summary['pairs']:
        a,b = original['cases']; errors, peaks, unmatched, raw, raw_unmatched = r.compare(histories[a], histories[b])
        for k,v in errors.items(): close(v, original['errors'][k], ('history',a,b,k))
        assert unmatched == original['unmatched_times']
        maturity = abs(max(e['time'] for e in reports[a]['events'])-max(e['time'] for e in reports[b]['events']))
        close(maturity, original['maturity_time_error'], 'maturity')
        ea = {e['parent']:e for e in reports[a]['events']}; eb = {e['parent']:e for e in reports[b]['events']}
        assert ea.keys() == eb.keys()
        for k in ea: assert ea[k]['children'] == eb[k]['children'] and ea[k]['draws'] == eb[k]['draws']
        la,lb = lineages[a],lineages[b]
        assert [(i['id'],i['parent']) for i in la] == [(i['id'],i['parent']) for i in lb]
        event_time = max(abs(ea[k]['time']-eb[k]['time']) for k in ea)
        angles = [float(np.degrees(np.arccos(np.clip(abs(np.dot(i['division_axis'],j['division_axis'])),0,1))))
            for i,j in zip(la,lb) if 'division_axis' in i]
        dense = []; counts_unmatched = []; topology = []; conductance = []
        for t,x in networks[a].items():
            y = networks[b][t]
            if len(x['ids']) != len(y['ids']): counts_unmatched.append(t); continue
            sx, sy = r.eig(x)[0], r.eig(y)[0]
            er = float(np.linalg.norm(sx-sy)/(np.linalg.norm(sy) or 1.))
            dense.append(dict(time=t, relative_spectrum=er, cells=len(x['ids']), edges=[x['edges'],y['edges']]))
            if x['ids'] == y['ids'] and x['parents'] == y['parents']:
                ga,gb = np.array(x['conductance']),np.array(y['conductance'])
                conductance.append(dict(time=t, relative_conductance=float(np.linalg.norm(ga-gb)/(np.linalg.norm(gb) or 1.))))
                if not np.array_equal(ga>0,gb>0): topology.append(t)
        pairs.append(dict(cases=[a,b],errors=errors,peak_times=peaks,unmatched_times=unmatched,
            homologous_errors=raw, homologous_unmatched_times=raw_unmatched, maturity_time_error=maturity,
            parent_abscission_time_max=event_time, axis_angle_degrees_max=max(angles), field_checks=[],
            dense_spectrum_peak=max(dense,key=lambda x:x['relative_spectrum']),
            dense_count_unmatched_times=counts_unmatched, dense_topology_mismatch_times=topology,
            dense_conductance_peak=max(conductance,key=lambda x:x['relative_conductance'])))
    samples = []; restarts = []
    for stamp in (0,3,6,9,12,15,18):
        densities = {}
        for job in p['jobs']:
            key = job['key']; z = load(ROOT/key/f'gpu-t{stamp:02d}.npz')
            check = r.checkpoint(z,job['config'],histories[key][int(round(stamp/.15))],stamp==18)
            assert json.loads(str(z['event_audits'])) == [e for e in reports[key]['events'] if e['time'] <= stamp+1e-12]
            samples.append(dict(case=key,**check)); densities[key] = r.density(z)
        for i,(a,b) in enumerate(zip(names,names[1:])):
            errors = r.field_error(densities[a],densities[b])
            expected = next(v for v in summary['pairs'][i]['fields'] if v['time'] == stamp)
            for k,v in errors.items(): close(v,expected['errors'][k], ('field',a,b,stamp,k))
            pairs[i]['field_checks'].append(dict(time=stamp,errors=errors))
        print('Independently reconstructed four native checkpoint fields at',stamp,flush=True)
    for job in p['jobs']:
        for stage in (2,4,8):
            z = load(ROOT/job['key']/f'restart-stage-{stage}.npz'); check = r.checkpoint(z,job['config'])
            assert stage < check['cells'] < 2*stage and check['active'] > 0
            restarts.append(dict(case=job['key'],stage=stage,**check))
    for i,row in enumerate(pairs):
        fm = max(v for f in row['field_checks'] for v in f['errors'].values())
        close(fm,summary['pairs'][i]['field_max'],'field max')
        failed = {k:v for k,v in row['errors'].items() if v > p['time_limits'][k]}
        if fm > p['time_limits']['field_max']: failed['field_max'] = fm
        if row['maturity_time_error'] > p['time_limits']['maturity_time_difference_max']:
            failed['maturity_time_error'] = row['maturity_time_error']
        row.update(field_max=fm, failed_gates=failed, independently_passed=not failed)
        assert row['independently_passed'] == summary['pairs'][i]['passed']
    preserved = read(p['preserved_manifest'])
    for file,expected in preserved['preserved_sha256'].items(): pin(file,expected)
    result = dict(review_completed=True, temporal_screen_passed=pairs[-1]['independently_passed'],
        full_development_qualified=False, qualification_accepted=False,
        reviewed_at_utc=datetime.now(timezone.utc).isoformat(),protocol_sha256=ph,
        source_count=len(p['source_sha256']),input_count=len(p['input_sha256']),
        preserved_count=len(preserved['preserved_sha256']),paths=paths,pairs=pairs,
        regular_checkpoints=samples, live_checkpoint_samples=restarts,
        independent_full_restart_replay_performed=False,
        recorded_checks_not_independently_replayed=['Every-step physical extrema and eight-step exact restarts; saved states inspected, original run evidence preserved.',
            'Read-only instrumentation gate and independent eight-step chemical prefix are recorded runtime checks.'],
        production_seconds=timing['seconds'], new_scientific_histories=0,
        consumed_sha256=r.WATCHED, script_sha256=pin(__file__),review_seconds=time.perf_counter()-started,
        scope='Adjacent timestep screen only: feedback off, one synthetic seed-7 zygote, grid72, t0–18, unchanged limits. No spatial/control/full-fine-CPU or fresh-history qualification.')
    target = HERE/'independent_review.json'
    if target.exists(): raise FileExistsError(target)
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('review_completed','temporal_screen_passed','full_development_qualified','preserved_count','review_seconds')}),flush=True)


if __name__ == '__main__': main()
