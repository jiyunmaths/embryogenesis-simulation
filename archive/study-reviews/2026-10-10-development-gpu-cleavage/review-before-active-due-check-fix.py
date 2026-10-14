"""Independent saved-array review. No embryo solver/assessment imports.

Runtime-only extrema are explicitly distinguished from literal reconstructions.
Acceptance limits are read from the prospectively pinned protocol.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.ndimage import label

ROOT = Path('output/studies/development-gpu-cleavage-qualification-v2')
HERE = Path(__file__).parent
WATCHED = {}


def sha(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def pin(path, expected=None):
    key = str(Path(path).resolve())
    value = sha(key)
    assert expected is None or value == expected, key
    assert key not in WATCHED or WATCHED[key] == value, key
    WATCHED[key] = value
    return value


def read(path):
    pin(path)
    return json.loads(Path(path).read_text())


def arrays(path):
    pin(path)
    with np.load(path, allow_pickle=False) as archive:
        return {key: archive[key].copy() for key in archive.files}


def compare(first, second):
    errors = {key: 0. for key in ('chemical_log_max', 'polarity_abs_max',
        'relative_axis_max', 'relative_volume_max', 'relative_transport_max',
        'ancestor_chemical_log_max', 'whole_amount_relative_max')}
    unmatched = []
    assert len(first) == len(second)
    for a, b in zip(first, second):
        assert abs(a['time'] - b['time']) < 1e-9
        ua, ub = np.array(a['chemistry']), np.array(b['chemistry'])
        va, vb = np.array(a['volumes']), np.array(b['volumes'])
        qa, qb = ua @ va, ub @ vb
        values = dict(relative_axis_max=abs(a['axis_ratio'] / b['axis_ratio'] - 1),
            whole_amount_relative_max=float(np.max(abs(qa / qb - 1))),
            ancestor_chemical_log_max=float(np.max(abs(np.log(qa / va.sum()) - np.log(qb / vb.sum())))))
        if a['ids'] == b['ids']:
            da, db = np.array(a['delta']), np.array(b['delta'])
            norm = np.linalg.norm(db)
            values.update(chemical_log_max=float(np.max(abs(np.log(ua) - np.log(ub)))),
                polarity_abs_max=float(np.max(abs(np.array(a['polarity']) - np.array(b['polarity'])))),
                relative_volume_max=float(np.max(abs(va / vb - 1))),
                relative_transport_max=float(np.linalg.norm(da - db) / norm) if norm else float(np.linalg.norm(da - db)))
        else:
            unmatched.append(a['time'])
        for key, value in values.items():
            errors[key] = max(errors[key], value)
    errors['maximum_unmatched_observations'] = len(unmatched)
    return errors, unmatched


def history_checks(history, config):
    assert len(history) == 101
    maximum_boundary = maximum_rowsum = maximum_mass = maximum_symmetry = 0.
    for index, row in enumerate(history):
        assert abs(row['time'] - index * .03) < 1e-9
        count = len(row['ids'])
        assert row['ids'] in ([0], [1, 2])
        u, v = np.array(row['chemistry']), np.array(row['volumes'])
        p, d = np.array(row['polarity']), np.array(row['delta'])
        assert u.shape == (2, count) and v.shape == (count,)
        assert p.shape == (count, 3) and d.shape == (count, count)
        assert all(np.isfinite(x).all() for x in (u, v, p, d))
        assert np.min(u) > 0 and np.min(v) > 0
        scale = max(1., np.linalg.norm(d))
        mass = v[:, None] * d
        maximum_rowsum = max(maximum_rowsum, float(abs(d.sum(axis=1)).max()) / scale)
        maximum_mass = max(maximum_mass, float(abs(v @ d).max()) / scale)
        maximum_symmetry = max(maximum_symmetry, float(abs(mass - mass.T).max()) / scale)
        off = d.copy(); np.fill_diagonal(off, 0)
        assert np.min(off) >= 0 and np.max(np.diag(d)) <= 0
        assert len(row['nondividing_components']) == count - len(row['dividing_ids'])
        assert all(c == 1 for c in row['nondividing_components'])
        assert row['boundary_occupancy'] < .01
        maximum_boundary = max(maximum_boundary, row['boundary_occupancy'])
    assert max(maximum_rowsum, maximum_mass, maximum_symmetry) < 1e-12
    return dict(maximum_boundary=maximum_boundary, normalized_row_sum_error=maximum_rowsum,
        normalized_amount_transport_error=maximum_mass, normalized_mass_symmetry_error=maximum_symmetry)


def checkpoint_checks(z, config, history=None, final=False):
    meta = json.loads(str(z['metadata']))
    n = len(z['ids'])
    assert str(z['chemical_method']) == 'compensated resident progressive development v1'
    assert int(z['development_schema']) == 1
    assert str(z['attribute_mode']) == ('direct' if config['feedback'] else 'no_feedback')
    assert int(z['attribute_schema']) == 1
    assert z['phi'].dtype == np.float32 and z['phi'].shape == (n, 72, 72, 72)
    assert z['phase_carry'].dtype == np.float64 and z['phase_carry'].shape == z['phi'].shape
    assert np.isfinite(z['phase_carry']).all()
    assert np.min(z['phi']) >= 0 and np.max(z['phi']) <= 1
    effective = z['phi'].astype(np.float64) + z['phase_carry']
    assert np.min(effective) >= -1e-15 and np.max(effective) <= 1 + 1e-15
    u = np.stack((z['activator'], z['inhibitor'])); e = z['chemical_correction']
    assert u.dtype == e.dtype == np.float64 and u.shape == e.shape == (2, n)
    assert np.isfinite(u).all() and np.min(u) > 0 and np.isfinite(e).all()
    assert np.all(abs(e) <= 2 * np.spacing(u))
    assert z['polarity'].shape == (n, 3) and np.isfinite(z['polarity']).all()
    assert z['target'].shape == (n,) and np.isfinite(z['target']).all() and np.min(z['target']) > 0
    assert np.all(z['fate'] == 0) and np.isfinite(z['due']).all()
    assert meta['config'] == config
    if final:
        assert meta['time'] == 3. and meta['step_number'] == config['steps']
        assert meta['next_id'] == 3 and meta['divisions'] == {}
        np.testing.assert_array_equal(z['ids'], [1, 2])
        np.testing.assert_array_equal(z['parents'], [0, 0])
        assert len(meta['lineage']) == 3
        assert [line['id'] for line in meta['lineage']] == [0, 1, 2]
        assert [line['parent'] for line in meta['lineage']] == [-1, 0, 0]
        assert all(label(phi >= .5)[1] == 1 for phi in z['phi'])
    if history is not None:
        np.testing.assert_array_equal(u, np.array(history[-1]['chemistry']))
        np.testing.assert_array_equal(z['polarity'], np.array(history[-1]['polarity']))
        # Match the specified visible float32 occupancy followed by float64 reduction.
        h = z['phi'] * z['phi'] * (np.float32(3) - np.float32(2) * z['phi'])
        v = h.sum(axis=(1, 2, 3), dtype=np.float64) * (2 * config['extent'] / config['grid']) ** 3
        np.testing.assert_allclose(v, history[-1]['volumes'], rtol=2e-14, atol=1e-15)
    return meta


def report_checks(report, history, limits):
    audit = report['audit']
    assert report['passed'] and len(report['events']) == 1
    quality = (audit['max_volume_error'] < limits['volume_error_max'] and
        audit['min_radius'] >= limits['radius_min'] and audit['max_clipping'] == 0 and
        audit['dilution_amount_error'] <= limits['dilution_amount_error_max'])
    event = report['events'][0]
    assert event['parent'] == 0 and event['children'] == [1, 2]
    assert event['daughter_components'] == [1, 1]
    quality = quality and event['neck'] <= .15 and event['overlap'] <= .002 and .1 <= event['fraction'] <= .9
    for key, limit in [('amount_relative_jump', 'event_amount_jump_max'),
        ('target_relative_jump', 'event_target_jump_max'), ('volume_relative_jump', 'event_occupancy_jump_max'),
        ('effective_occupancy_absolute_jump', 'effective_occupancy_absolute_jump_max')]:
        quality = quality and event[key] <= limits[limit]
    return bool(quality)


def main():
    pin(__file__)
    protocol = read(ROOT / 'protocol.json'); protocol_hash = sha(ROOT / 'protocol.json')
    for path, expected in {**protocol['source_sha256'], **protocol['input_sha256']}.items():
        pin(path, expected)
    assert read(ROOT / 'status.json')['state'] in ('completed', 'completed_with_unresolved_checks')
    production = read(ROOT / 'summary.json')
    assert production['protocol_sha256'] == protocol_hash
    assert len(production['backend_rows']) == 16 and len(production['timestep_rows']) == 8
    runtime = read(ROOT / 'runtime.json'); pin(runtime['binary'], runtime['binary_sha256'])
    assert runtime['gpu'] == 'NVIDIA GeForce GTX 1080 Ti'
    for path, expected in production['evidence_sha256'].items():
        pin(path, expected)
    input_audit = read(HERE / 'input-audit.json')
    restart = read(HERE / 'restart-replay.json')
    for evidence in (input_audit, restart):
        assert evidence['passed'] and evidence['cases'] == 16 and evidence['protocol_sha256'] == protocol_hash
        pin(HERE / ('audit_inputs.py' if evidence is input_audit else 'replay_restarts.py'), evidence['script_sha256'])
    for path, expected in restart['evidence_sha256'].items():
        pin(path, expected)
    mature = read(ROOT / 'mature-regression.json')
    assert mature['passed'] and mature['exact_steps'] == 8
    pin(mature['source'], mature['source_sha256'])
    records, time_records = [], []
    histories, finals, reports = {}, {}, {}
    for job in protocol['jobs']:
        folder = ROOT / job['key']; hs = []; zs = []; rs = []; gates = []
        for arm in ('gpu', 'cpu'):
            report = read(folder / (arm + '-result.json'))
            assert report['job'] == job and report['protocol_sha256'] == protocol_hash
            for path, expected in report['evidence_sha256'].items(): pin(path, expected)
            history = read(folder / (arm + '-history.json'))
            graph_checks = history_checks(history, job['config'])
            z = arrays(folder / (arm + '-final.npz'))
            checkpoint_checks(z, job['config'], history, final=True)
            assert json.loads(str(z['event_audits'])) == report['events']
            gates.append(dict(arm=arm, runtime_quality_passed=report_checks(report, history, protocol['limits']),
                graph_checks=graph_checks))
            hs.append(history); zs.append(z); rs.append(report)
        for key in ('rng', 'fate_rng', 'signal_rng'):
            assert json.loads(str(zs[0]['metadata']))[key] == json.loads(str(zs[1]['metadata']))[key]
        np.testing.assert_array_equal(zs[0]['due'], zs[1]['due'])
        errors, unmatched = compare(*hs)
        errors['final_phi_abs_max'] = float(np.max(abs(zs[0]['phi'] - zs[1]['phi'])))
        time_error = abs(rs[0]['events'][0]['time'] - rs[1]['events'][0]['time'])
        native = read(folder / 'native-prefix.json'); prefix_errors, _ = compare(hs[0][:len(native)], native)
        prefix_phi = float(np.max(abs(arrays(folder / 'gpu-prefix-phi.npz')['phi'] - arrays(folder / 'native-prefix-phi.npz')['phi'])))
        assert rs[0]['exact_restart_steps'] == {'active': 8, 'post': 8}
        assert rs[0]['independent_chemical_steps'] == 8
        for name, count in (('gpu-started.npz', 1), ('active-restart.npz', 1), ('post-restart.npz', 2)):
            z = arrays(folder / name); m = checkpoint_checks(z, job['config'])
            assert len(z['ids']) == count
            if count == 1: assert set(m['divisions']) == {'0'}
            else:
                assert m['divisions'] == {} and all(label(f >= .5)[1] == 1 for f in z['phi'])
                source_target = arrays(job['source'])['target'].sum()
                assert abs(z['target'].sum() / source_target - 1) <= protocol['limits']['event_target_jump_max']
                np.testing.assert_array_equal(np.stack((z['activator'], z['inhibitor'])), rs[0]['events'][0]['chemistry_after'])
        ev = rs[0]['events'][0]
        vb, va = np.array(ev['volumes_before']), np.array(ev['volumes_after'])
        ub, ua = np.array(ev['chemistry_before']), np.array(ev['chemistry_after'])
        qb, qa = ub @ vb, ua @ va
        literal_amount_jump = float(np.max(abs(qa / qb - 1)))
        corrected_before = (ub.astype(np.longdouble) - np.array(ev['chemical_correction_before'], dtype=np.longdouble)) @ vb.astype(np.longdouble)
        corrected_after = (ua.astype(np.longdouble) - np.array(ev['chemical_correction_after'], dtype=np.longdouble)) @ va.astype(np.longdouble)
        represented_amount_jump = float(np.max(abs(corrected_after / corrected_before - 1)))
        literal_volume_jump = float(abs(va.sum() / vb.sum() - 1))
        passed = (all(value <= protocol['limits'][key] for key, value in errors.items()) and
            time_error <= 2 * job['dt'] + 1e-12 and all(gate['runtime_quality_passed'] for gate in gates) and
            prefix_phi <= protocol['limits']['native_prefix_phi_max'] and
            prefix_errors['chemical_log_max'] <= protocol['limits']['native_prefix_chemical_max'] and
            rs[0]['independent_chemical_log_max'] <= protocol['limits']['numerical_chemical_max'] and
            literal_amount_jump <= protocol['limits']['event_amount_jump_max'] and
            represented_amount_jump <= protocol['limits']['event_amount_jump_max'] and
            literal_volume_jump <= protocol['limits']['event_occupancy_jump_max'])
        records.append(dict(case=job['key'],passed=bool(passed),errors=errors,unmatched_times=unmatched,
            event_time_error=time_error,native_prefix_phi_error=prefix_phi,
            native_prefix_chemical_error=prefix_errors['chemical_log_max'],
            literal_gpu_event_amount_jump=literal_amount_jump,represented_gpu_event_amount_jump=represented_amount_jump,
            literal_gpu_event_volume_jump=literal_volume_jump,gates=gates))
        histories[job['key']], finals[job['key']], reports[job['key']] = hs[0], zs[0], rs[0]
    for axis in ('axial', 'oblique'):
        for branch in ('off', 'tension_adhesion', 'polarity', 'full'):
            names = [axis + '_' + branch + '_' + level for level in ('coarse', 'fine')]
            errors, unmatched = compare(*[histories[name] for name in names])
            errors['final_phi_abs_max'] = float(np.max(abs(finals[names[0]]['phi'] - finals[names[1]]['phi'])))
            event_error = abs(reports[names[0]]['events'][0]['time'] - reports[names[1]]['events'][0]['time'])
            passed = all(value <= protocol['time_limits'][key] for key, value in errors.items()) and event_error <= .0075 + 1e-12
            time_records.append(dict(axis=axis,branch=branch,passed=bool(passed),errors=errors,
                unmatched_times=unmatched,event_time_error=event_error))
    assert [r['passed'] for r in production['backend_rows']] == [r['passed'] for r in records]
    assert [r['passed'] for r in production['timestep_rows']] == [r['passed'] for r in time_records]
    preservation = read('archive/study-document-snapshots/2026-10-10-before-development-gpu/manifest.json')
    for path, expected in preservation['preserved_sha256'].items(): pin(path, expected)
    for filename, fragment in [('tests.log', '28 passed'), ('gpu-regression-tests.log', '5 passed'),
        ('development-cuda-software-tests.log', '2 passed')]:
        path = Path('archive/study-launches/2026-10-10-development-gpu') / filename
        pin(path); assert fragment in path.read_text()
    pin('tests/test_development_gpu_cuda.py')
    for path, expected in WATCHED.copy().items(): assert sha(path) == expected, path
    result = dict(passed=all(r['passed'] for r in records + time_records),
        independent_review_pending=False,protocol_sha256=protocol_hash,new_scientific_histories=0,
        scope='Prescribed single cleavage only, grid 72, two timesteps, through t=3.',
        backend_rows=records,timestep_rows=time_records,source_pins=len(protocol['source_sha256']),
        input_pins=len(protocol['input_sha256']),previous_evidence_pins=len(preservation['preserved_sha256']),
        consumed_sha256=WATCHED,script_sha256=sha(__file__),
        runtime_evidence_limits=['Every-step clipping/dilution/volume extrema and CPU event jumps are recorded runtime evidence.',
            'Intermediate daughter connectivity uses recorded component counts; immediate GPU and final CPU/GPU fields were relabeled independently.',
            'Restart replays use the candidate solver and qualify restoration, not independent physics.',
            'Two GPU software tests of shape-aligned/repeated cleavage do not qualify developmental spatial/time accuracy.'])
    (HERE / 'independent_review.json').write_text(json.dumps(result,indent=2,allow_nan=False) + '\n')
    print(json.dumps(dict(passed=result['passed'],backend_cases=len(records),timestep_pairs=len(time_records),
        source_pins=result['source_pins'],input_pins=result['input_pins'],previous_evidence_pins=result['previous_evidence_pins'],
        consumed_pins=len(WATCHED)),indent=2))


if __name__ == '__main__':
    main()
