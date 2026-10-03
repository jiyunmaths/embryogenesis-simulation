import numpy as np
import pytest

pytest.importorskip('torch')
from embryo import exchange_response_histories as study
from embryo.attribute_exchange import exchange
from embryo.cell_response_moving import name
from embryo.feedback_long import digest
from embryo.resolution import write_json


def test_prespecified_selection_and_conservative_exchange_with_unequal_volumes():
    base = np.array([[1., 4., 2.], [2., 8., 3.]])
    masses, ids = np.array([1., 1.4, .7]), np.array([19, 20, 21])
    fresh, pair, error = study.exchange_states(base, base.copy(), ids, masses, [0, 1], [19, 20])
    assert pair == [0, 1] and error < 2e-14
    np.testing.assert_allclose(fresh@masses, base@masses, rtol=2e-14)
    np.testing.assert_array_equal(fresh[:, 2], base[:, 2])
    np.testing.assert_array_equal(base, [[1., 4., 2.], [2., 8., 3.]])
    with pytest.raises(ValueError, match='selection'):
        study.exchange_states(base, base, ids, masses, [1, 0], [20, 19])
    # Activity ties use cell IDs, never a subsequent response measurement.
    tied = np.ones((2, 3))
    _, pair, _ = study.exchange_states(tied, tied, ids[::-1], masses, [2, 0], [19, 21])
    assert pair == [2, 0]


@pytest.mark.parametrize('field,bad', [
    ('base', np.array([[1., 0.], [1., 1.]])),
    ('relaxed', np.array([[1., np.nan], [1., 1.]])),
    ('ids', np.array([19, 19])), ('masses', np.array([1., -1.])),
    ('masses', np.ones(3)), ('base', np.ones((3, 2))),
])
def test_invalid_exchange_inputs_are_rejected(field, bad):
    args = dict(base=np.array([[1., 4.], [1., 8.]]), relaxed=np.ones((2, 2)),
                ids=np.array([19, 20]), masses=np.ones(2), expected_pair=[0, 1], expected_ids=[19, 20])
    args[field] = bad
    with pytest.raises(ValueError, match='compartments'):
        study.exchange_states(**args)


def protocol():
    return dict(seeds=[8, 9], histories=[dict(seed=s, pair=[0, 1], targets=[19, 20]) for s in (8, 9)],
                backgrounds=study.BACKGROUNDS, factors=[.9, 1.1], formation_start=150.,
                formation_duration=4., response_start=154., response_duration=4., interval=.5,
                criteria=dict(initial_pair_informative_min=.1, pair_return_ratio=.1,
                              late_duration=1., log_activator_sd_min=.1), scope='test', limits='pilot',
                source_sha256={}, input_sha256={})


def rows(path, times, start):
    return [dict(elapsed=float(t), time=float(start+t), ids=[19, 20], chemistry=x.tolist(),
                 volumes=[1., 1.], polarity=[[.1, 0., 0.], [.1, 0., 0.]],
                 delta=[[-1., 1.], [1., -1.]], axis_ratio=1.2, boundary_occupancy=0.)
            for t, x in zip(times, path)]


def write_trial(child, job, path, times, start):
    folder = child/name(job); folder.mkdir()
    write_json(folder/'result.json', dict(quality_pass=True, protocol_sha256=digest(child/'protocol.json')))
    write_json(folder/'history.json', rows(path, times, start))
    return folder


def formation_fixture(root, seed=8, outcome='donor', uniform=False):
    p = protocol(); times = np.arange(9)*p['interval']
    source = root/f'seed-{seed}'/'formation-source'; source.mkdir(parents=True)
    child = root/f'seed-{seed}'/'formation'; child.mkdir()
    masses, ids = np.ones(2), np.array([19, 20])
    base = np.ones((2, 2)) if uniform else np.array([[1., 4.], [1., 8.]])
    fresh = exchange(base, masses, 0, 1)
    initial = dict(unexchanged=base, fresh_exchange=fresh, relaxed_exchange=fresh.copy())
    np.savez(source/'initial_states.npz', ids=ids, masses=masses, **initial)
    np.savez(child/'initial_states.npz', ids=ids, masses=masses, **initial)
    frozen = {k: np.repeat(v[None], len(times), axis=0) for k, v in initial.items()}
    np.savez(source/'frozen_reference.npz', times=times, **frozen)
    write_json(child/'protocol.json', dict(start=p['formation_start'], source_sha256={}, input_sha256={}))
    for family in study.BACKGROUNDS:
        path = frozen[family].copy()
        if family != 'unexchanged':
            if outcome == 'destination': path[1:] = base
            elif outcome == 'reorganized': path[1:] = np.sqrt(base*fresh)
        write_trial(child, dict(family=family, target=None, factor=1.), path, times, p['formation_start'])
    return p, child, source, times


@pytest.mark.parametrize('outcome', ['donor', 'destination', 'reorganized'])
def test_formation_distinguishes_transfer_return_and_reorganization(tmp_path, outcome):
    p, _, _, _ = formation_fixture(tmp_path, outcome=outcome)
    report = study.formation_assessment(tmp_path, 8, p)
    untouched, fresh, relaxed = report['trials']
    assert untouched['destination_like'] and untouched['contrast_retained']
    for row in (fresh, relaxed):
        assert row['destination_like'] == (outcome == 'destination')
        assert row['transferred_like'] == (outcome == 'donor')
        assert row['late_pair_order_reversed'] == (outcome == 'donor')
        if outcome == 'reorganized': assert not row['contrast_retained']


def test_uninformative_uniform_pair_is_reported_without_dividing_by_zero(tmp_path):
    p, _, _, _ = formation_fixture(tmp_path, uniform=True)
    for row in study.formation_assessment(tmp_path, 8, p)['trials']:
        assert not row['informative_pair'] and not row['contrast_retained']
        assert row['late_destination_ratio'] is None and row['late_transferred_ratio'] is None
        assert not row['destination_like'] and not row['transferred_like']


@pytest.mark.parametrize('fault', ['partial', 'elapsed', 'clock', 'ids', 'nan', 'zero', 'shape', 'volume', 'boundary', 'failed', 'hash'])
def test_claimed_completed_formation_rejects_invalid_evidence(tmp_path, fault):
    p, child, _, times = formation_fixture(tmp_path)
    job = dict(family='fresh_exchange', target=None, factor=1.)
    folder = child/name(job); history = study.read(folder/'history.json')
    if fault == 'partial': history.pop()
    elif fault == 'elapsed': history[-1]['elapsed'] += .01
    elif fault == 'clock': history[-1]['time'] += .01
    elif fault == 'ids': history[-1]['ids'].reverse()
    elif fault == 'nan': history[-1]['chemistry'][0][0] = float('nan')
    elif fault == 'zero': history[-1]['chemistry'][0][0] = 0.
    elif fault == 'shape': history[-1]['chemistry'] = [[1., 1.]]
    elif fault == 'volume': history[-1]['volumes'][0] = -1.
    elif fault == 'boundary': history[-1]['boundary_occupancy'] = float('nan')
    else:
        result = study.read(folder/'result.json')
        if fault == 'failed': result['quality_pass'] = False
        else: result['protocol_sha256'] = 'stale'
        write_json(folder/'result.json', result)
    # Deliberately inject malformed external evidence; production JSON disallows NaN.
    (folder/'history.json').write_text(json.dumps(history))
    with pytest.raises(ValueError):
        study.completed_history(child, job, times)


@pytest.mark.parametrize('fault', ['times', 'chemistry', 'initial'])
def test_formation_rejects_changed_frozen_or_initial_reference(tmp_path, fault):
    p, child, source, _ = formation_fixture(tmp_path)
    file = source/('initial_states.npz' if fault == 'initial' else 'frozen_reference.npz')
    with np.load(file) as data: changed = {k: data[k].copy() for k in data.files}
    if fault == 'times': changed['times'][-1] += .1
    elif fault == 'chemistry': changed['fresh_exchange'][-1, 0, 0] = -1.
    else: changed['fresh_exchange'][0, 0] *= 1.1
    np.savez(file, **changed)
    with pytest.raises(ValueError, match='mismatch'):
        study.formation_assessment(tmp_path, 8, p)


def response_fixture(root, seed, transferred):
    p = protocol(); times = np.arange(9)*p['interval']
    for family in study.BACKGROUNDS:
        child = root/f'seed-{seed}'/'response'/family; child.mkdir(parents=True)
        jobs = [dict(family=family, target=None, factor=1.)]+[
            dict(family=family, target=cell, factor=factor) for cell in [19, 20] for factor in [.9, 1.1]]
        write_json(child/'protocol.json', dict(start=p['response_start'], jobs=jobs))
        np.savez(child/'initial_states.npz', ids=[19, 20], masses=np.ones(2))
        control = np.exp(.003*times[:, None, None])*np.ones((len(times), 2, 2))
        for job in jobs:
            path = control.copy()
            if job['target'] is not None:
                i = [19, 20].index(job['target'])
                ref = 1-i if transferred and family != 'unexchanged' else i
                tau = [1., 3.][ref]
                path[:, 0, i] *= np.exp(np.log(job['factor'])*np.exp(-times/tau))
            write_trial(child, job, path, times, p['response_start'])


def test_aggregate_keeps_independent_histories_separate_from_nested_interventions(tmp_path):
    p, _, _, _ = formation_fixture(tmp_path, seed=8)
    formation_fixture(tmp_path, seed=9, outcome='destination')
    response_fixture(tmp_path, 8, transferred=True)
    response_fixture(tmp_path, 9, transferred=False)
    anchor = tmp_path/'historical-seed7'; anchor.mkdir()
    write_json(anchor/'comparison.json', dict(comparisons=[dict(nearest_reference='donor')]*8))
    p['historical_seed7'] = str(anchor)
    write_json(tmp_path/'protocol.json', p)
    report = study.assess(tmp_path)
    assert report['independent_histories'] == 3
    assert report['new_histories'][0]['descriptive_counts'] == dict(donor=8, destination=0, unresolved=0)
    assert report['new_histories'][1]['descriptive_counts'] == dict(donor=0, destination=8, unresolved=0)
    assert len(report['historical_seed7_comparisons']) == 8
    # A worker's quality flag cannot conceal malformed geometry during aggregation.
    folder = tmp_path/'seed-9/response/fresh_exchange'/name(dict(family='fresh_exchange', target=19, factor=.9))
    h = study.read(folder/'history.json'); h[-1]['delta'][0][0] = float('nan')
    (folder/'history.json').write_text(json.dumps(h))
    with pytest.raises(ValueError, match='Nonfinite'):
        study.assess(tmp_path)


def test_failed_gpu_gate_leaves_no_prepared_study(tmp_path, monkeypatch):
    def failed(_): raise ValueError('GPU validation failed')
    monkeypatch.setattr(study, 'require_validation', failed)
    root = tmp_path/'study'
    with pytest.raises(ValueError, match='GPU validation'):
        study.prepare(root)
    assert not root.exists()


def test_new_history_backend_disagreement_blocks_all_scientific_runs(tmp_path, monkeypatch):
    p = protocol(); p['validation'] = 'gate'
    write_json(tmp_path/'protocol.json', p)
    monkeypatch.setattr(study, 'require_validation', lambda _: (None, []))
    calls = []
    def check(root, seed, family, _):
        calls.append((seed, family))
        if family == 'fresh_exchange': raise RuntimeError('backend disagreement')
    def forbidden(*args): raise AssertionError('Scientific run bypassed failed starting-state check')
    monkeypatch.setattr(study, 'prefix_check', check)
    monkeypatch.setattr(study, 'gpu_prepare', forbidden)
    monkeypatch.setattr(study, 'gpu_run', forbidden)
    with pytest.raises(RuntimeError, match='backend disagreement'):
        study.run(tmp_path)
    status = study.read(tmp_path/'status.json')
    assert status['state'] == 'failed' and status['completed'] == 0 and status['prefix_completed'] == 1
    assert calls == [(8, 'unexchanged'), (8, 'fresh_exchange')]


def test_saved_prefix_evidence_cannot_be_reused_after_mutation(tmp_path):
    source = tmp_path/'seed-8/formation-source'; source.mkdir(parents=True)
    write_json(source/'protocol.json', {})
    folder = tmp_path/'seed-8/prefix/fresh_exchange'; folder.mkdir(parents=True)
    evidence = folder/'gpu-history.json'; write_json(evidence, [])
    result = dict(passed=True, source_protocol_sha256=digest(source/'protocol.json'),
                  evidence_sha256={str(evidence): digest(evidence)})
    write_json(folder/'comparison.json', result)
    assert study.prefix_check(tmp_path, 8, 'fresh_exchange', {}) == result
    write_json(evidence, [dict(changed=True)])
    with pytest.raises(ValueError, match='evidence changed'):
        study.prefix_check(tmp_path, 8, 'fresh_exchange', {})


def test_response_handoff_preserves_actual_moving_endpoint_and_rejects_mismatch(tmp_path):
    from embryo.attribute_development import AttributeSimulation
    from embryo.model import Config
    p = protocol()
    sim = AttributeSimulation(Config(grid=24, max_cells=2, interface_width=.12, dt=.005), 'direct')
    sim.phi = np.repeat(sim.phi, 2, axis=0)
    sim.ids = np.array([19, 20]); sim.parents = np.array([-1, -1]); sim.due = np.full(2, np.inf)
    sim.fate = np.zeros(2); sim.activator = np.array([1., 4.]); sim.inhibitor = np.array([1., 8.])
    sim.polarity = np.array([[.1, 0., 0.], [0., .2, 0.]])
    sim.target = sim.volumes(); sim.next_id = 21; sim.step_number = 30800; sim.time = 154.
    formation = tmp_path/'seed-8/formation'; formation.mkdir(parents=True)
    write_json(tmp_path/'protocol.json', p)
    write_json(formation/'protocol.json', {})
    folder = formation/'fresh_exchange_control'; folder.mkdir()
    sim.checkpoint(folder/'latest_state.npz')
    write_json(folder/'result.json', dict(quality_pass=True, protocol_sha256=digest(formation/'protocol.json')))
    history = rows(np.array([[sim.activator, sim.inhibitor]]), [4.], 150.)
    write_json(folder/'history.json', history)
    p['dt'] = .005
    child = study.materialize_response(tmp_path, 8, 'fresh_exchange', p)
    assert digest(child/'source.npz') == digest(folder/'latest_state.npz')
    cp = study.read(child/'protocol.json'); assert len(cp['jobs']) == 5
    restarted = AttributeSimulation.restore(child/'source.npz')
    assert restarted.time == 154. and restarted.rng.bit_generator.state == sim.rng.bit_generator.state
    for key in ('phi', 'activator', 'inhibitor', 'polarity', 'ids', 'target'):
        np.testing.assert_array_equal(getattr(restarted, key), getattr(sim, key))
    assert study.materialize_response(tmp_path, 8, 'fresh_exchange', p) == child
    # Prepare a different background from a deliberately inconsistent endpoint.
    bad = formation/'relaxed_exchange_control'; bad.mkdir()
    sim.checkpoint(bad/'latest_state.npz')
    write_json(bad/'result.json', dict(quality_pass=True, protocol_sha256=digest(formation/'protocol.json')))
    history[0]['chemistry'][0][0] *= 1.1; write_json(bad/'history.json', history)
    with pytest.raises(ValueError, match='mismatch'):
        study.materialize_response(tmp_path, 8, 'relaxed_exchange', p)
    assert not (tmp_path/'seed-8/response-source/relaxed_exchange').exists()
import json
