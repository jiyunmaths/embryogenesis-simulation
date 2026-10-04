import json
from pathlib import Path

import numpy as np
import pytest

from embryo import neighbor_context as study
from embryo.cell_response import pulse, response_metrics
from embryo.feedback_long import digest
from embryo.resolution import write_json
from embryo.transport import conservative_transport


def graph():
    masses = np.array([1., 2., 3., 4.])
    g = np.array([[0., 1., 2., 0.], [1., 0., 0., 1.],
                  [2., 0., 0., 0.], [0., 1., 0., 0.]])
    delta = conservative_transport(g, masses).delta.toarray()
    initial = np.array([[.2, .8, 2., 1.2], [.4, 1., 3., 1.6]])
    ids = np.array([10, 11, 12, 13])
    return initial, masses, delta, ids


def test_neighbor_mix_conserves_unequal_volume_amounts_and_preserves_target():
    x, masses, delta, ids = graph(); before = x.copy(); old = delta.copy()
    changed, meta = study.challenge(x, masses, delta, ids, 10, 'neighbor_mix')
    np.testing.assert_allclose(changed@masses, x@masses, rtol=1e-14)
    np.testing.assert_allclose(changed[:, 1:3], np.repeat(((x[:, 1:3]@masses[1:3])/5.)[:, None], 2, axis=1))
    np.testing.assert_array_equal(changed[:, [0, 3]], x[:, [0, 3]])
    np.testing.assert_array_equal(x, before); np.testing.assert_array_equal(delta, old)
    assert meta['neighbors'] == [11, 12] and meta['neighbor_log_rms'] > .05
    assert meta['target_unchanged'] and meta['nonneighbors_unchanged']


@pytest.mark.parametrize('arm,column', [('neighbor_low', 0), ('neighbor_high', 2)])
def test_reset_changes_only_neighbors_and_records_external_amounts(arm, column):
    x, masses, delta, ids = graph()
    changed, meta = study.challenge(x, masses, delta, ids, 10, arm)
    np.testing.assert_array_equal(changed[:, 1], x[:, column])
    np.testing.assert_array_equal(changed[:, 2], x[:, column])
    np.testing.assert_array_equal(changed[:, [0, 3]], x[:, [0, 3]])
    np.testing.assert_allclose(meta['amount_change'], (changed-x)@masses)
    assert meta['relative_amount_change'] > 0 and meta['template_cell'] == ids[column]
    # Reactions are initially identical at the target: the derivative change
    # must be entirely the exact conservative transport change.
    fun = study.rhs(delta, 2., .02, .4)
    actual = (fun(0., changed.ravel())-fun(0., x.ravel())).reshape(2, -1)[:, 0]
    np.testing.assert_allclose(actual, np.array([.02, .4])*meta['target_transport_change'], atol=1e-15)


def test_untouched_sham_and_uniform_neighbor_mix_are_exact_noops():
    x, m, delta, ids = graph(); x[:, 2] = x[:, 1]
    for arm in ('untouched', 'sham', 'neighbor_mix'):
        changed, meta = study.challenge(x, m, delta, ids, 10, arm)
        np.testing.assert_array_equal(changed, x)
        assert meta['neighbor_log_rms'] == 0


def test_template_tie_break_uses_ids_and_never_assigns_target():
    x, m, delta, ids = graph(); x[0, 0] = x[0, 3]; ids[3] = 2
    _, meta = study.challenge(x, m, delta, ids, 10, 'neighbor_low')
    assert meta['template_cell'] == 11
    x[0, 1] = x[0, 0]
    _, meta = study.challenge(x, m, delta, ids, 10, 'neighbor_low')
    assert meta['template_cell'] == 2


@pytest.mark.parametrize('kind', ['negative', 'nonfinite', 'bad_volume', 'wrong_ids', 'missing_target', 'bad_arm', 'nonconservative', 'nonreciprocal'])
def test_invalid_inputs_and_nonconservative_graphs_are_rejected(kind):
    x, m, delta, ids = graph(); target, arm = 10, 'neighbor_mix'
    if kind == 'negative': x[0, 1] = 0.
    elif kind == 'nonfinite': x[0, 1] = np.nan
    elif kind == 'bad_volume': m[0] = 0.
    elif kind == 'wrong_ids': ids[1] = ids[0]
    elif kind == 'missing_target': target = 99
    elif kind == 'bad_arm': arm = 'unknown'
    elif kind == 'nonconservative': delta[0, 0] *= 2.
    else:
        # Row sums still zero; unequal-volume reciprocity does not hold.
        delta = np.array([[-2., 1., 1., 0.], [1., -2., 0., 1.], [1., 0., -1., 0.], [0., 1., 0., -1.]])
    with pytest.raises(ValueError): study.challenge(x, m, delta, ids, target, arm)


def fixture(root, uniform_neighbors=False):
    root.mkdir(); x, m, delta, ids = graph()
    if uniform_neighbors: x[:, 2] = x[:, 1]
    source = root/'source.npz'; np.savez(source, initial=x, masses=m, delta=delta, ids=ids)
    jobs = [dict(target=t, arm=a) for t in (10, 12) for a in study.ARMS]
    c = dict(key='context', seed=7, background='unexchanged', source=str(source),
             targets=[10, 12], beta=2., da=.02, db=.4, jobs=jobs)
    times = np.array([0., 1., 24., 60.])
    p = dict(contexts=[c], times=times.tolist(), factors=study.FACTORS, arms=study.ARMS, criteria=study.CRITERIA,
        source_sha256={}, input_sha256={str(source):digest(source)}, independent_histories=1,
        interpretation='fixture', limits='fixture', total_jobs=len(jobs))
    write_json(root/'protocol.json', p); ph = digest(root/'protocol.json')
    (root/'context').mkdir()
    for job in jobs:
        folder = root/'context'/study.job_name(job); folder.mkdir()
        changed, meta = study.challenge(x, m, delta, ids, **job)
        control = np.broadcast_to(changed, (len(times), *changed.shape)).copy()
        i = meta['target_index']; paths = []
        for factor in study.FACTORS:
            path = control.copy()
            path[:, 0, i] *= np.exp(np.log(factor)*np.exp(-times*(.4 if i == 0 else .2)))
            path[0] = pulse(changed, i, factor)
            paths.append(path)
        paths = np.asarray(paths)
        responses = np.asarray([np.log(paths[f, :, :, i]/control[:, :, i])/abs(np.log(factor))
                                for f, factor in enumerate(study.FACTORS)])
        filename = folder/'trajectories.npz'
        np.savez(filename, initial=changed, control=control, pulse_paths=paths, responses=responses, times=times)
        write_json(folder/'result.json', dict(protocol_sha256=ph, trajectory_sha256=digest(filename),
            quality_pass=True, max_solver_log_error=0., job=job, context=c['key'], intervention=meta,
            stationary_stable=True, trials=[dict(factor=factor,
                injected_activator_amount=float((factor-1)*changed[0, i]*m[i]),
                **response_metrics(paths[f], control, times, m, i, factor)) for f, factor in enumerate(study.FACTORS)]))
    return p, c


def test_matched_assessment_does_not_invent_identity_from_neighbor_change(tmp_path):
    root = tmp_path/'run'; p, c = fixture(root)
    report = study.compare_context(root, c, p)
    assert report['numerical_pass'] and report['sham_max_error'] == 0
    row = next(r for r in report['comparisons'] if r['target'] == 10 and r['arm'] == 'neighbor_low')
    assert row['informative'] and row['state_outcome'] == 'near_baseline'
    assert all(r['own_response_distance'] == 0 and not r['response_shift'] for r in row['responses'])
    assert all(r['nearest_baseline_response'] == 'own' for r in row['responses'])


def test_weak_neighbor_challenge_is_explicitly_uninformative(tmp_path):
    root = tmp_path/'run'; p, c = fixture(root, uniform_neighbors=True)
    report = study.compare_context(root, c, p)
    row = next(r for r in report['comparisons'] if r['target'] == 10 and r['arm'] == 'neighbor_mix')
    assert not row['informative'] and row['state_outcome'] == 'uninformative'
    assert all(r['nearest_baseline_response'] == 'unresolved' for r in row['responses'])


@pytest.mark.parametrize('kind', ['wrong_target', 'partial', 'nonpositive', 'wrong_waveform', 'pulse_amount'])
def test_assessment_rejects_corrupt_evidence_despite_passing_worker_flag(tmp_path, kind):
    root = tmp_path/'run'; p, c = fixture(root)
    folder = root/'context/cell-10_neighbor_low'; file = folder/'trajectories.npz'
    with np.load(file) as d: arrays = {k:d[k].copy() for k in d.files}
    if kind == 'wrong_target': arrays['initial'][0, 0] *= 2.
    elif kind == 'partial': arrays['control'] = arrays['control'][:-1]
    elif kind == 'nonpositive': arrays['control'][1, 0, 0] = 0.
    elif kind == 'wrong_waveform': arrays['responses'][0, 1, 0] += .1
    else: arrays['pulse_paths'][0, 0, 0, 0] *= .9
    np.savez(file, **arrays)
    record = study.read(folder/'result.json'); record['trajectory_sha256'] = digest(file)
    write_json(folder/'result.json', record)
    with pytest.raises(ValueError): study.compare_context(root, c, p)


def test_real_worker_checks_independent_solvers_and_reuses_unchanged_results(tmp_path):
    root = tmp_path/'run'; p, c = fixture(root)
    job = dict(target=10, arm='neighbor_mix'); folder = root/'context'/study.job_name(job)
    (folder/'result.json').unlink()
    result = study.worker((root, c, job))
    assert result['quality_pass'] and result['max_solver_log_error'] < 1e-8
    assert result['intervention']['target_unchanged']
    for row in result['trials']:
        assert row['injected_activator_amount'] == pytest.approx((row['factor']-1)*.2)
    before = digest(folder/'trajectories.npz')
    assert study.worker((root, c, job)) == result
    assert digest(folder/'trajectories.npz') == before


def test_changed_recovery_flag_is_rejected_even_if_arrays_are_unchanged(tmp_path):
    root = tmp_path/'run'; p, c = fixture(root)
    file = root/'context/cell-10_neighbor_low/result.json'
    record = study.read(file); record['trials'][0]['target_recovery_time'] = None
    write_json(file, record)
    with pytest.raises(ValueError, match='recovery'): study.compare_context(root, c, p)


def test_parallel_runner_finishes_real_solutions_and_assesses_all_arms(tmp_path):
    root = tmp_path/'run'; p, c = fixture(root)
    for result in (root/'context').glob('*/result.json'): result.unlink()
    study.run(root, workers=2)
    status = study.read(root/'status.json'); result = study.read(root/'results.json')
    assert status['state'] == 'completed' and status['completed'] == 10 and status['numerical_pass']
    assert result['numerical_pass'] and len(result['contexts'][0]['comparisons']) == 6
    assert result['contexts'][0]['sham_max_error'] == 0.


def test_prepare_cannot_bypass_refinement_or_overwrite_results(tmp_path):
    with pytest.raises(FileExistsError): study.prepare(tmp_path)
    refined = tmp_path/'refined'; refined.mkdir()
    write_json(refined/'protocol.json', dict(input_sha256={}, source_sha256={}))
    write_json(refined/'status.json', dict(state='completed', completed=18, passed=False))
    with pytest.raises(ValueError, match='refinement'):
        study.prepare(tmp_path/'new', refined=refined)
    assert not (tmp_path/'new').exists()


@pytest.mark.parametrize('seed,status,quality,accepted', [
    (7, None, True, True), (8, None, True, False), (7, None, False, False),
    (9, 'running', True, False), (9, 'completed', True, True)])
def test_legacy_source_requires_valid_control_and_new_sources_require_completion(tmp_path, seed, status, quality, accepted):
    cp = dict(source_sha256={}, input_sha256={})
    write_json(tmp_path/'protocol.json', cp)
    if status is not None: write_json(tmp_path/'status.json', dict(state=status))
    control = tmp_path/'unexchanged_control'; control.mkdir()
    write_json(control/'result.json', dict(quality_pass=quality, protocol_sha256=digest(tmp_path/'protocol.json')))
    if accepted: assert study.require_source(tmp_path, seed, 'unexchanged') == cp
    else:
        with pytest.raises(ValueError): study.require_source(tmp_path, seed, 'unexchanged')
