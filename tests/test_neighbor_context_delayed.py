import json

import numpy as np
import pytest

from embryo import neighbor_context_delayed as study
from embryo.cell_response import pulse
from embryo.feedback_long import digest
from embryo.resolution import write_json
from embryo.transport import conservative_transport


def fixture(root, kinds=('challenged',)):
    root.mkdir()
    masses = np.array([1., 2., 3.])
    conductance = np.array([[0., .1, .2], [.1, 0., .3], [.2, .3, 0.]])
    delta = conservative_transport(conductance, masses).delta.toarray()
    initial = np.ones((2, 3))
    original = np.array([[.5, .8, 1.2], [.7, .9, 1.1]])
    jobs = []
    for kind in kinds:
        folder = root/kind; folder.mkdir()
        source = folder/'source.npz'
        arrays = dict(initial=initial, original_initial=original, ids=np.array([19, 20, 21]),
                      delta=delta, masses=masses)
        np.savez(source, **arrays)
        job = dict(key=kind, kind=kind, case='case', target=19, target_index=0,
                   beta=2., da=.02, db=.04, source=str(source))
        fun, jac = study.dynamics(job, arrays)
        job['initial_stability'] = study.endpoint_metrics(initial, fun, jac, study.CRITERIA)
        jobs.append(job)
    p = dict(jobs=jobs, times=[0., .1, .5, 1., 2., 3., 5., 10., 24., 60.], criteria=study.CRITERIA,
             total_jobs=len(jobs), trajectories=5*len(jobs), independent_histories=1, cases=[],
             interpretation='fixture', limits='fixture', source_sha256={},
             input_sha256={j['source']:digest(j['source']) for j in jobs})
    write_json(root/'protocol.json', p)
    return p, jobs


def test_matching_absolute_dose_does_not_confuse_fraction_with_amount():
    conditions = study.pulse_conditions(1.28, .05, 3.7)
    assert [c['factor'] for c in conditions[:2]] == [.9, 1.1]
    for c in conditions[2:]:
        assert c['injected_activator_amount'] == pytest.approx(c['original_activator_amount'], abs=1e-15)
        assert abs(c['factor']-1) < .004
    assert abs(conditions[0]['injected_activator_amount']) > 25*abs(conditions[2]['injected_activator_amount'])


@pytest.mark.parametrize('current,original,volume', [(0., 1., 1.), (np.nan, 1., 1.),
    (1., 0., 1.), (1., 1., -1.), (.05, 1., 1.)])
def test_invalid_or_nonpositive_matched_doses_are_rejected(current, original, volume):
    with pytest.raises(ValueError): study.pulse_conditions(current, original, volume)


def test_copying_full_endpoint_preserves_collective_state_and_target_alignment(tmp_path):
    p, jobs = fixture(tmp_path/'study'); job = jobs[0]
    first = study.checked_snapshot(job)
    first['initial'][0, 1] = 1.2
    second = study.checked_snapshot(job)
    np.testing.assert_array_equal(second['initial'], np.ones((2, 3)))
    assert second['original_initial'][0, 0] == .5
    assert job['target_index'] == list(second['ids']).index(job['target'])
    job['target_index'] = 1
    with pytest.raises(ValueError, match='alignment'): study.checked_snapshot(job)


def test_real_worker_independent_solvers_target_only_doses_and_resume(tmp_path):
    root = tmp_path/'study'; p, jobs = fixture(root); job = jobs[0]
    record = study.worker((root, job)); saved, arrays = study.load_checked(root, p, job)
    assert saved == record and record['max_solver_log_error'] < 1e-8
    assert record['control_drift_log_max'] < 1e-10
    assert all(t['target_recovery_time'] is not None and t['network_recovery_time'] is not None for t in record['trials'])
    assert all(t['classification'] == 'returned' for t in record['trials'])
    for f, condition in enumerate(record['trials']):
        np.testing.assert_array_equal(arrays['pulse_paths'][f, 0], pulse(arrays['initial'], 0, condition['factor']))
        if condition['dose'] == 'matched_amount':
            assert condition['injected_activator_amount'] == pytest.approx(condition['original_activator_amount'], abs=1e-15)
    before = digest(root/job['key']/'trajectories.npz')
    assert study.worker((root, job)) == record
    assert digest(root/job['key']/'trajectories.npz') == before


@pytest.mark.parametrize('kind', ['truncated', 'nonpositive', 'wrong_target', 'wrong_amount',
                                 'wrong_wave', 'wrong_recovery', 'wrong_classification', 'solver_nan'])
def test_assessment_rejects_corruption_even_with_pass_flag_and_new_file_hash(tmp_path, kind):
    root = tmp_path/'study'; p, jobs = fixture(root); job = jobs[0]; study.worker((root, job))
    folder = root/job['key']; filename = folder/'trajectories.npz'
    record = study.read(folder/'result.json')
    with np.load(filename) as d: arrays = {k:d[k].copy() for k in d.files}
    if kind == 'truncated': arrays['pulse_paths'] = arrays['pulse_paths'][:, :-1]
    elif kind == 'nonpositive': arrays['control'][1, 0, 0] = 0.
    elif kind == 'wrong_target': arrays['pulse_paths'][0, 0, 0, 1] *= .9
    elif kind == 'wrong_amount': record['trials'][2]['injected_activator_amount'] *= 2.
    elif kind == 'wrong_wave': arrays['responses'][0, 1, 0] += 1.
    elif kind == 'wrong_recovery': record['trials'][0]['network_recovery_time'] = None
    elif kind == 'wrong_classification': record['trials'][0]['classification'] = 'different_stable_endpoint'
    else: record['max_solver_log_error'] = float('nan')
    np.savez(filename, **arrays); record['trajectory_sha256'] = digest(filename)
    # Deliberately bypass the production writer's own NaN rejection to verify
    # that the reader also rejects a nonstandard/corrupt JSON value.
    (folder/'result.json').write_text(json.dumps(record))
    with pytest.raises(ValueError): study.assess(root)


def test_source_mutation_is_rejected_before_reusing_a_completed_job(tmp_path):
    root = tmp_path/'study'; _, jobs = fixture(root); job = jobs[0]; study.worker((root, job))
    with open(job['source'], 'ab') as file: file.write(b'changed')
    with pytest.raises(ValueError, match='Frozen input changed'): study.worker((root, job))


def test_unsettled_source_cannot_be_used_as_a_delayed_equilibrium(tmp_path):
    root = tmp_path/'study'; p, jobs = fixture(root); job = jobs[0]
    with np.load(job['source']) as d: arrays = {k:d[k].copy() for k in d.files}
    arrays['initial'][0, 0] = .5; np.savez(job['source'], **arrays)
    p['input_sha256'][job['source']] = digest(job['source']); write_json(root/'protocol.json', p)
    with pytest.raises(ValueError, match='stable equilibrium'): study.worker((root, job))


def test_parallel_run_reconstructs_all_endpoint_backgrounds_and_readonly_assessment(tmp_path):
    root = tmp_path/'study'; p, _ = fixture(root, study.KINDS)
    study.run(root, workers=2)
    status = study.read(root/'status.json'); result = study.read(root/'results.json')
    assert status['state'] == 'completed' and status['completed'] == 3
    assert result['numerical_pass'] and result['trajectories'] == 15
    assert all(result['counts'][kind]['returned'] == 4 for kind in study.KINDS)
    before = digest(root/'results.json'); assert study.assess(root, write=False) == result
    assert digest(root/'results.json') == before


def test_incomplete_screen_cannot_be_followed_and_existing_outputs_are_preserved(tmp_path):
    with pytest.raises(FileExistsError): study.prepare(tmp_path)
    source = tmp_path/'source'; source.mkdir()
    write_json(source/'protocol.json', dict(source_sha256={}, input_sha256={}))
    write_json(source/'status.json', dict(state='running', completed=59, numerical_pass=True))
    write_json(source/'results.json', dict(numerical_pass=True))
    with pytest.raises(ValueError, match='complete accepted'): study.prepare(tmp_path/'new', source)
    assert not (tmp_path/'new').exists()
