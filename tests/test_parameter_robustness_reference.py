import numpy as np
import pytest

from embryo import parameter_robustness_reference as study
from embryo.feedback_long import digest
from embryo.resolution import write_json


def fixture(root):
    root.mkdir()
    graph = root/'graph.npz'
    delta = np.array([[-3., 3.], [1.5, -1.5]])
    initial = np.array([[.999, 1.0005], [1., 1.]])
    np.savez(graph, delta=delta, masses=[1., 2.], ids=[3, 4], uniform=initial, pattern=initial)
    moving = root/'moving'; moving.mkdir()
    write_json(moving/'protocol.json', dict(fixture=True))
    jobs = [dict(key=f'chi-{chi}', seed=7, point=dict(ratio=20., chi=chi), family='uniform', graph=str(graph))
            for chi in (0., .7)]
    p = dict(jobs=jobs, duration=12., interval=2., late_window=4., late_log_sd_min=.1,
             moving=str(moving), parent_protocol_sha256=digest(moving/'protocol.json'),
             source_sha256={}, input_sha256={str(graph):digest(graph)})
    write_json(root/'protocol.json', p)
    return p, initial


def test_equal_horizon_reference_is_independent_of_frozen_polarity_contrast(tmp_path):
    root = tmp_path/'reference'; p, initial = fixture(root)
    result = study.run(root)
    assert result['completed'] == 2 and result['max_solver_log_error'] < 1e-8
    assert result['horizon_limited_initiation'] == []
    with np.load(root/'chi-0.0.npz') as first, np.load(root/'chi-0.7.npz') as second:
        np.testing.assert_array_equal(first['trajectory'], second['trajectory'])
        np.testing.assert_array_equal(first['times'], np.arange(0.,13.,2.))
        np.testing.assert_array_equal(first['trajectory'][0], initial)
    assert study.assess(root)['completed_pairs'] == 0


def test_matched_comparison_refuses_changed_start_or_observation_times(tmp_path):
    root = tmp_path/'reference'; p, initial = fixture(root)
    study.run(root)
    folder = root/'moving'/'chi-0.0'; folder.mkdir()
    with np.load(root/'chi-0.0.npz') as z:
        history = [dict(elapsed=float(t), chemistry=s.tolist()) for t,s in zip(z['times'],z['trajectory'])]
    def save():
        write_json(folder/'history.json', history)
        write_json(folder/'result.json', dict(protocol_sha256=p['parent_protocol_sha256'],
            history_sha256=digest(folder/'history.json'), late_min_log_activator_sd=0., persistent_contrast=False))
    save()
    assert study.assess(root)['completed_pairs'] == 1
    history[0]['chemistry'][0][0] += .01; save()
    with pytest.raises(ValueError, match='initial chemistry'): study.assess(root)
    history[0]['chemistry'] = initial.tolist(); history[-1]['elapsed'] += .1; save()
    with pytest.raises(ValueError, match='horizon mismatch'): study.assess(root)
