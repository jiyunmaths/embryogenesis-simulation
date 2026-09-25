import numpy as np
import pytest

from embryo import Config
from embryo.graph_analysis import cycle
from embryo.signaling import gm_jacobian
from embryo.timescales import (apply_inheritance, capture, load_capture, replay,
                               frozen, first_crossing, matched_followup)


def test_partition_replay_reuses_relative_draw_and_preserves_weighted_mean():
    # Unequal daughter volumes: .3 and .7. Original parent was 2, new parent is 5.
    event = {"prolongation": [[1], [1]], "activator_before": [2.],
             "activator_partition_jump": [.14, -.06]}
    inherited = apply_inheritance(np.array([5.]), event, "activator")
    np.testing.assert_allclose(inherited, [5.35, 4.85])
    np.testing.assert_allclose(np.array([.3, .7]) @ inherited, 5.)


@pytest.fixture(scope="module")
def recorded(tmp_path_factory):
    path = tmp_path_factory.mktemp("timescales") / "capture"
    capture(Config(grid=20, interface_width=.14, division_interval=.12,
                   max_cells=4, steps=300), path)
    return path, load_capture(path)


def test_native_replay_reconstructs_full_signals_across_real_cleavages(recorded):
    _, (config, trace, events) = recorded
    assert len(events) == 3
    assert trace["counts"][-1] == 4
    result = replay(config, trace, events)
    for key in ("activator", "inhibitor", "time"):
        np.testing.assert_allclose(result[key], trace[key], atol=1e-12, rtol=0)


def test_stretch_changes_precleavage_time_but_keeps_equal_postcleavage_followup(recorded):
    _, (config, trace, events) = recorded
    result = replay(config, trace, events, 2.)
    last = round(events[-1]["time"] / config.dt)
    np.testing.assert_allclose(result["time"][:last+1], 2 * trace["time"][:last+1], atol=1e-12)
    np.testing.assert_allclose(result["time"][last:] - result["time"][last],
                               trace["time"][last:] - trace["time"][last], atol=1e-12)
    np.testing.assert_array_equal(result["counts"], trace["counts"])
    assert np.all(result["inhibitor"][result["inhibitor"] != 0] > 0)


def test_capture_cannot_overwrite_existing_experiment(recorded):
    path, (config, _, _) = recorded
    with pytest.raises(FileExistsError):
        capture(config, path)


def test_frozen_linear_solution_matches_analytic_eigenmode():
    graph = cycle(8)
    lam = graph.eigenvalues[1]
    values, vectors = np.linalg.eig(gm_jacobian() - lam * np.diag([1., 20.]))
    index = np.argmax(values.real)
    chemical = vectors[:, index].real
    spatial = np.cos(2 * np.pi * np.arange(8) / 8)
    a, b = 1 + 1e-6 * chemical[0] * spatial, 1 + 1e-6 * chemical[1] * spatial
    _, result = frozen(Config(), graph.weights, a, b, duration=2.)
    expected = (a - 1).std() * np.exp(values[index].real * result["time"])
    np.testing.assert_allclose(result["linear_activator_deviation"].std(axis=1), expected, rtol=1e-9)
    np.testing.assert_allclose(result["activator"].std(axis=1), expected, rtol=2e-5)


def test_onset_reports_absence_and_initial_crossing_explicitly():
    assert first_crossing([0, 1, 2], [0, .02, .2], .1) == 2
    assert first_crossing([0, 1], [.2, .3], .1) == 0
    assert first_crossing([0, 1], [0, .02], .1) is None


def test_matched_followup_compares_age_since_final_cleavage_not_absolute_time():
    first = {"time": np.array([1., 2., 3.]), "counts": np.array([16, 16, 16]),
             "activator": np.array([[0., 2.] * 8, [0., 4.] * 8, [0., 6.] * 8])}
    second = {"time": np.array([2., 3.]), "counts": np.array([16, 16]),
              "activator": np.array([[0., 2.] * 8, [0., 4.] * 8])}
    result = matched_followup({"first": first, "second": second})
    assert result["time_since_16_cells"] == 1.
    assert result["activator_std"] == {"first": 2., "second": 2.}


@pytest.mark.parametrize("factor", [0, -1, float("nan")])
def test_invalid_replay_timescale_is_rejected(recorded, factor):
    _, (config, trace, events) = recorded
    with pytest.raises(ValueError, match="finite and positive"):
        replay(config, trace, events, factor)
