import numpy as np
import pytest

from embryo.signaling import (normalized_graph, stability, integrate, mode_transfer,
                              cleavage_prolongation, gm_reaction, gm_jacobian)
from embryo.graph_analysis import cycle, verify_mode


def test_normalization_preserves_uniform_state_on_irregular_graph_and_isolate():
    w = np.array([[0, 2, 1, 0], [2, 0, .2, 0], [1, .2, 0, 0], [0, 0, 0, 0]])
    graph = normalized_graph(w)
    np.testing.assert_allclose(graph.delta @ np.ones(4), 0, atol=1e-15)
    np.testing.assert_allclose(graph.degree @ graph.delta, 0, atol=1e-15)
    np.testing.assert_allclose(np.sort(np.linalg.eigvals(-graph.delta)), graph.eigenvalues, atol=1e-12)
    assert stability(graph)["components"] == 2
    assert stability(graph)["spectral_gap"] == 0
    a, h = integrate(np.ones(4), np.ones(4), graph, 1.)
    np.testing.assert_allclose(a, 1, atol=1e-14)
    np.testing.assert_allclose(h, 1, atol=1e-14)


def test_discrete_gap_can_miss_continuous_unstable_band():
    small, large = stability(cycle(4)), stability(cycle(8))
    assert small["local_stable"] and large["local_stable"]
    np.testing.assert_allclose(small["continuous_lambda_band"], [.1298437881283576, .7701562118716424])
    assert not small["diffusion_driven_instability"]
    assert large["diffusion_driven_instability"]
    assert len(large["unstable_modes"]) == 2
    report = stability(cycle(16))
    assert report["growth_rates"][1] < 0  # Refined fundamental is no longer unstable.
    assert len(report["unstable_modes"]) == 4


def test_equal_diffusivities_do_not_destabilize_stable_local_kinetics():
    assert not stability(cycle(16), da=1., dh=1.)["diffusion_driven_instability"]


def test_locally_unstable_kinetics_are_not_labeled_turing():
    report = stability(cycle(8), beta=.5)
    assert not report["local_stable"]
    assert not report["diffusion_driven_instability"]


def test_cleavage_preserves_target_weighted_signals_and_inherits_polarity():
    from embryo import Config, Simulation
    sim = Simulation(Config(grid=20, interface_width=.14))
    sim.activator[:] = 1.7
    sim.inhibitor[:] = .8
    sim.polarity[0] = [.2, .1, 0]
    before = [float(sim.target @ sim.activator), float(sim.target @ sim.inhibitor)]
    sim.divide(0, [1, 0, 0])
    daughters, neck, overlap, fraction = sim._cleavage_fields(0, sim.divisions[0])
    # Unit test of inheritance only; the public scheduler enforces neck gating.
    sim._complete_division(0, daughters, fraction, neck, overlap)
    np.testing.assert_allclose([sim.target @ sim.activator, sim.target @ sim.inhibitor], before)
    np.testing.assert_allclose(sim.polarity, [[.2, .1, 0], [.2, .1, 0]])
    assert len(sim.graph_events) == 1


@pytest.mark.parametrize("n", [4, 8])
def test_measured_linear_mode_matches_graph_prediction(n):
    result = verify_mode(cycle(n))
    assert result["absolute_error"] < 2e-5


def test_jacobian_matches_kinetics():
    eps = 1e-6
    numerical = np.column_stack([(np.array(gm_reaction(*(np.ones(2) + eps * v)))
                                  - np.array(gm_reaction(*(np.ones(2) - eps * v)))) / (2 * eps)
                                 for v in np.eye(2)])
    np.testing.assert_allclose(numerical, gm_jacobian(), rtol=1e-6)


def test_modal_prediction_matches_full_irregular_graph_jacobian():
    from scipy.optimize import linear_sum_assignment
    graph = normalized_graph([[0, 2, .3], [2, 0, .7], [.3, .7, 0]])
    j, eye = gm_jacobian(), np.eye(3)
    full = np.block([[j[0, 0] * eye + graph.delta, j[0, 1] * eye],
                     [j[1, 0] * eye, j[1, 1] * eye + 20 * graph.delta]])
    predicted = np.concatenate([np.linalg.eigvals(j - value * np.diag([1., 20.]))
                                for value in graph.eigenvalues])
    actual = np.linalg.eigvals(full)
    distances = abs(predicted[:, None] - actual[None, :])
    rows, columns = linear_sum_assignment(distances)
    assert distances[rows, columns].max() < 1e-10


def test_degenerate_mode_transfer_is_basis_invariant():
    before, after = cycle(4), cycle(8)
    p = np.repeat(np.eye(4), 2, axis=0)
    first = mode_transfer(before, after, p)
    angle = .37
    rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    before.eigenvectors[:, 1:3] = before.eigenvectors[:, 1:3] @ rotation
    after.eigenvectors[:, 1:3] = after.eigenvectors[:, 1:3] @ rotation
    second = mode_transfer(before, after, p)
    np.testing.assert_allclose(first["energy_fraction_new_by_old"], second["energy_fraction_new_by_old"], atol=1e-12)
    np.testing.assert_allclose(np.array(first["energy_fraction_new_by_old"]).sum(axis=0), 1)


def test_prolongation_copies_mother_without_imposing_new_signal_peak():
    p = cleavage_prolongation(4, 1)
    np.testing.assert_array_equal(p @ [1, 2, 3, 4], [1, 3, 4, 2, 2])
    np.testing.assert_array_equal(p @ np.ones(4), np.ones(5))


def test_signals_stay_positive_through_nonlinear_pattern_growth():
    graph = cycle(8)
    a = 1 + .01 * np.cos(2 * np.pi * np.arange(8) / 8)
    h = np.ones(8)
    for _ in range(400):
        a, h = integrate(a, h, graph, .1)
    assert np.all(a > 0) and np.all(h > 0)
    assert a.std() > .1
