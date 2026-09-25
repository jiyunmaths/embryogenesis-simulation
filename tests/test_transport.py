"""Conservative bulk transport is distinct from normalized activity exchange."""

import numpy as np
import pytest
from scipy import linalg, sparse
from scipy.integrate import solve_ivp
from scipy.sparse.linalg import eigsh, expm_multiply

from embryo.signaling import gm_reaction, normalized_graph
from embryo.transport import cartesian_transport, conservative_transport, integrate_gm


def irregular_transport():
    conductance = np.array([[0., 2., 1.], [2., 0., 3.], [1., 3., 0.]])
    return conservative_transport(conductance, [1., 2., 4.])


def test_irregular_transport_preserves_constants_and_volume_weighted_amount():
    transport = irregular_transport()
    values = np.array([2., 1., .5])
    np.testing.assert_allclose(transport.delta @ np.ones(3), 0, atol=1e-14)
    np.testing.assert_allclose(transport.volumes @ transport.delta, 0, atol=1e-14)
    derivative = transport.delta @ values
    assert abs(transport.volumes @ derivative) < 1e-14
    np.testing.assert_allclose(derivative, [-3.5, .25, .75])
    assert sparse.isspmatrix_csr(transport.delta)
    assert sparse.isspmatrix_csr(transport.symmetric)
    assert sparse.isspmatrix_csr(transport.conductance)
    # This is the mismatch that motivates the new transport interpretation.
    normalized = normalized_graph(transport.conductance.toarray(), 0)
    assert abs(transport.volumes @ (normalized.delta @ values)) > .1



def test_roundoff_asymmetry_is_reconciled_without_mutating_inputs():
    matrix = sparse.csr_matrix([[0., 2.], [2. + 1e-13, 0.]])
    volumes = np.array([.1, .9])
    original = matrix.copy()
    transport = conservative_transport(matrix, volumes)
    assert (transport.conductance - transport.conductance.T).nnz == 0
    np.testing.assert_allclose(volumes @ transport.delta, 0, atol=1e-15)
    np.testing.assert_array_equal(matrix.toarray(), original.toarray())
    transport.volumes[0] = .2
    np.testing.assert_array_equal(volumes, [.1, .9])


def test_generalized_spectrum_and_physical_modes_match_weighted_operator():
    transport = irregular_transport()
    volumes = transport.volumes
    stiffness = -np.diag(volumes) @ transport.delta.toarray()
    expected, _ = linalg.eigh(stiffness, np.diag(volumes))
    observed, symmetric_modes = np.linalg.eigh(transport.symmetric.toarray())
    np.testing.assert_allclose(observed, expected, atol=1e-14)
    np.testing.assert_allclose(np.sort(-np.linalg.eigvals(transport.delta.toarray())), expected, atol=1e-14)
    modes = symmetric_modes / np.sqrt(volumes)[:, None]
    np.testing.assert_allclose(transport.delta @ modes, -modes * observed, atol=1e-14)
    np.testing.assert_allclose(modes.T @ (volumes[:, None] * modes), np.eye(3), atol=1e-14)
    assert expected[0] >= -1e-14
    # Physical concentration eigenvectors are not raw symmetric eigenvectors.
    np.testing.assert_allclose(np.diff(modes[:, 0]), 0, atol=1e-14)


def test_diffusion_exponential_is_positive_and_conserves_amount():
    transport = irregular_transport()
    initial = np.array([0., 2., 0.])
    final = expm_multiply(.3 * transport.delta, initial)
    assert np.all(final > 0)
    np.testing.assert_allclose(transport.volumes @ final, transport.volumes @ initial, rtol=1e-14)
    assert final.max() < initial.max()
    assert np.ptp(expm_multiply(20 * transport.delta, initial)) < 1e-13


def test_isolated_compartments_and_disconnected_components():
    conductance = sparse.csr_matrix([[0., 1., 0.], [1., 0., 0.], [0., 0., 0.]])
    transport = conservative_transport(conductance, [2., 3., 7.])
    final = expm_multiply(transport.delta, np.array([1., 2., 3.]))
    np.testing.assert_allclose(final[-1], 3, rtol=1e-15)
    assert np.count_nonzero(np.linalg.eigvalsh(transport.symmetric.toarray()) < 1e-12) == 2
    single = cartesian_transport(1)
    assert single.delta.nnz == 0
    np.testing.assert_array_equal(single.centers, [[.5, .5, .5]])


def test_uniform_cube_geometry_and_no_flux_boundaries():
    transport = cartesian_transport(4, length=2)
    assert transport.shape == (4, 4, 4)
    np.testing.assert_allclose(transport.volumes, .5**3)
    np.testing.assert_allclose(transport.conductance.data, .5)
    assert transport.volumes.sum() == 8
    np.testing.assert_allclose(transport.centers[0], [.25, .25, .25])
    np.testing.assert_allclose(transport.centers[-1], [1.75, 1.75, 1.75])
    np.testing.assert_array_equal(transport.conductance[0].indices, [1, 4, 16])
    assert transport.conductance[0, 48] == 0
    assert transport.delta[0, 0] == -12
    assert transport.delta[21, 21] == -24
    assert transport.conductance.nnz == 2 * 3 * 3 * 4**2
    np.testing.assert_allclose(transport.delta @ np.ones(64), 0)


def test_nonuniform_cartesian_areas_distances_volumes_and_flattening():
    edges = (np.array([0., .1, .4, 1.]), np.array([0., .2, 1.]), np.array([0., .3, .6, 1.]))
    transport = cartesian_transport(edges=edges)
    assert transport.shape == (3, 2, 3)
    np.testing.assert_allclose(transport.volumes.sum(), 1)
    np.testing.assert_allclose(transport.centers[5], [.05, .6, .8])
    np.testing.assert_allclose(transport.volumes[5], .1 * .8 * .4)
    np.testing.assert_allclose(transport.conductance[5, 11], .8 * .4 / .2)
    np.testing.assert_allclose(transport.conductance[5, 2], .1 * .4 / .5)
    np.testing.assert_allclose(transport.conductance[5, 4], .1 * .8 / .35)
    np.testing.assert_allclose(transport.delta @ np.ones(18), 0, atol=3e-14)
    np.testing.assert_allclose(transport.volumes @ transport.delta, 0, atol=1e-15)
    initial = np.arange(18, dtype=float) / 18
    final = expm_multiply(.015 * transport.delta, initial)
    np.testing.assert_allclose(transport.volumes @ final, transport.volumes @ initial, rtol=1e-13)


def test_uniform_3d_cosine_modes_converge_to_neumann_spectrum():
    errors = []
    frequency = np.array([1, 2, 1])
    continuum = np.pi**2 * np.sum(frequency**2)
    for n in (4, 8, 16):
        transport = cartesian_transport(n)
        mode = np.cos(np.pi * transport.centers * frequency).prod(axis=1)
        exact_discrete = 4 * n**2 * np.sum(np.sin(np.pi * frequency / (2*n))**2)
        np.testing.assert_allclose(transport.delta @ mode, -exact_discrete * mode, rtol=3e-13, atol=2e-12)
        errors.append(abs(exact_discrete - continuum))
    assert 3.5 < errors[0] / errors[1] < 4.1
    assert 3.8 < errors[1] / errors[2] < 4.1
    low = np.sort(eigsh(cartesian_transport(8).symmetric, k=5, which="SM", return_eigenvectors=False))
    first = 4 * 8**2 * np.sin(np.pi / 16)**2
    np.testing.assert_allclose(low, [0, first, first, first, 2*first], atol=2e-12)


def test_length_scaling_and_rectangular_resolution():
    unit = cartesian_transport((2, 3, 4))
    scaled = cartesian_transport((2, 3, 4), length=3)
    np.testing.assert_allclose(scaled.volumes, 27 * unit.volumes)
    np.testing.assert_allclose(scaled.conductance.toarray(), 3 * unit.conductance.toarray())
    np.testing.assert_allclose(scaled.delta.toarray(), unit.delta.toarray() / 9)
    np.testing.assert_allclose(scaled.centers, 3 * unit.centers)


def test_sparse_storage_stays_linear_in_compartment_count():
    transport = cartesian_transport(24)
    count = 24**3
    assert transport.delta.shape == (count, count)
    assert transport.delta.nnz < 7 * count
    assert transport.symmetric.nnz < 7 * count
    assert transport.conductance.nnz < 6 * count


def test_gm_homogeneous_equilibrium_and_input_immutability():
    transport = irregular_transport()
    a, h = np.ones(3), np.ones(3)
    result_a, result_h = integrate_gm(a, h, transport, .1)
    np.testing.assert_allclose(result_a, 1, atol=1e-15)
    np.testing.assert_allclose(result_h, 1, atol=1e-15)
    np.testing.assert_array_equal(a, np.ones(3))
    np.testing.assert_array_equal(h, np.ones(3))
    assert not np.shares_memory(result_a, a)
    assert not np.shares_memory(result_h, h)


def test_gm_substeps_preserve_positivity_for_strong_transport():
    transport = cartesian_transport(5)
    rng = np.random.default_rng(7)
    a, h = rng.uniform(0, 1, 125), rng.uniform(.1, 1, 125)
    original_a, original_h = a.copy(), h.copy()
    result_a, result_h = integrate_gm(a, h, transport, .3, da=.3, dh=2.)
    assert np.isfinite(result_a).all() and np.isfinite(result_h).all()
    assert np.all(result_a >= 0) and np.all(result_h > 0)
    np.testing.assert_array_equal(a, original_a)
    np.testing.assert_array_equal(h, original_h)


def test_gm_second_order_accuracy_against_independent_ode_solver():
    transport = irregular_transport()
    initial_a, initial_h = np.array([.9, 1.05, 1.02]), np.array([1.01, .99, 1.03])

    def rhs(t, values):
        a, h = values[:3], values[3:]
        fa, fh = gm_reaction(a, h)
        return np.r_[fa + .02*(transport.delta @ a), fh + .4*(transport.delta @ h)]

    reference = solve_ivp(rhs, [0, .2], np.r_[initial_a, initial_h], rtol=1e-12, atol=1e-14).y[:, -1]
    errors = []
    for dt in (.02, .01, .005):
        a, h = initial_a.copy(), initial_h.copy()
        for _ in range(round(.2/dt)):
            a, h = integrate_gm(a, h, transport, dt)
        errors.append(np.linalg.norm(np.r_[a, h] - reference))
    assert 3.8 < errors[0] / errors[1] < 4.3
    assert 3.8 < errors[1] / errors[2] < 4.3


def test_gm_zero_diffusion_is_per_compartment_reaction_only():
    transport = irregular_transport()
    a, h = integrate_gm(np.ones(3)*.9, np.ones(3)*1.1, transport, .02, da=0, dh=0)
    np.testing.assert_allclose(a, a[0], atol=1e-15)
    np.testing.assert_allclose(h, h[0], atol=1e-15)


@pytest.mark.parametrize("matrix, volumes, message", [
    (np.zeros((0, 0)), [], "nonempty square"),
    (np.zeros((2, 3)), [1, 1], "square"),
    ([0], [1], "square"),
    ([[0, 1], [1, 0]], [1], "volumes"),
    ([[0, 1], [1, 0]], [0, 1], "volumes"),
    ([[0, 1], [1, 0]], [float("inf"), 1], "volumes"),
    ([[0, -1], [-1, 0]], [1, 1], "nonnegative"),
    ([[0, float("nan")], [1, 0]], [1, 1], "finite"),
    ([[0, 1], [2, 0]], [1, 1], "symmetric"),
    ([[1, 1], [1, 0]], [1, 1], "diagonal"),
    ([[0, 1j], [1j, 0]], [1, 1], "real"),
    ([[0, 1], [1, 0]], [1e-320, 1], "nonfinite"),
])
def test_general_transport_rejects_invalid_inputs(matrix, volumes, message):
    with pytest.raises(ValueError, match=message):
        conservative_transport(matrix, volumes)


@pytest.mark.parametrize("kwargs", [
    {}, {"n": 0}, {"n": True}, {"n": 2.5}, {"n": (2, 2)},
    {"n": (2, 2, 0)}, {"n": 2, "length": 0}, {"n": 2, "length": float("nan")},
    {"edges": ([0, 1], [0, 1])},
    {"edges": ([0, 1], [0, 0], [0, 1])},
    {"edges": ([0, 1], [0, float("inf")], [0, 1])},
    {"n": 2, "edges": ([0, 1], [0, 1], [0, 1])},
])
def test_cartesian_transport_rejects_invalid_geometry(kwargs):
    with pytest.raises(ValueError):
        cartesian_transport(**kwargs)


@pytest.mark.parametrize("a,h,parameters", [
    ([1, 1], [1, 1], {}), ([1, -1, 1], [1, 1, 1], {}),
    ([1, 1, 1], [1, 0, 1], {}), ([1, np.nan, 1], [1, 1, 1], {}),
    ([1, 1, 1], [1, np.inf, 1], {}),
    ([1, 1, 1], [1, 1, 1], {"dt": 0}),
    ([1, 1, 1], [1, 1, 1], {"beta": -1}),
    ([1, 1, 1], [1, 1, 1], {"da": -1}),
    ([1, 1, 1], [1, 1, 1], {"dh": np.inf}),
])
def test_gm_rejects_invalid_states_and_parameters(a, h, parameters):
    with pytest.raises(ValueError):
        integrate_gm(a, h, irregular_transport(), **{"dt": .1, **parameters})
