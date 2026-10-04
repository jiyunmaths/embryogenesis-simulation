import numpy as np
from scipy.integrate import solve_ivp

from embryo.model import Config
from embryo.mechanics_float64_reference import (
    coefficients, flux_divergence, force_from_coefficients, frozen_gamma_energy, geometry, rhs)


def test_no_flux_operator_matches_face_balance_and_conserves_integral():
    rng = np.random.default_rng(41)
    u, gamma = rng.uniform(size=(5, 4, 3)), rng.uniform(.2, 1.5, size=(5, 4, 3))
    expected = np.zeros_like(u); dx = .17
    for index in np.ndindex(u.shape):
        for axis in range(3):
            if index[axis]+1 < u.shape[axis]:
                neighbor = list(index); neighbor[axis] += 1; neighbor = tuple(neighbor)
                flow = (gamma[index]+gamma[neighbor])*.5*(u[neighbor]-u[index])/dx**2
                expected[index] += flow; expected[neighbor] -= flow
    actual = flux_divergence(u, gamma, dx)
    np.testing.assert_allclose(actual, expected, rtol=2e-14, atol=7e-14)
    assert abs(actual.sum()) < 1e-12
    np.testing.assert_allclose(flux_divergence(np.ones_like(u), gamma, dx), 0., atol=1e-13)


def test_force_is_negative_energy_gradient_with_frozen_gamma():
    rng = np.random.default_rng(3)
    phi = rng.uniform(.1, .85, size=(2, 5, 5, 5))
    dx, width, stiffness, repulsion = .2, .08, 1.7, .6
    gamma = rng.uniform(.5, 1.2, phi.shape)
    target = np.array([.7, .85])
    adhesion = np.array([[0., .24], [.24, 0.]])
    q = phi*phi*(1.-phi)**2
    attraction = (adhesion@q.reshape(2, -1)).reshape(phi.shape)
    volume, _ = geometry(phi, dx, .5)
    vf = stiffness*(target-volume)/target
    squares = phi*phi; excluded = squares.sum(0)[None]-squares
    force = force_from_coefficients(phi, gamma, vf, excluded, attraction, dx, width, repulsion)
    direction = rng.normal(size=phi.shape); direction /= np.linalg.norm(direction)
    energy = lambda field:frozen_gamma_energy(field, gamma, adhesion, target, dx, width, stiffness, repulsion)
    eps = 1e-5
    derivative = (energy(phi+eps*direction)-energy(phi-eps*direction))/(2*eps)
    np.testing.assert_allclose(derivative, -np.sum(force*direction)*dx**3, rtol=1e-7, atol=1e-10)


def test_geometry_and_directional_tension_keep_double_precision():
    c = Config(grid=4, extent=1., interface_width=.12, polarity_tension=.7)
    phi = np.full((1, 4, 4, 4), .5, dtype=np.float32)
    volume, center = geometry(phi, .5, 1.)
    np.testing.assert_allclose(volume, 4.)
    np.testing.assert_allclose(center, 0., atol=1e-16)
    gamma, *_ = coefficients(phi, np.ones(1), np.array([[1., 0., 0.]]), volume, c)
    assert gamma.dtype == np.float64
    np.testing.assert_allclose(gamma[0, 0]+gamma[0, -1], 2*c.surface_tension, atol=1e-15)
    assert gamma[0, 0, 0, 0] > gamma[0, -1, 0, 0]


def test_independent_euler_update_converges_to_scalar_ode_reference():
    c = Config(grid=4, extent=1., interface_width=.12, surface_tension=1.,
               polarity_tension=0., volume_stiffness=0., repulsion=0., adhesion=0.)
    initial = .35; duration = .12
    truth = solve_ivp(lambda t, y:-(2*y-6*y*y+4*y*y*y), (0., duration), [initial],
                      method='DOP853', rtol=2e-13, atol=2e-14).y[0, -1]
    errors = []
    for dt in (.00375, .001875, .0009375):
        phi = np.full((1, 4, 4, 4), initial, dtype=np.float64)
        for _ in range(round(duration/dt)):
            phi += dt*rhs(phi, np.ones(1), np.zeros((1, 3)), np.ones(1), c)
        errors.append(abs(phi[0, 0, 0, 0]-truth))
        assert phi.dtype == np.float64
    np.testing.assert_allclose(np.array(errors[:-1])/errors[1:], 2., rtol=.01)
