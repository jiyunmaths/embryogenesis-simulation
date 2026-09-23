import numpy as np
from scipy.ndimage import laplace

from embryo.polarity import exposure_cue, tension_field, flux_divergence
from embryo import Config, Simulation


def test_free_sphere_does_not_have_a_grid_selected_polarity():
    sim = Simulation(Config(grid=20, interface_width=.14))
    np.testing.assert_array_equal(exposure_cue(sim.phi, sim.dx), 0)
    for _ in range(10):
        sim.step()
    np.testing.assert_array_equal(sim.polarity, 0)


def test_contact_directs_apical_cue_away_from_neighbor():
    sim = Simulation(Config(grid=24, interface_width=.12))
    fields = []
    for x in (-.4, .4):
        r = np.sqrt(np.sum((sim.xyz - np.array([x, 0, 0])[:, None, None, None])**2, axis=0))
        fields.append(.5 * (1 - np.tanh((r - .55) / .12)))
    cues = exposure_cue(np.array(fields), sim.dx)
    assert cues[0, 0] < -.01 and cues[1, 0] > .01
    np.testing.assert_allclose(cues[:, 1:], 0, atol=1e-7)


def test_flux_includes_coefficient_gradient_and_has_no_boundary_source():
    field = np.random.default_rng(1).uniform(size=(8, 8, 8))
    constant = np.full_like(field, 1.3)
    np.testing.assert_allclose(flux_divergence(field, constant, .1), 1.3 * laplace(field, mode="nearest") / .1**2, atol=1e-12)
    coefficient = np.linspace(.5, 1.5, 8)[:, None, None] * np.ones_like(field)
    force = flux_divergence(field, coefficient, .1)
    assert abs(force.sum()) < 1e-10
    assert not np.allclose(force, coefficient * laplace(field, mode="nearest") / .1**2)


def test_polarity_changes_cortical_mechanics_and_apical_basal_tension():
    sim = Simulation(Config(grid=20, interface_width=.14, max_cells=1))
    gamma = tension_field(sim.xyz, np.array([0, 0, .8]), 1., .35, .14)
    assert gamma[:, :, -1].mean() < gamma[:, :, 0].mean()
    initial = sim.phi.copy()
    sim.mechanical_step()
    isotropic = sim.phi.copy()
    sim.phi = initial
    sim.polarity[0] = [0, 0, .8]
    sim.mechanical_step()
    assert np.max(abs(sim.phi - isotropic)) > 1e-5
