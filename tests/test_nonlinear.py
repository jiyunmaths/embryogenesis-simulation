"""Independent solver, conservation, and cross-mesh initialization checks."""

import json

import numpy as np
import pytest
from scipy import sparse
from scipy.integrate import solve_ivp
from scipy.sparse.linalg import spsolve

from embryo.irregular import graded_edges, l_prism_mask
from embryo.irregular_signaling import prism_spectrum
from embryo.nonlinear import (IMEXGM, PrismResolvent, initial_fields, restrict_field,
                              field_difference, run_trajectory, run_benchmark)
from embryo.signaling import gm_reaction
from embryo.transport import conservative_transport, integrate_gm, masked_cartesian_transport


@pytest.mark.parametrize("n,length,alpha,diffusivity", [(4, 1., 1., .02), (8, 2.3, 1.5, .4)])
def test_product_resolvent_matches_independent_full_3d_solve(n, length, alpha, diffusivity):
    spectrum = prism_spectrum(n, length=length)
    rhs = np.random.default_rng(17).normal(size=len(spectrum.transport.volumes))
    dt = .07
    matrix = alpha * sparse.eye(len(rhs)) - dt * diffusivity * spectrum.transport.delta
    expected = spsolve(matrix.tocsc(), rhs)
    result = PrismResolvent(spectrum, dt, diffusivity, alpha).solve(rhs)
    np.testing.assert_allclose(result, expected, rtol=2e-12, atol=2e-13)
    np.testing.assert_allclose(matrix @ result, rhs, atol=2e-12)


def small_transport():
    return conservative_transport([[0., 1., .2], [1., 0., .5], [.2, .5, 0.]], [.2, .3, .5])


def test_sbdf2_second_order_against_independent_adaptive_ode():
    transport = small_transport()
    initial = np.array([.8, 1.3, .95, 1.1, .9, 1.05])

    def rhs(t, state):
        a, h = state[:3], state[3:]
        fa, fh = gm_reaction(a, h)
        return np.r_[fa + .02 * (transport.delta @ a), fh + .4 * (transport.delta @ h)]

    reference = solve_ivp(rhs, (0., 1.), initial, method="DOP853", rtol=1e-12, atol=1e-14)
    assert reference.success
    errors = []
    for dt in (.02, .01, .005):
        engine = IMEXGM(transport, initial[:3], initial[3:], dt)
        for _ in range(round(1 / dt)):
            engine.step()
        errors.append(np.linalg.norm(np.r_[engine.a, engine.h] - reference.y[:, -1]))
        assert engine.max_amount_balance_residual < 1e-13
    assert all(3.5 < a / b < 4.5 for a, b in zip(errors, errors[1:]))


def test_imex_matches_existing_positive_explicit_solver_on_small_prism():
    spectrum = prism_spectrum(4)
    a, h = initial_fields(spectrum.edges, spectrum.mask)
    engine = IMEXGM(spectrum.transport, a, h, .0005, spectrum=spectrum)
    for _ in range(200):
        engine.step()
        a, h = integrate_gm(a, h, spectrum.transport, .0005)
    np.testing.assert_allclose(engine.a, a, rtol=0, atol=1e-7)
    np.testing.assert_allclose(engine.h, h, rtol=0, atol=1e-7)


def test_pure_diffusion_preserves_amount_with_unequal_capacities():
    transport = small_transport()
    engine = IMEXGM(transport, [.7, 1.1, 1.4], [1.3, 1., .8], .05,
                    reaction=lambda a, h: (np.zeros_like(a), np.zeros_like(h)))
    amounts = transport.volumes @ np.column_stack((engine.a, engine.h))
    for _ in range(50):
        engine.step()
    np.testing.assert_allclose(transport.volumes @ np.column_stack((engine.a, engine.h)), amounts, atol=1e-14)
    assert engine.max_amount_balance_residual < 1e-14


def test_integer_initial_input_copied_and_homogeneous_equilibrium_preserved():
    transport = small_transport()
    a = np.ones(3, dtype=int)
    engine = IMEXGM(transport, a, a, .05)
    a[0] = 9
    for _ in range(20):
        engine.step()
    assert engine.a.dtype == np.float64
    np.testing.assert_allclose(engine.a, 1., atol=1e-14)
    np.testing.assert_allclose(engine.h, 1., atol=1e-14)


@pytest.mark.parametrize("a,h", [([0, 1, 1], [1, 1, 1]), ([1, 1, 1], [-1, 1, 1]),
                                  ([np.nan, 1, 1], [1, 1, 1]), ([1, 1, 1], [np.inf, 1, 1])])
def test_invalid_initial_concentrations_rejected(a, h):
    with pytest.raises(FloatingPointError):
        IMEXGM(small_transport(), a, h, .05)


def test_unstable_step_raises_without_clipping_or_advancing():
    engine = IMEXGM(small_transport(), np.ones(3), np.ones(3), 1.,
                    reaction=lambda a, h: (-2 * a, -2 * h))
    with pytest.raises(FloatingPointError, match="reduce dt"):
        engine.step()
    assert engine.steps == 0
    np.testing.assert_array_equal(engine.a, 1.)


@pytest.mark.parametrize("length", [1., 2.3])
def test_initial_cell_averages_match_independent_gauss_quadrature(length):
    edges, mask = graded_edges(4, length, .35), l_prism_mask(4)
    a, h = initial_fields(edges, mask)
    nodes, weights = np.polynomial.legendre.leggauss(8)
    values = []
    for i, j, k in np.argwhere(mask):
        xyz = [(edge[index] + edge[index+1]) / 2 + nodes * (edge[index+1] - edge[index]) / 2
               for edge, index in zip(edges, (i, j, k))]
        x, y, z = np.meshgrid(*xyz, indexing="ij")
        q = (np.cos(np.pi*x/length)*np.cos(np.pi*z/length) + .6*np.cos(np.pi*y/length)
             + .3*np.cos(2*np.pi*x/length)*np.cos(np.pi*y/length) - .4/np.pi)
        values.append(np.einsum("i,j,k,ijk", weights/2, weights/2, weights/2, q))
    np.testing.assert_allclose(a, 1 + .01*np.array(values), atol=2e-15)
    np.testing.assert_allclose(h, 1 - .007*np.array(values), atol=2e-15)
    transport = masked_cartesian_transport(edges=edges, mask=mask)
    assert transport.volumes @ (a - 1) == pytest.approx(0., abs=1e-15)


def test_nested_restriction_preserves_amount_and_same_continuous_initial_field():
    coarse, fine = prism_spectrum(4), prism_spectrum(8)
    initial_coarse = initial_fields(coarse.edges, coarse.mask)
    initial_fine = initial_fields(fine.edges, fine.mask)
    for c, f in zip(initial_coarse, initial_fine):
        np.testing.assert_allclose(restrict_field(f, fine.transport.volumes, 8, 4), c, atol=3e-16)
    field = np.random.default_rng(2).uniform(size=len(fine.transport.volumes))
    restricted = restrict_field(field, fine.transport.volumes, 8, 4)
    assert coarse.transport.volumes @ restricted == pytest.approx(fine.transport.volumes @ field, abs=1e-15)


def test_undefined_uniform_relative_error_is_not_reported_as_zero():
    result = field_difference(np.ones(3)*2, np.ones(3), np.ones(3), np.ones(3), np.array([1., 2., 3.]))
    assert result["absolute_joint_rms"] == pytest.approx(1.)
    assert result["relative_joint_rms"] is None
    assert result["activator_correlation"] is None


def test_trajectory_preserves_requested_times_and_reports_late_window():
    report, fields = run_trajectory(prism_spectrum(4), duration=1., dt=.025, sample_interval=.25)
    np.testing.assert_allclose(fields["time"], [0., .2, .4, .6, .75, .8, .9, 1.])
    assert report["steps"] == 40
    assert report["late_window_start"] == .75
    assert report["late_samples"] == 2
    assert report["minimum_late_std_a"] > 0
    assert report["onset_std_0_1"] is None
    assert report["minimum_a_all_steps"] > 0


def test_short_run_records_failed_scientific_checks_and_reproducible_artifacts(tmp_path):
    output = tmp_path / "short"
    report = run_benchmark(output, resolutions=(4, 8, 16), duration=1., dt=.05, sample_interval=.25)
    assert json.loads((output / "analysis.json").read_text()) == report
    assert report["checks"]["initial_fields_agree_under_restriction"]
    assert report["checks"]["discrete_reaction_amount_balance"]
    assert not report["checks"]["persistent_fine_mesh_contrast"]
    assert not report["checks"]["equal_diffusivity_control_loses_contrast"]
    for name in ("fields.npz", "convergence.png", "patterns.png"):
        assert (output / name).stat().st_size > 0
    with np.load(output / "fields.npz", allow_pickle=False) as data:
        assert data["n16_a"].shape == (8, 3072)
        assert all(np.isfinite(data[key]).all() for key in data.files)
    original = (output / "analysis.json").read_bytes()
    with pytest.raises(FileExistsError):
        run_benchmark(output)
    assert (output / "analysis.json").read_bytes() == original


def test_unresolved_perturbation_retains_report_with_failed_relative_checks(tmp_path):
    report = run_benchmark(tmp_path / "unresolved", resolutions=(4, 8, 16), duration=.2,
                           dt=.05, sample_interval=.05, amplitude=1e-18)
    assert report["observed_temporal_order"] is None
    assert not report["checks"]["fine_temporal_error_below_0_5_percent"]
    assert not report["checks"]["fine_spatial_final_error_below_5_percent"]
    assert not report["checks"]["fine_mesh_late_pattern_stationary"]


@pytest.mark.parametrize("kwargs", [
    {"resolutions": (4, 8)}, {"resolutions": (4, 8, 12)}, {"resolutions": (4, 8, 64)},
    {"resolutions": (4, 8, 8)}, {"resolutions": (4, 8, 16.)},
    {"duration": 0}, {"dt": True}, {"dt": .03}, {"grading": 1.},
    {"amplitude": .1}, {"length": -1}, {"sample_interval": 100.},
])
def test_invalid_benchmark_settings_do_not_create_output(tmp_path, kwargs):
    output = tmp_path / "bad"
    with pytest.raises(ValueError):
        run_benchmark(output, **kwargs)
    assert not output.exists()
