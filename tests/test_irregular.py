"""Independent manufactured-solution checks on a graded, nonconvex 3D domain."""

import json

import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss

from embryo.irregular import (arm_pulse_check, cosine_average, diffusion_check, graded_edges,
                              l_prism_mask, run_benchmark)
from embryo.transport import masked_cartesian_transport


def quadrature_average(edges, mode, mask):
    """Integrate cellwise cosines with Gauss quadrature, not antiderivatives."""
    nodes, weights = leggauss(12)
    length = edges[0][-1]
    averages = []
    for axis, frequency in zip(edges, mode):
        low, high = axis[:-1, None], axis[1:, None]
        positions = .5 * (low + high) + .5 * (high - low) * nodes
        averages.append(.5 * (np.cos(np.pi * frequency * positions / length) @ weights))
    values = (averages[0][:, None, None] * averages[1][None, :, None]
              * averages[2][None, None, :])
    return values[mask]


def reference_fields(edges, mask, diffusivity, duration):
    """An independently integrated exact no-flux solution on the L prism."""
    initial = np.ones(np.count_nonzero(mask))
    exact = initial.copy()
    for amplitude, mode in ((.10, (2, 0, 0)), (.07, (0, 2, 1)), (.04, (2, 2, 2))):
        average = quadrature_average(edges, mode, mask)
        decay = np.exp(-diffusivity * duration * np.pi**2
                       * sum(frequency**2 for frequency in mode) / edges[0][-1]**2)
        initial += amplitude * average
        exact += amplitude * decay * average
    return initial, exact


def test_graded_mesh_is_nested_and_preserves_exact_domain_corners():
    coarse, fine = graded_edges(8), graded_edges(16)
    for low, high in zip(coarse, fine):
        assert len(low) == 9 and len(high) == 17
        assert low[0] == 0. and low[-1] == 1. and low[4] == .5
        assert np.all(np.diff(low) > 0)
        assert np.ptp(np.diff(low)) > .001
        np.testing.assert_allclose(high[::2], low, rtol=0, atol=2e-15)
    assert not np.allclose(coarse[0], coarse[1])
    assert not np.allclose(coarse[1], coarse[2])
    for uniform in graded_edges(8, grading=0.):
        np.testing.assert_allclose(uniform, np.linspace(0, 1, 9), rtol=0, atol=2e-15)
    for unit, scaled in zip(coarse, graded_edges(8, length=2.3)):
        np.testing.assert_allclose(scaled, 2.3 * unit, rtol=2e-15, atol=2e-15)


def test_l_prism_mask_selects_three_quadrants_in_c_order():
    mask = l_prism_mask(8)
    assert mask.shape == (8, 8, 8) and mask.dtype == np.bool_
    assert np.count_nonzero(mask) == 3 * 8**3 // 4
    assert mask[:4, :, :].all() and mask[:, :4, :].all()
    assert not mask[4:, 4:, :].any()
    edges = graded_edges(8, length=2.)
    transport = masked_cartesian_transport(edges=edges, mask=mask)
    assert transport.volumes.sum() == pytest.approx(.75 * 2.**3, abs=2e-14)
    expected_centers = np.stack(np.meshgrid(
        *[.5 * (axis[:-1] + axis[1:]) for axis in edges], indexing="ij"), axis=-1)[mask]
    np.testing.assert_allclose(transport.centers, expected_centers, atol=1e-15)


@pytest.mark.parametrize("mode", [(0, 0, 0), (2, 0, 0), (0, 2, 1), (2, 4, 3)])
def test_cosine_cell_averages_match_independent_quadrature(mode):
    edges, mask = graded_edges(6, length=1.7), l_prism_mask(6)
    result = cosine_average(edges, mode, mask)
    assert result.shape == (np.count_nonzero(mask),)
    np.testing.assert_allclose(result, quadrature_average(edges, mode, mask),
                               rtol=2e-13, atol=2e-14)
    unit = cosine_average(graded_edges(6), mode, mask)
    np.testing.assert_allclose(result, unit, rtol=2e-13, atol=2e-14)


@pytest.mark.parametrize("mode", [(1, 0, 0), (0, 1, 0), (-2, 0, 0), (2, 0, -1),
                                  (2., 0, 0), (2, 0), (True, 0, 0)])
def test_reference_modes_reject_non_neumann_or_noninteger_indices(mode):
    with pytest.raises(ValueError):
        cosine_average(graded_edges(4), mode, l_prism_mask(4))


def test_diffusion_report_uses_exact_cell_averages_and_volume_weighted_error():
    edges, mask = graded_edges(8, length=1.7), l_prism_mask(8)
    transport = masked_cartesian_transport(edges=edges, mask=mask)
    diffusivity, duration = .03, .7
    report, fields = diffusion_check(transport, edges, mask, diffusivity=diffusivity,
                                     duration=duration)
    initial, exact = reference_fields(edges, mask, diffusivity, duration)
    np.testing.assert_allclose(fields["initial"], initial, atol=2e-14)
    np.testing.assert_allclose(fields["exact"], exact, atol=2e-14)
    assert fields["numerical"].shape == initial.shape
    np.testing.assert_allclose(fields["error"], fields["numerical"] - exact, atol=2e-14)
    weighted_error = np.sqrt(transport.volumes @ (fields["numerical"] - exact)**2
                             / transport.volumes.sum())
    assert report["l2_error"] == pytest.approx(weighted_error, abs=2e-14)
    assert report["initial_mass"] == pytest.approx(transport.volumes @ initial, abs=2e-14)
    assert report["final_mass"] == pytest.approx(transport.volumes @ fields["numerical"], abs=2e-14)
    assert report["minimum_concentration"] == pytest.approx(fields["numerical"].min())
    # Positive reciprocal diffusion obeys a discrete maximum principle.
    assert fields["numerical"].min() >= initial.min() - 2e-13
    assert fields["numerical"].max() <= initial.max() + 2e-13


@pytest.mark.parametrize("grading", [0., .35])
def test_l_prism_refinement_converges_and_preserves_amount(grading):
    errors = []
    for n in (4, 8, 16):
        edges, mask = graded_edges(n, grading=grading), l_prism_mask(n)
        transport = masked_cartesian_transport(edges=edges, mask=mask)
        report, _ = diffusion_check(transport, edges, mask)
        errors.append(report["l2_error"])
        assert report["relative_mass_drift"] < 1e-11
        assert report["minimum_concentration"] > 0
        assert report["initial_mass"] == pytest.approx(.75, abs=2e-14)
    orders = np.log2(np.asarray(errors[:-1]) / errors[1:])
    # The coarsest graded mesh is preasymptotic: require monotone error
    # reduction throughout and the predeclared rate only on the finest pair.
    assert np.all(np.diff(errors) < 0)
    assert orders[-1] > 1.7, orders


def test_diffusion_length_scaling_preserves_nondimensional_solution():
    mask = l_prism_mask(6)
    unit_edges, scaled_edges = graded_edges(6), graded_edges(6, length=2.)
    unit_transport = masked_cartesian_transport(edges=unit_edges, mask=mask)
    scaled_transport = masked_cartesian_transport(edges=scaled_edges, mask=mask)
    unit_report, unit_fields = diffusion_check(unit_transport, unit_edges, mask, duration=.3)
    scaled_report, scaled_fields = diffusion_check(scaled_transport, scaled_edges, mask, duration=1.2)
    for name in ("initial", "exact", "numerical"):
        np.testing.assert_allclose(scaled_fields[name], unit_fields[name], atol=2e-13)
    assert scaled_report["initial_mass"] == pytest.approx(8 * unit_report["initial_mass"])
    assert scaled_report["l2_error"] == pytest.approx(unit_report["l2_error"], abs=2e-13)


def test_arm_pulse_crosses_junction_and_approaches_volume_weighted_equilibrium():
    report, fields = arm_pulse_check()
    transport = masked_cartesian_transport(edges=fields["edges"], mask=fields["mask"])
    volumes = transport.volumes
    destination = transport.centers[:, 0] > .5
    initial = fields["pulse_0"]
    mass = volumes @ initial
    mean = mass / volumes.sum()
    initial_rms = np.sqrt(volumes @ (initial - mean)**2 / volumes.sum())
    assert report["connected_components"] == 1
    assert mean == pytest.approx(1 / 3)
    assert report["initial_mass"] == pytest.approx(.25)
    assert report["equilibrium_concentration"] == pytest.approx(mean)
    assert report["destination_equilibrium_mass"] == pytest.approx(1 / 12)
    assert not np.any(initial[destination])
    traces = report["trace"]
    assert [entry["diffusion_time"] for entry in traces] == [0., .05, .5, 5.]
    for entry in traces:
        values = fields[f'pulse_{entry["diffusion_time"]:g}']
        assert values.shape == volumes.shape
        assert volumes @ values == pytest.approx(mass, rel=1e-11)
        assert values.min() >= -1e-13 and values.max() <= 1 + 1e-13
        expected_destination_mass = volumes[destination] @ values[destination]
        assert entry["destination_mass"] == pytest.approx(expected_destination_mass, abs=1e-14)
        assert entry["relative_mass_drift"] < 1e-11
        relative_rms = np.sqrt(volumes @ (values - mean)**2 / volumes.sum()) / initial_rms
        assert entry["relative_equilibrium_rms"] == pytest.approx(relative_rms, abs=1e-13)
        assert entry["physical_time"] == pytest.approx(entry["diffusion_time"] / .02)
    assert traces[1]["destination_mass"] > 0
    assert traces[-1]["destination_mass"] == pytest.approx(1 / 12, abs=1e-7)
    assert traces[-1]["relative_equilibrium_rms"] < 1e-6
    assert all(later["relative_equilibrium_rms"] < earlier["relative_equilibrium_rms"]
               for earlier, later in zip(traces, traces[1:]))


def test_benchmark_report_and_artifacts_round_trip_without_overwrite(tmp_path):
    output = tmp_path / "irregular"
    report = run_benchmark(output, resolutions=(4, 8))
    assert json.loads((output / "analysis.json").read_text()) == report
    assert set(report["families"]) == {"uniform_l_prism", "graded_l_prism"}
    # A short, coarse graded study must report its failed rate criterion.
    assert not report["checks"]["last_order_above_1_7_in_both_families"]
    for family in report["families"].values():
        assert [row["n"] for row in family["rows"]] == [4, 8]
        assert len(family["observed_orders"]) == 1
        assert all(row["relative_mass_drift"] < 1e-11 for row in family["rows"])
    for name in ("convergence.png", "fields.png", "mesh.png", "fields.npz"):
        assert (output / name).stat().st_size > 0
    with np.load(output / "fields.npz", allow_pickle=False) as archive:
        assert archive.files
        assert all(np.isfinite(archive[key]).all() for key in archive.files)
    original = (output / "analysis.json").read_bytes()
    with pytest.raises(FileExistsError):
        run_benchmark(output, resolutions=(4, 8))
    assert (output / "analysis.json").read_bytes() == original


@pytest.mark.parametrize("parameters", [
    {"resolutions": (4,)}, {"resolutions": (8, 4)}, {"resolutions": (4, 4)},
    {"resolutions": (3, 6)}, {"resolutions": (4, 8.)},
    {"length": 0.}, {"length": -1.}, {"length": float("nan")},
    {"length": float("inf")}, {"length": True},
    {"grading": float("nan")}, {"grading": float("inf")},
    {"diffusivity": -.02}, {"diffusivity": float("nan")},
    {"duration": 0.}, {"duration": -1.}, {"duration": float("inf")},
])
def test_invalid_benchmark_inputs_leave_no_output(tmp_path, parameters):
    output = tmp_path / "invalid"
    with pytest.raises(ValueError):
        run_benchmark(output, **{"resolutions": (4, 8), **parameters})
    assert not output.exists()


def test_roundoff_limited_convergence_is_reported_as_unresolved():
    from embryo.irregular import observed_orders
    assert observed_orders([.04, .01, .0025], [4, 8, 16]) == pytest.approx([2., 2.])
    assert observed_orders([0., 0., 1e-16], [4, 8, 16]) == [None, None]
    assert observed_orders([1e-5, 0.], [4, 8]) == [None]
    assert json.loads(json.dumps(observed_orders([0., 0.], [4, 8]), allow_nan=False)) == [None]
