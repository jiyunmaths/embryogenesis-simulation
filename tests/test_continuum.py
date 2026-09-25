"""Analytic and cross-model checks for the fixed-domain continuum bridge."""

import json
import numpy as np
import pytest

from embryo.continuum import (continuum_eigenvalue, continuum_supported_modes,
                              cosine_average, diffusion_check, fv_eigenvalue,
                              normalized_operator, run_benchmark, supported_modes,
                              verify_growth)
from embryo.signaling import normalized_graph
from embryo.transport import cartesian_transport


def test_full_small_cube_spectrum_matches_analytic_neumann_modes():
    transport = cartesian_transport(4, length=2.)
    actual = np.linalg.eigvalsh(transport.symmetric.toarray())
    exact = sorted(fv_eigenvalue(4, (i, j, k), length=2.)
                   for i in range(4) for j in range(4) for k in range(4))
    np.testing.assert_allclose(actual, exact, atol=1e-12)
    assert continuum_eigenvalue((1, 1, 0), 2.) == pytest.approx(np.pi**2 / 2)


def test_prescribed_cell_averages_preserve_mean_and_resolve_operator_mode():
    transport = cartesian_transport(7, length=1.7)
    mode = (1, 2, 0)
    field = cosine_average(7, mode)
    assert abs(transport.volumes @ field) < 1e-14
    np.testing.assert_allclose(transport.delta @ field,
                               -fv_eigenvalue(7, mode, 1.7) * field, atol=1e-12)
    np.testing.assert_allclose(cosine_average(7, (0, 0, 0)), 1)
    with pytest.raises(ValueError, match="resolved"):
        cosine_average(7, (7, 0, 0))


def test_sparse_normalized_control_matches_existing_graph_model():
    transport = cartesian_transport(4)
    sparse, degree = normalized_operator(transport)
    current = normalized_graph(transport.conductance.toarray(), relative_cutoff=0.)
    np.testing.assert_allclose(sparse.toarray(), current.delta, atol=1e-15)
    np.testing.assert_allclose(degree, current.degree)


def test_diffusion_converges_while_fixed_exchange_changes_physical_process():
    rows = [diffusion_check(cartesian_transport(n), n, 1., .02, 1., 4)[0]
            for n in (4, 8, 16)]
    for coarse, fine in zip(rows, rows[1:]):
        assert coarse["conservative_l2_error"] / fine["conservative_l2_error"] > 3.5
        assert fine["fixed_rate_l2_error"] > coarse["fixed_rate_l2_error"]
        assert fine["interior_effective_diffusivity"] == pytest.approx(coarse["interior_effective_diffusivity"] / 4)
    for row in rows:
        assert row["conservative_relative_mass_drift"] < 1e-12
        assert row["fixed_rate_relative_mass_drift"] > 1e-5
        assert row["fixed_rate_relative_degree_sum_drift"] < 1e-12
        assert row["capacity_factorization_max_error"] < 1e-12


def test_refinement_corrects_spurious_near_edge_instability():
    exact = continuum_supported_modes(1., 2., .02, .4)
    assert len(exact) == 7
    assert supported_modes(4, 1., 2., .02, .4)["unstable_mode_count"] == 10
    fine = supported_modes(16, 1., 2., .02, .4)
    assert fine["unstable_mode_count"] == 7
    assert {tuple(item["mode"]) for item in fine["unstable_modes"]} == {tuple(item["mode"]) for item in exact}


@pytest.mark.parametrize("mode, sign", [((1, 1, 0), 1), ((2, 2, 0), -1)])
def test_nonlinear_solver_grows_and_decays_at_predicted_discrete_rate(mode, sign):
    report = verify_growth(cartesian_transport(6), 6, mode, duration=.1)
    assert sign * report["measured_growth"] > 0
    assert report["growth_error"] < 1e-4
    assert min(report["minimum_activator"], report["minimum_inhibitor"]) > 0


def test_cli_report_retains_failed_resolution_criterion_and_refuses_overwrite(tmp_path):
    output = tmp_path / "bridge"
    report = run_benchmark(output, resolutions=[4, 8], growth_time=.04)
    saved = json.loads((output / "analysis.json").read_text())
    assert saved == report
    assert not saved["checks"]["finest_unstable_count_matches_continuum"]
    assert saved["checks"]["conservative_mass_drift_below_1e_10"]
    for name in ("convergence.png", "diffusion.png", "diffusion_slices.npz"):
        assert (output / name).stat().st_size > 0
    before = (output / "analysis.json").read_bytes()
    with pytest.raises(FileExistsError):
        run_benchmark(output, resolutions=[4, 8])
    assert (output / "analysis.json").read_bytes() == before


@pytest.mark.parametrize("resolutions", [[8, 4], [4, 4], [2, 4], [4], [4, 8.0]])
def test_bad_refinement_sequences_do_not_create_output(tmp_path, resolutions):
    output = tmp_path / "invalid"
    with pytest.raises(ValueError, match="resolutions"):
        run_benchmark(output, resolutions=resolutions)
    assert not output.exists()
