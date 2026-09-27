"""Independent 3D and time-integration checks of the L-prism spectral study."""

import json
import numpy as np
import pytest
from scipy.integrate import solve_ivp

from embryo.irregular_signaling import (prism_spectrum, stability_report, select_probes,
                                         verify_growth, run_benchmark)
from embryo.signaling import gm_reaction, gm_jacobian


@pytest.mark.parametrize("grading", [0., .35])
def test_factor_spectrum_matches_entire_independent_dense_3d_eigenproblem(grading):
    spectrum = prism_spectrum(4, grading=grading)
    actual = np.linalg.eigvalsh(spectrum.transport.symmetric.toarray())
    np.testing.assert_allclose(spectrum.eigenvalues, actual, atol=4e-13)
    assert len(actual) == 48
    assert np.count_nonzero(np.abs(actual) < 1e-10) == 1
    assert spectrum.factorization_relative_residual < 1e-14
    # All physical eigenvectors, including degeneracies, are orthonormal in
    # the volume inner product. Euclidean normalization alone is insufficient.
    modes = np.column_stack([spectrum.mode(i) for i in range(len(actual))])
    weights = spectrum.transport.volumes / spectrum.transport.volumes.sum()
    np.testing.assert_allclose(modes.T @ (weights[:, None] * modes), np.eye(48), atol=2e-13)
    np.testing.assert_allclose(spectrum.transport.delta @ modes,
                               -modes * spectrum.eigenvalues, atol=3e-12)
    np.testing.assert_allclose(modes[:, 0], 1, atol=1e-13)
    np.testing.assert_allclose(weights @ modes[:, 1:], 0, atol=1e-13)


def test_spatial_spectrum_scales_as_inverse_length_squared():
    unit, scaled = prism_spectrum(6), prism_spectrum(6, length=2.5)
    np.testing.assert_allclose(scaled.eigenvalues, unit.eigenvalues / 2.5**2, atol=3e-13)
    # The z factor also has an analytic continuum reference, but graded
    # finite-volume eigenvalues are not replaced by a uniform-grid formula.
    assert abs(prism_spectrum(16).z_values[1] - np.pi**2) < abs(unit.z_values[1] - np.pi**2)


def test_mode_blocks_match_full_linearized_chemical_jacobian():
    spectrum = prism_spectrum(4)
    delta = spectrum.transport.delta.toarray()
    identity = np.eye(len(delta))
    jacobian = np.block([[identity + .02 * delta, -identity],
                         [4 * identity, -2 * identity + .4 * delta]])
    independent = np.linalg.eigvals(jacobian)
    blocks = np.concatenate([np.linalg.eigvals(gm_jacobian() - value * np.diag([.02, .4]))
                              for value in spectrum.eigenvalues])
    np.testing.assert_allclose(np.sort(independent.real), np.sort(blocks.real), atol=5e-12)
    np.testing.assert_allclose(np.sort(abs(independent.imag)), np.sort(abs(blocks.imag)), atol=5e-12)
    report = stability_report(spectrum)
    assert report["local_stable"]
    assert report["unstable_mode_count"] == np.count_nonzero(independent.real > 1e-10)
    lower, upper = report["stationary_instability_band"]
    assert report["unstable_indices"] == np.flatnonzero((spectrum.eigenvalues > lower) &
                                                        (spectrum.eigenvalues < upper)).tolist()


def test_equal_diffusivities_are_a_stable_control():
    spectrum = prism_spectrum(4)
    report = stability_report(spectrum, da=.02, dh=.02)
    assert report["local_stable"]
    assert report["unstable_mode_count"] == 0
    assert not report["diffusion_driven_instability"]
    np.testing.assert_allclose(report["growth_rates"], -.5 - .02 * spectrum.eigenvalues, atol=1e-12)


@pytest.mark.parametrize("label", ["growing", "decaying"])
def test_tiny_nonlinear_mode_matches_growth_and_independent_ode_solution(label):
    spectrum = prism_spectrum(4)
    selected = dict(select_probes(spectrum, stability_report(spectrum)))[label]
    report, fields = verify_growth(spectrum, selected, duration=.1, dt=.005)
    assert (report["predicted_growth"] > 0) == (label == "growing")
    assert report["growth_error"] < 1e-4
    assert report["relative_linear_field_error"] < 1e-3
    assert report["relative_eigen_residual"] < 1e-10
    assert report["volume_weighted_mode_rms"] == pytest.approx(1.)
    assert abs(report["volume_weighted_mode_mean"]) < 1e-12
    chemical = np.array(report["chemical_eigenvector"])
    initial = np.concatenate([1 + 1e-5 * chemical[i] * fields["spatial_mode"] for i in range(2)])
    count = len(spectrum.eigenvalues)

    def rhs(t, state):
        a, h = state[:count], state[count:]
        fa, fh = gm_reaction(a, h)
        return np.r_[fa + .02 * (spectrum.transport.delta @ a),
                     fh + .4 * (spectrum.transport.delta @ h)]

    result = solve_ivp(rhs, [0, .1], initial, method="DOP853", rtol=1e-12, atol=1e-14)
    assert result.success
    np.testing.assert_allclose(np.r_[fields["activator_final"], fields["inhibitor_final"]],
                               result.y[:, -1], rtol=0, atol=5e-10)


def test_refinement_tracks_complete_counts_and_spectrum_changes():
    spectra = [prism_spectrum(n) for n in (4, 8, 16)]
    reports = [stability_report(spectrum) for spectrum in spectra]
    assert [report["unstable_mode_count"] for report in reports] == [9, 7, 4]
    assert spectra[0].eigenvalues[1] < spectra[1].eigenvalues[1] < spectra[2].eigenvalues[1]
    # The true first mode is absent from the even-cosine reference family.
    assert spectra[-1].eigenvalues[1] < np.pi**2


def test_cli_artifacts_include_failed_coarse_refinement_criteria(tmp_path):
    output = tmp_path / "spectra"
    report = run_benchmark(output, resolutions=(4, 8), duration=.05)
    assert json.loads((output / "analysis.json").read_text()) == report
    assert report["checks"]["measured_growth_matches_prediction"]
    assert not report["checks"]["finest_unstable_counts_agree"]
    assert not report["checks"]["finest_low_spectrum_change_below_3_percent"]
    for name in ("spectrum.png", "modes.png", "modes.npz"):
        assert (output / name).stat().st_size > 0
    with np.load(output / "modes.npz", allow_pickle=False) as data:
        assert all(np.isfinite(data[key]).all() for key in data.files)
        assert data["n8_growing_spatial_mode"].shape == (384,)
    original = (output / "analysis.json").read_bytes()
    with pytest.raises(FileExistsError):
        run_benchmark(output, resolutions=(4, 8))
    assert (output / "analysis.json").read_bytes() == original


@pytest.mark.parametrize("kwargs", [
    {"resolutions": (4,)}, {"resolutions": (8, 4)}, {"resolutions": (4, 4)},
    {"resolutions": (4, 5)}, {"resolutions": (4, 8.)}, {"resolutions": (4, 64)},
    {"beta": 1.}, {"beta": float("nan")}, {"da": 0.}, {"dh": -1.},
    {"duration": 0.}, {"dt": True}, {"amplitude": .01},
    {"length": 0.}, {"grading": .9},
])
def test_invalid_parameters_leave_no_output(tmp_path, kwargs):
    output = tmp_path / "invalid"
    with pytest.raises(ValueError):
        run_benchmark(output, **{"resolutions": (4, 8), **kwargs})
    assert not output.exists()
