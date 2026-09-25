"""A fixed-domain transport bridge; not yet a continuum embryo model.

Compare sparse conservative finite volumes with fixed-rate normalized exchange
on exactly the same Cartesian compartments. Analytic Neumann modes separate
spatial discretization error from early Gierer–Meinhardt growth-rate error.
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply

from .signaling import gm_jacobian, mode_growth
from .transport import cartesian_transport, integrate_gm


def _mode(mode, n=None):
    values = np.asarray(mode)
    if values.shape != (3,) or not np.issubdtype(values.dtype, np.integer) or np.any(values < 0):
        raise ValueError("mode must contain three nonnegative integers")
    if n is not None and np.any(values >= n):
        raise ValueError("mode must be resolved: each index must be smaller than n")
    return values


def continuum_eigenvalue(mode, length=1.):
    values = _mode(mode)
    if not np.isfinite(length) or length <= 0:
        raise ValueError("length must be positive and finite")
    return float((np.pi / length)**2 * (values @ values))


def fv_eigenvalue(n, mode, length=1.):
    values = _mode(mode, n)
    continuum_eigenvalue(mode, length)  # Shared length validation.
    return float(4 * (n / length)**2 * np.sum(np.sin(values * np.pi / (2 * n))**2))


def cosine_average(n, mode):
    """Exact cell averages of a Neumann cosine product on a uniform cube.

    Scaling the cube changes its eigenvalue, but not these fractional samples.
    Matching averages, rather than point samples, isolates finite-volume error.
    """
    values = _mode(mode, n)
    coordinates = (np.arange(n) + .5) / n
    axes = [np.cos(np.pi * index * coordinates) * np.sinc(index / (2 * n)) for index in values]
    return (axes[0][:, None, None] * axes[1][None, :, None] * axes[2][None, None, :]).ravel()


def normalized_operator(transport):
    """Same random-walk rule as signaling.normalized_graph, kept sparse.

    All present faces have positive conductance; no relative contact pruning is
    applied to this ideal mesh. Degrees are conductances, not cell capacities.
    """
    degree = np.asarray(transport.conductance.sum(axis=1)).ravel()
    active = degree > 0
    inverse = np.divide(1., degree, out=np.zeros_like(degree), where=active)
    delta = diags(inverse) @ transport.conductance - diags(active.astype(float))
    return delta.tocsr(), degree


def diffusion_check(transport, n, length, diffusivity, duration, reference_n):
    modes = [(1, 1, 0), (0, 0, 2)]
    amplitudes = [.10, .06]
    spatial = [cosine_average(n, mode) for mode in modes]
    initial = 1 + sum(amplitude * vector for amplitude, vector in zip(amplitudes, spatial))
    exact = 1 + sum(amplitude * vector * np.exp(-diffusivity * continuum_eigenvalue(mode, length) * duration)
                    for amplitude, vector, mode in zip(amplitudes, spatial, modes))
    conservative = expm_multiply(diffusivity * duration * transport.delta, initial)
    normalized, degree = normalized_operator(transport)
    # One calibration at the reference mesh, using an interior six-face voxel.
    # Boundary rows have fewer faces; a single scalar is not exact even there.
    rate = 6 * diffusivity * (reference_n / length)**2
    fixed_rate = expm_multiply(rate * duration * normalized, initial)
    mass = float(transport.volumes @ initial)
    norm = lambda values: float(np.sqrt(transport.volumes @ values**2 / transport.volumes.sum()))
    factorized = diags(degree / transport.volumes) @ normalized
    difference = factorized - transport.delta
    weighted_sum = float(degree @ initial)
    report = {
        "time": duration, "initial_mass": mass,
        "conservative_l2_error": norm(conservative - exact),
        "fixed_rate_l2_error": norm(fixed_rate - exact),
        "conservative_relative_mass_drift": abs(float(transport.volumes @ conservative) - mass) / mass,
        "fixed_rate_relative_mass_drift": abs(float(transport.volumes @ fixed_rate) - mass) / mass,
        "fixed_rate_relative_degree_sum_drift": abs(float(degree @ fixed_rate) - weighted_sum) / weighted_sum,
        "fixed_exchange_rate": rate,
        "interior_effective_diffusivity": rate * (length / n)**2 / 6,
        "capacity_factorization_max_error": float(np.max(np.abs(difference.data))) if difference.nnz else 0.,
        "minimum_conservative_concentration": float(conservative.min()),
    }
    slices = {"initial": initial.reshape((n, n, n))[:, :, n // 2],
              "continuum": exact.reshape((n, n, n))[:, :, n // 2],
              "conservative": conservative.reshape((n, n, n))[:, :, n // 2],
              "fixed_rate": fixed_rate.reshape((n, n, n))[:, :, n // 2]}
    return report, slices


def verify_growth(transport, n, mode, length=1., beta=2., da=.02, dh=.4,
                  duration=.25, dt=.01, amplitude=1e-5):
    """Measure a small mode in the nonlinear equations, independently of its prediction."""
    if min(duration, dt, amplitude) <= 0 or not np.isfinite([duration, dt, amplitude]).all():
        raise ValueError("duration, dt, and amplitude must be finite and positive")
    if amplitude >= .01:
        raise ValueError("growth verification requires amplitude < 0.01")
    eigenvalue = fv_eigenvalue(n, mode, length)
    block = gm_jacobian(beta) - eigenvalue * np.diag([da, dh])
    values, vectors = np.linalg.eig(block)
    index = np.argmax(values.real)
    if abs(values[index].imag) > 1e-10:
        raise ValueError("selected mode is oscillatory; a scalar exponential fit is inappropriate")
    chemical = vectors[:, index].real
    spatial = cosine_average(n, mode)
    spatial /= np.sqrt(transport.volumes @ spatial**2 / transport.volumes.sum())
    a, inhibitor = 1 + amplitude * chemical[0] * spatial, 1 + amplitude * chemical[1] * spatial
    count = max(4, int(np.ceil(duration / dt)))
    times = np.linspace(0, duration, count + 1)
    trace = []
    for step in range(count + 1):
        trace.append(float(transport.volumes @ (spatial * (a - 1)) / transport.volumes.sum()))
        if step < count:
            a, inhibitor = integrate_gm(a, inhibitor, transport, duration / count, beta, da, dh)
    measured = float(np.polyfit(times, np.log(np.abs(trace)), 1)[0])
    predicted = float(values[index].real)
    return {"mode": list(mode), "predicted_discrete_growth": predicted,
            "predicted_continuum_growth": float(mode_growth([continuum_eigenvalue(mode, length)], beta, da, dh)[0]),
            "measured_growth": measured, "growth_error": abs(measured - predicted),
            "time": times.tolist(), "projected_activator_amplitude": trace,
            "minimum_activator": float(a.min()), "minimum_inhibitor": float(inhibitor.min())}


def mode_report(transport, n, mode, length, beta, da, dh):
    discrete, continuous = fv_eigenvalue(n, mode, length), continuum_eigenvalue(mode, length)
    spatial = cosine_average(n, mode)
    residual = transport.delta @ spatial + discrete * spatial
    return {"mode": list(mode), "continuum_eigenvalue": continuous,
            "fv_eigenvalue": discrete, "relative_eigenvalue_error": abs(discrete / continuous - 1),
            "operator_eigenpair_relative_residual": float(np.linalg.norm(residual) / np.linalg.norm(spatial) / max(discrete, 1)),
            "continuum_growth": float(mode_growth([continuous], beta, da, dh)[0]),
            "fv_growth": float(mode_growth([discrete], beta, da, dh)[0])}


def unstable_band(beta, da, dh):
    if mode_growth([0], beta, da, dh)[0] >= 0:
        return None
    roots = np.roots([da * dh, beta * da - dh, beta])
    if np.any(np.abs(roots.imag) > 1e-10) or np.any(roots.real <= 0):
        return None
    return sorted(roots.real.tolist())


def supported_modes(n, length, beta, da, dh):
    """All n³ separable FV modes; no dense n³-by-n³ eigendecomposition."""
    axis = 4 * (n / length)**2 * np.sin(np.pi * np.arange(n) / (2 * n))**2
    spectrum = (axis[:, None, None] + axis[None, :, None] + axis[None, None, :]).ravel()
    band = unstable_band(beta, da, dh)
    if band is None:
        return {"unstable_mode_count": 0, "unstable_eigenvalues": [], "unstable_modes": []}
    mask = (spectrum > band[0]) & (spectrum < band[1])
    selected = spectrum[mask]
    indices = np.argwhere(mask.reshape((n, n, n)))
    return {"unstable_mode_count": len(selected), "unstable_eigenvalues": np.sort(selected).tolist(),
            "unstable_modes": [{"mode": index.tolist(), "eigenvalue": float(spectrum[np.ravel_multi_index(index, (n, n, n))])}
                               for index in indices]}


def continuum_supported_modes(length, beta, da, dh):
    band = unstable_band(beta, da, dh)
    if band is None:
        return []
    bound = int(np.ceil(length * np.sqrt(band[1]) / np.pi))
    if bound > 100:
        raise ValueError("too many continuum modes for this benchmark; adjust domain or diffusivities")
    result = []
    for i in range(bound + 1):
        for j in range(bound + 1):
            for k in range(bound + 1):
                value = continuum_eigenvalue((i, j, k), length)
                if band[0] < value < band[1]:
                    result.append({"mode": [i, j, k], "eigenvalue": value})
    return result


def run_benchmark(output, resolutions=(4, 8, 16, 32), length=1., beta=2., da=.02, dh=.4,
                  diffusion_time=1., growth_time=.25, growth_dt=.01, amplitude=1e-5):
    output = Path(output)
    if output.exists():
        raise FileExistsError("output directory already exists; choose a new path")
    if len(resolutions) < 2 or any(type(n) is not int or not 4 <= n <= 64 for n in resolutions):
        raise ValueError("supply at least two resolutions, each an integer between 4 and 64")
    if list(resolutions) != sorted(set(resolutions)):
        raise ValueError("resolutions must be distinct and increasing")
    if min(length, beta, da, dh, diffusion_time, growth_time, growth_dt, amplitude) <= 0 or not np.isfinite(
            [length, beta, da, dh, diffusion_time, growth_time, growth_dt, amplitude]).all():
        raise ValueError("all physical and numerical parameters must be finite and positive")
    if beta <= 1 or amplitude >= .01:
        raise ValueError("require beta > 1 for local stability and amplitude < 0.01 for linear verification")
    exact_modes = continuum_supported_modes(length, beta, da, dh)
    # Check the selected chemical eigensystems before creating output/long work.
    for n in resolutions:
        for mode in ((1, 1, 0), (2, 2, 0)):
            values = np.linalg.eigvals(gm_jacobian(beta) - fv_eigenvalue(n, mode, length) * np.diag([da, dh]))
            if np.max(abs(values.imag)) > 1e-10:
                raise ValueError("selected parameters give oscillatory verification modes")
    output.mkdir(parents=True)
    started = time.perf_counter()
    report = {"schema": 1, "parameters": {"resolutions": list(resolutions), "length": length,
              "beta": beta, "da": da, "dh": dh, "reference_n": resolutions[0],
              "diffusion_time": diffusion_time, "growth_time": growth_time,
              "growth_dt": growth_dt, "amplitude": amplitude},
              "interpretation": "Fixed orthogonal bulk-diffusion control volumes, not dividing biological cells. "
              "Convergence here does not establish a continuum limit of the embryo contact graph, fate, or mechanics.",
              "continuum_unstable_band": unstable_band(beta, da, dh),
              "continuum_unstable_modes": exact_modes, "resolutions": []}
    all_slices = {}
    for n in resolutions:
        transport = cartesian_transport(n, length)
        print(f"n={n}: {n**3:,} compartments; diffusion and early GM growth", flush=True)
        diffusion, slices = diffusion_check(transport, n, length, da, diffusion_time, resolutions[0])
        row = {"n": n, "compartments": n**3, "spacing": length / n,
               "operator_nonzeros": transport.delta.nnz,
               "constant_preservation_residual": float(np.max(np.abs(transport.delta @ np.ones(n**3)))),
               "mass_preservation_residual": float(np.max(np.abs(transport.volumes @ transport.delta))),
               "diffusion": diffusion, **supported_modes(n, length, beta, da, dh),
               "eigenmodes": [mode_report(transport, n, mode, length, beta, da, dh)
                              for mode in ((1, 0, 0), (1, 1, 0), (2, 0, 0), (2, 2, 0))],
               "growth_checks": [verify_growth(transport, n, mode, length, beta, da, dh,
                                                 growth_time, growth_dt, amplitude)
                                 for mode in ((1, 1, 0), (2, 2, 0))]}
        report["resolutions"].append(row)
        all_slices.update({f"n{n}_{key}": value for key, value in slices.items()})
        print(f'  diffusion error={diffusion["conservative_l2_error"]:.3e}; '
              f'fixed-rate error={diffusion["fixed_rate_l2_error"]:.3e}; '
              f'unstable modes={row["unstable_mode_count"]}', flush=True)
        (output / "analysis.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    rows = report["resolutions"]
    errors = np.array([row["diffusion"]["conservative_l2_error"] for row in rows])
    report["observed_diffusion_orders"] = (np.log(errors[:-1] / errors[1:]) /
                                             np.log(np.array(resolutions[1:]) / resolutions[:-1])).tolist()
    report["checks"] = {
        "conservative_mass_drift_below_1e_10": all(row["diffusion"]["conservative_relative_mass_drift"] < 1e-10 for row in rows),
        "diffusion_error_decreases": bool(np.all(np.diff(errors) < 0)),
        "last_diffusion_order_above_1_8": report["observed_diffusion_orders"][-1] > 1.8,
        "measured_growth_error_below_1e_4": all(check["growth_error"] < 1e-4 for row in rows for check in row["growth_checks"]),
        "finest_unstable_count_matches_continuum": rows[-1]["unstable_mode_count"] == len(exact_modes),
        "finest_unstable_mode_indices_match_continuum":
            {tuple(item["mode"]) for item in rows[-1]["unstable_modes"]} ==
            {tuple(item["mode"]) for item in exact_modes},
    }
    report["elapsed_seconds"] = time.perf_counter() - started
    (output / "analysis.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    np.savez_compressed(output / "diffusion_slices.npz", **all_slices)
    plot_report(report, all_slices, output)
    return report


def plot_report(report, slices, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = report["resolutions"]
    ns = np.array([row["n"] for row in rows])
    green, coral = "#237969", "#bc7158"
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for key, label, color in [("conservative_l2_error", "Conservative finite volume", green),
                               ("fixed_rate_l2_error", "Fixed normalized exchange", coral)]:
        axes[0, 0].loglog(ns, [row["diffusion"][key] for row in rows], "o-", color=color, label=label)
    error = rows[0]["diffusion"]["conservative_l2_error"]
    axes[0, 0].loglog(ns, error * (ns[0] / ns)**2, "k:", alpha=.5, label="Second-order reference")
    axes[0, 0].set(title="Approach to the same diffusion equation", xlabel="Compartments per axis, n", ylabel=f"Volume-weighted L2 error at t = {report['parameters']['diffusion_time']:g}")
    axes[0, 0].legend(fontsize=8)
    axes[0, 1].plot(ns, [row["diffusion"]["interior_effective_diffusivity"] for row in rows], "o-", color=coral, label="Fixed normalized rate")
    axes[0, 1].axhline(report["parameters"]["da"], color=green, linestyle="--", label="Prescribed bulk diffusivity")
    axes[0, 1].set(title="Fixed per-cell rates change the physical scale", xlabel="Compartments per axis, n", ylabel="Interior effective diffusivity")
    axes[0, 1].legend(fontsize=8)
    for index, color, label in [(1, green, "Spatial mode (1,1,0)"), (2, coral, "Spatial mode (2,0,0)")]:
        axes[1, 0].plot(ns, [row["eigenmodes"][index]["fv_growth"] for row in rows], "o-", color=color, label=label)
        axes[1, 0].axhline(rows[0]["eigenmodes"][index]["continuum_growth"], color=color, linestyle=":")
    axes[1, 0].axhline(0, color="gray", linewidth=.7)
    axes[1, 0].set(title="Finite spatial modes change growth predictions", xlabel="Compartments per axis, n", ylabel="Early GM growth rate (dotted: continuum)")
    axes[1, 0].legend(fontsize=8)
    for index, color, label in [(0, green, "(1,1,0)"), (1, coral, "(2,2,0)")]:
        axes[1, 1].semilogy(ns, [max(row["growth_checks"][index]["growth_error"], 1e-16) for row in rows], "o-", color=color, label=label)
    axes[1, 1].set(title="Nonlinear solver versus discrete linear prediction", xlabel="Compartments per axis, n", ylabel="Absolute growth-rate error")
    axes[1, 1].legend(fontsize=8)
    fig.suptitle("Discrete-to-continuum bridge · fixed 3D cube, no-flux boundaries")
    fig.savefig(output / "convergence.png", dpi=160)
    plt.close(fig)
    n = rows[-1]["n"]
    keys = ["initial", "continuum", "conservative", "fixed_rate"]
    titles = ["Initial field", "Continuum cell averages", "Conservative exchange", "Fixed normalized rate"]
    low = min(slices[f"n{n}_{key}"].min() for key in keys)
    high = max(slices[f"n{n}_{key}"].max() for key in keys)
    fig, axes = plt.subplots(1, 4, figsize=(12, 3.8), constrained_layout=True)
    for ax, key, title in zip(axes, keys, titles):
        picture = ax.imshow(slices[f"n{n}_{key}"].T, origin="lower", extent=(0, report["parameters"]["length"], 0, report["parameters"]["length"]),
                            vmin=low, vmax=high, cmap="viridis")
        ax.set(title=title, xlabel="x", ylabel="y")
    fig.colorbar(picture, ax=axes, shrink=.7, label="Concentration (shared scale)")
    fig.suptitle(f"Same prescribed smooth field · n={n} · slice z={(n // 2 + .5) * report['parameters']['length'] / n:.3f}")
    fig.savefig(output / "diffusion.png", dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--resolutions", nargs="+", type=int, default=[4, 8, 16, 32])
    parser.add_argument("--length", type=float, default=1.)
    parser.add_argument("--beta", type=float, default=2.)
    parser.add_argument("--da", type=float, default=.02)
    parser.add_argument("--dh", type=float, default=.4)
    parser.add_argument("--diffusion-time", type=float, default=1.)
    parser.add_argument("--growth-time", type=float, default=.25)
    parser.add_argument("--growth-dt", type=float, default=.01)
    parser.add_argument("--amplitude", type=float, default=1e-5)
    args = parser.parse_args()
    try:
        result = run_benchmark(**vars(args))
    except (ValueError, FileExistsError) as error:
        parser.error(str(error))
    print(json.dumps(result["checks"], indent=2), flush=True)
    if not all(result["checks"].values()):
        print("Some benchmark criteria were not met. Inspect analysis.json before drawing conclusions.", flush=True)


if __name__ == "__main__":
    main()
