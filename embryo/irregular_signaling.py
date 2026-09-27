"""Full discrete L-prism spectra and early Gierer-Meinhardt growth checks.

The orthogonal mesh is a product of an L-shaped cross-section and a z interval.
Its weighted Laplacian is a Kronecker sum. Only the 2D and 1D factors are dense;
no dense 3D operator or full 3D eigenvector matrix is allocated.
"""

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import time

import numpy as np
from scipy import linalg, sparse

from .continuum import unstable_band
from .irregular import graded_edges, l_prism_mask, _positive
from .signaling import gm_jacobian, mode_growth
from .transport import Transport, cartesian_transport, masked_cartesian_transport, integrate_gm


@dataclass
class PrismSpectrum:
    transport: Transport
    edges: tuple
    mask: np.ndarray
    eigenvalues: np.ndarray
    pairs: np.ndarray
    xy_values: np.ndarray
    z_values: np.ndarray
    xy_vectors: np.ndarray
    z_vectors: np.ndarray
    factorization_relative_residual: float

    def mode(self, index):
        """A physical concentration mode with unit volume-weighted RMS."""
        if type(index) is not int or not 0 <= index < len(self.eigenvalues):
            raise ValueError("mode index is outside the discrete spectrum")
        i, j = self.pairs[index]
        q = (self.xy_vectors[:, i, None] * self.z_vectors[None, :, j]).ravel()
        u = q / np.sqrt(self.transport.volumes / self.transport.volumes.sum())
        if u[np.argmax(np.abs(u))] < 0:
            u = -u
        return u


def prism_spectrum(n, length=1., grading=.35):
    """Complete discrete spectrum, including multiplicities, on the L-prism.

    This product construction is specific to the present extruded geometry.
    Limit n to 32 here because the complete 2D factor is diagonalized densely.
    Orthogonality uses compartment capacity, not normalized graph degree.
    """
    if type(n) is not int or n < 4 or n > 32 or n % 2:
        raise ValueError("spectral resolutions must be even integers between 4 and 32")
    edges = graded_edges(n, length, grading)
    mask = l_prism_mask(n)
    transport = masked_cartesian_transport(edges=edges, mask=mask)
    xy = masked_cartesian_transport(edges=(edges[0], edges[1], np.array([0., 1.])),
                                    mask=mask[:, :, :1])
    z = cartesian_transport(edges=(np.array([0., 1.]), np.array([0., 1.]), edges[2]))
    factored = (sparse.kron(xy.symmetric, sparse.eye(n), format="csr")
                + sparse.kron(sparse.eye(len(xy.volumes)), z.symmetric, format="csr"))
    difference = transport.symmetric - factored
    scale = float(np.max(np.abs(transport.symmetric.data)))
    residual = float(np.max(np.abs(difference.data))) / scale if difference.nnz else 0.
    if residual > 1e-12:
        raise ValueError("3D operator does not match the proposed product factorization")
    xy_values, xy_vectors = linalg.eigh(xy.symmetric.toarray())
    z_values, z_vectors = linalg.eigh(z.symmetric.toarray())
    for values in (xy_values, z_values):
        tolerance = 1e-10 * max(1., float(values[-1]))
        if abs(values[0]) > tolerance or values[1] <= 0:
            raise ValueError("expected one constant mode in each connected factor")
        values[0] = 0.  # Remove only the known constant-mode roundoff.
    sums = (xy_values[:, None] + z_values[None, :]).ravel()
    order = np.argsort(sums, kind="stable")
    pairs = np.column_stack(np.unravel_index(order, (len(xy_values), len(z_values))))
    return PrismSpectrum(transport, edges, mask, sums[order], pairs, xy_values, z_values,
                         xy_vectors, z_vectors, residual)


def stability_report(spectrum, beta=2., da=.02, dh=.4):
    beta, da, dh = (_positive(name, value) for name, value in (("beta", beta), ("da", da), ("dh", dh)))
    rates = mode_growth(spectrum.eigenvalues, beta, da, dh)
    unstable = np.flatnonzero((spectrum.eigenvalues > 0) & (rates > 1e-10))
    return {"local_growth_rate": float(rates[0]), "local_stable": bool(rates[0] < -1e-10),
            "stationary_instability_band": unstable_band(beta, da, dh),
            "eigenvalues": spectrum.eigenvalues.tolist(), "growth_rates": rates.tolist(),
            "unstable_indices": unstable.tolist(), "unstable_mode_count": len(unstable),
            "diffusion_driven_instability": bool(rates[0] < -1e-10 and len(unstable)),
            "spectral_gap": float(spectrum.eigenvalues[1]),
            "maximum_spatial_growth": float(np.max(rates[1:]))}


def select_probes(spectrum, report, beta=2., da=.02, dh=.4):
    """Fastest growing and slowest decaying nonconstant real chemical modes."""
    rates = np.asarray(report["growth_rates"])
    probes = []
    if report["unstable_indices"]:
        index = int(max(report["unstable_indices"], key=lambda i: rates[i]))
        probes.append(("growing", index))
    candidates = np.flatnonzero((spectrum.eigenvalues > 0) & (rates < -1e-8))
    for index in candidates[np.argsort(rates[candidates])[::-1]]:
        values = np.linalg.eigvals(gm_jacobian(beta) - spectrum.eigenvalues[index] * np.diag([da, dh]))
        if np.max(np.abs(values.imag)) < 1e-10:
            probes.append(("decaying", int(index)))
            break
    return probes


def verify_growth(spectrum, index, beta=2., da=.02, dh=.4, duration=.25, dt=.01, amplitude=1e-5):
    """Compare nonlinear SSP-RK2 dynamics with an independently diagonalized mode."""
    duration, dt, amplitude = (_positive(name, value) for name, value in
                               (("duration", duration), ("dt", dt), ("amplitude", amplitude)))
    beta, da, dh = (_positive(name, value) for name, value in (("beta", beta), ("da", da), ("dh", dh)))
    if amplitude >= .01:
        raise ValueError("early-growth probes require amplitude < 0.01")
    spatial = spectrum.mode(index)
    value = float(spectrum.eigenvalues[index])
    chemical_values, chemical_vectors = np.linalg.eig(gm_jacobian(beta) - value * np.diag([da, dh]))
    selected = int(np.argmax(chemical_values.real))
    if abs(chemical_values[selected].imag) > 1e-10:
        raise ValueError("a scalar growth fit requires a real chemical eigenvalue")
    rate = float(chemical_values[selected].real)
    chemical = chemical_vectors[:, selected].real
    chemical /= np.linalg.norm(chemical)
    if chemical[0] < 0:
        chemical = -chemical
    if abs(chemical[0]) < 1e-10:
        raise ValueError("activator projection is too small for a scalar fit")
    a, h = 1 + amplitude * chemical[0] * spatial, 1 + amplitude * chemical[1] * spatial
    if min(a.min(), h.min()) <= 0:
        raise ValueError("perturbation must retain positive concentrations")
    volume = spectrum.transport.volumes
    weights = volume / volume.sum()
    count = max(4, int(np.ceil(duration / dt)))
    times = np.linspace(0., duration, count + 1)
    projections, inhibitor_projections = [], []
    for step in range(count + 1):
        projections.append(float(weights @ (spatial * (a - 1))))
        inhibitor_projections.append(float(weights @ (spatial * (h - 1))))
        if step < count:
            a, h = integrate_gm(a, h, spectrum.transport, duration / count, beta, da, dh)
    projections = np.array(projections)
    if np.any(projections <= 0):
        raise FloatingPointError("probe no longer supports a positive exponential amplitude fit")
    measured = float(np.polyfit(times, np.log(projections), 1)[0])
    linear_a = amplitude * chemical[0] * np.exp(rate * duration) * spatial
    linear_h = amplitude * chemical[1] * np.exp(rate * duration) * spatial
    relative_linear_error = float(np.sqrt(weights @ ((a - 1 - linear_a)**2 + (h - 1 - linear_h)**2)) / amplitude)
    residual = spectrum.transport.delta @ spatial + value * spatial
    residual_rms = float(np.sqrt(weights @ residual**2))
    eigen_residual = residual_rms / max(value, 1. / spectrum.edges[0][-1]**2)
    report = {"mode_index": index, "factor_indices": spectrum.pairs[index].tolist(),
              "spatial_eigenvalue": value, "predicted_growth": rate, "measured_growth": measured,
              "growth_error": abs(measured - rate), "relative_linear_field_error": relative_linear_error,
              "relative_eigen_residual": eigen_residual,
              "volume_weighted_mode_mean": float(weights @ spatial),
              "volume_weighted_mode_rms": float(np.sqrt(weights @ spatial**2)),
              "chemical_eigenvector": chemical.tolist(), "time": times.tolist(),
              "activator_projection": projections.tolist(), "inhibitor_projection": inhibitor_projections,
              "linear_activator_projection": (amplitude * chemical[0] * np.exp(rate * times)).tolist(),
              "minimum_activator": float(a.min()), "minimum_inhibitor": float(h.min())}
    return report, {"spatial_mode": spatial, "activator_final": a, "inhibitor_final": h}


def run_benchmark(output, resolutions=(4, 8, 16, 32), length=1., grading=.35,
                  beta=2., da=.02, dh=.4, duration=.25, dt=.01, amplitude=1e-5):
    output = Path(output)
    if output.exists():
        raise FileExistsError("output directory already exists; choose a new path")
    if not isinstance(resolutions, (tuple, list)) or len(resolutions) < 2:
        raise ValueError("supply at least two increasing resolutions")
    if any(type(n) is not int or n < 4 or n > 32 or n % 2 for n in resolutions):
        raise ValueError("resolutions must be even integers between 4 and 32")
    if list(resolutions) != sorted(set(resolutions)):
        raise ValueError("resolutions must be distinct and increasing")
    length, beta, da, dh, duration, dt, amplitude = (_positive(name, value) for name, value in
        (("length", length), ("beta", beta), ("da", da), ("dh", dh), ("duration", duration), ("dt", dt), ("amplitude", amplitude)))
    graded_edges(resolutions[0], length, grading)
    if beta <= 1 or amplitude >= .01:
        raise ValueError("require beta > 1 and amplitude < 0.01 for diffusion-driven linear verification")
    output.mkdir(parents=True)
    started = time.perf_counter()
    report = {"schema": 1, "parameters": {"resolutions": list(resolutions), "length": length,
              "grading": float(grading), "beta": beta, "da": da, "dh": dh,
              "duration": duration, "dt": dt, "amplitude": amplitude},
              "scope": "Complete discrete spectra on a fixed orthogonal L-prism, followed by tiny-mode "
                       "nonlinear probes. Refinement comparisons are not an exact continuum spectrum, "
                       "persistent nonlinear pattern, identity, or shape validation.",
              "acceptance": {"factorization_residual": 1e-12, "relative_eigen_residual": 1e-9,
                             "growth_error": 1e-4, "relative_linear_field_error": 1e-3,
                             "last_low_spectrum_relative_change": .03, "low_modes_compared": 12},
              "rows": []}
    archive = {}
    for n in resolutions:
        print(f"n={n}: computing complete weighted spectrum", flush=True)
        spectrum = prism_spectrum(n, length, grading)
        linear = stability_report(spectrum, beta, da, dh)
        row = {"n": n, "compartments": len(spectrum.eigenvalues),
               "factorization_relative_residual": spectrum.factorization_relative_residual,
               "xy_eigenvalues": spectrum.xy_values.tolist(), "z_eigenvalues": spectrum.z_values.tolist(),
               "factor_indices": spectrum.pairs.tolist(), **linear, "probes": []}
        for label, index in select_probes(spectrum, linear, beta, da, dh):
            probe, fields = verify_growth(spectrum, index, beta, da, dh, duration, dt, amplitude)
            row["probes"].append({"label": label, **probe})
            for key, values in fields.items():
                archive[f"n{n}_{label}_{key}"] = values
        archive[f"n{n}_mask"] = spectrum.mask
        archive[f"n{n}_volumes"] = spectrum.transport.volumes
        for axis, edges in zip("xyz", spectrum.edges):
            archive[f"n{n}_{axis}_edges"] = edges
        report["rows"].append(row)
        print(f'  gap={row["spectral_gap"]:.7g}, unstable={row["unstable_mode_count"]}, '
              f'maximum growth={row["maximum_spatial_growth"]:.7g}', flush=True)
    rows = report["rows"]
    comparisons = []
    for coarse, fine in zip(rows, rows[1:]):
        a, b = np.array(coarse["eigenvalues"][1:13]), np.array(fine["eigenvalues"][1:13])
        comparisons.append({"coarse_n": coarse["n"], "fine_n": fine["n"],
                            "max_relative_low_eigenvalue_change": float(np.max(np.abs(a - b) / b)),
                            "unstable_count_agrees": coarse["unstable_mode_count"] == fine["unstable_mode_count"]})
    report["refinement"] = comparisons
    probes = [probe for row in rows for probe in row["probes"]]
    report["checks"] = {
        "product_matches_3d_operator": all(row["factorization_relative_residual"] < 1e-12 for row in rows),
        "local_equilibrium_stable": all(row["local_stable"] for row in rows),
        "growing_and_decaying_probes_present": all({probe["label"] for probe in row["probes"]} == {"growing", "decaying"} for row in rows),
        "physical_modes_satisfy_3d_eigenproblem": bool(probes) and all(probe["relative_eigen_residual"] < 1e-9 for probe in probes),
        "measured_growth_matches_prediction": bool(probes) and all(probe["growth_error"] < 1e-4 for probe in probes),
        "nonlinear_fields_remain_near_linear_prediction": bool(probes) and all(probe["relative_linear_field_error"] < 1e-3 for probe in probes),
        "probe_concentrations_positive": bool(probes) and all(min(probe["minimum_activator"], probe["minimum_inhibitor"]) > 0 for probe in probes),
        "finest_low_spectrum_change_below_3_percent": comparisons[-1]["max_relative_low_eigenvalue_change"] < .03,
        "finest_unstable_counts_agree": comparisons[-1]["unstable_count_agrees"],
    }
    report["elapsed_seconds"] = time.perf_counter() - started
    (output / "analysis.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    np.savez_compressed(output / "modes.npz", **archive)
    plot_report(report, archive, output)
    return report


def plot_report(report, archive, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = report["rows"]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    colors = plt.cm.viridis(np.linspace(.1, .85, len(rows)))
    for row, color in zip(rows, colors):
        values, rates = np.array(row["eigenvalues"]), np.array(row["growth_rates"])
        band = row["stationary_instability_band"]
        bound = band[1] * 1.35 if band is not None else values[min(15, len(values) - 1)]
        chosen = values <= bound
        axes[0].plot(values[chosen], rates[chosen], "o", color=color, label=f'n={row["n"]}', markersize=4)
    axes[0].axhline(0, color="gray", linewidth=.7)
    axes[0].set(title="Actual discrete spatial modes", xlabel="Weighted Laplacian eigenvalue", ylabel="Predicted GM growth rate")
    axes[0].legend(fontsize=8)
    ns = [row["n"] for row in rows]
    for i in range(1, 7):
        axes[1].plot(ns, [row["eigenvalues"][i] for row in rows], "o-", label=f"rank {i}")
    axes[1].set(title="Low spatial spectrum under refinement", xlabel="Bounding-box resolution n", ylabel="Eigenvalue")
    axes[1].legend(fontsize=8, ncol=2)
    for probe in rows[-1]["probes"]:
        observed = np.array(probe["activator_projection"])
        exact = np.array(probe["linear_activator_projection"])
        axes[2].plot(probe["time"], observed / observed[0], "o", markersize=3, label=f'{probe["label"]}: nonlinear')
        axes[2].plot(probe["time"], exact / exact[0], "--", label=f'{probe["label"]}: linear')
    axes[2].set(title=f'Early growth, finest n={ns[-1]}', xlabel="Time", ylabel="Activator projection / initial")
    axes[2].legend(fontsize=8)
    fig.suptitle("Gierer-Meinhardt instability on a fixed graded L-prism")
    fig.savefig(output / "spectrum.png", dpi=160)
    plt.close(fig)
    n = ns[-1]
    probes = rows[-1]["probes"]
    if not probes:
        return
    fig, axes = plt.subplots(1, len(probes), figsize=(5 * len(probes), 4.5), squeeze=False, constrained_layout=True)
    mask = archive[f"n{n}_mask"]
    x, y, z = (archive[f"n{n}_{axis}_edges"] for axis in "xyz")
    for ax, probe in zip(axes[0], probes):
        field = archive[f'n{n}_{probe["label"]}_spatial_mode']
        bound = float(np.max(np.abs(field)))
        cube = np.full(mask.shape, np.nan)
        cube[mask] = field
        slab = int(np.argmax(np.nanmean(cube**2, axis=(0, 1))))
        picture = ax.pcolormesh(x, y, np.ma.masked_invalid(cube[:, :, slab].T), cmap="RdBu_r", vmin=-bound, vmax=bound, shading="flat")
        ax.set(aspect="equal", xlabel="x", ylabel="y", title=f'{probe["label"].capitalize()} mode: lambda={probe["spatial_eigenvalue"]:.3f}\nr={probe["predicted_growth"]:+.4f}, z=[{z[slab]:.3f}, {z[slab+1]:.3f}]')
        fig.colorbar(picture, ax=ax, shrink=.8, label="Mode, unit volume RMS")
    fig.suptitle(f'Physical concentration modes, n={n}; each panel shows its highest-energy z slab\nSign and basis within degenerate eigenspaces are arbitrary')
    fig.savefig(output / "modes.png", dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--resolutions", nargs="+", type=int, default=[4, 8, 16, 32])
    for name, default in (("length", 1.), ("grading", .35), ("beta", 2.), ("da", .02),
                          ("dh", .4), ("duration", .25), ("dt", .01), ("amplitude", 1e-5)):
        parser.add_argument("--" + name, type=float, default=default)
    args = parser.parse_args()
    try:
        report = run_benchmark(**vars(args))
    except (ValueError, FileExistsError) as error:
        parser.error(str(error))
    print(json.dumps(report["checks"], indent=2), flush=True)
    if not all(report["checks"].values()):
        raise SystemExit("Some spectral/growth criteria were not met; inspect analysis.json.")


if __name__ == "__main__":
    main()
