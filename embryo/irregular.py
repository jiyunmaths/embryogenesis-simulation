"""Conservative diffusion on a graded 3D L-shaped domain.

This is an irregular-boundary transport benchmark with exact smooth Neumann
solutions, not a model of emergent shape or arbitrary nonorthogonal cells.
"""

import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import expm_multiply

from .transport import masked_cartesian_transport


MODES = ((2, 0, 0), (0, 2, 1), (2, 2, 2))
AMPLITUDES = (.10, .07, .04)


def _resolution(n):
    if type(n) is not int or n < 4 or n > 64 or n % 2:
        raise ValueError("resolutions must be even integers between 4 and 64")


def _positive(name, value):
    if (isinstance(value, (bool, np.bool_)) or not np.isscalar(value)
            or not np.isrealobj(value)):
        raise ValueError(f"{name} must be finite and positive")
    try:
        value = float(value)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError(f"{name} must be finite and positive") from error
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def graded_edges(n, length=1., grading=.35):
    """Smooth orthogonal meshes; integer-factor refinements are nested.

    x/L = s + g sin(2 pi s)/(2 pi), with axis coefficients
    (grading, -.6*grading, .4*grading). The map is monotone for |g|<1.
    This benchmark limits grading to [0,.8) to avoid nearly collapsed cells.
    """
    _resolution(n)
    length = _positive("length", length)
    if (isinstance(grading, (bool, np.bool_)) or not np.isscalar(grading)
            or not np.isrealobj(grading)):
        raise ValueError("grading must be finite and in [0, 0.8)")
    try:
        grading = float(grading)
    except (ValueError, TypeError, OverflowError) as error:
        raise ValueError("grading must be finite and in [0, 0.8)") from error
    if not np.isfinite(grading) or not 0 <= grading < .8:
        raise ValueError("grading must be finite and in [0, 0.8)")
    s = np.linspace(0., 1., n + 1)
    edges = []
    for coefficient in (grading, -.6 * grading, .4 * grading):
        axis = length * (s + coefficient * np.sin(2 * np.pi * s) / (2 * np.pi))
        axis[[0, n // 2, n]] = (0., length / 2, length)
        edges.append(axis)
    return tuple(edges)


def l_prism_mask(n):
    """Retain [0,L]^3 except x>L/2 AND y>L/2, for every z layer."""
    _resolution(n)
    mask = np.ones((n, n, n), dtype=bool)
    mask[n // 2:, n // 2:, :] = False
    return mask


def cosine_average(edges, mode, mask):
    """Exact averages of a smooth L-prism Neumann mode on rectangular cells.

    Even x/y indices ensure zero flux at the reentrant walls. These modes are
    only a subset of the true domain spectrum, and their averages are not
    discrete eigenvectors on a graded mesh. No quadrature is used here.
    """
    indices = np.asarray(mode)
    if (indices.shape != (3,) or not np.issubdtype(indices.dtype, np.integer)
            or np.any(indices < 0) or np.any(indices[:2] % 2)):
        raise ValueError("mode needs nonnegative integer indices, with even x/y indices")
    if not isinstance(edges, (tuple, list)) or len(edges) != 3:
        raise ValueError("edges must contain three axis vectors")
    if any(np.iscomplexobj(axis) for axis in edges):
        raise ValueError("edges must be real")
    axes = tuple(np.asarray(axis, dtype=float) for axis in edges)
    if any(axis.ndim != 1 or len(axis) < 2 or not np.isfinite(axis).all()
           or axis[0] != 0 or np.any(np.diff(axis) <= 0) for axis in axes):
        raise ValueError("edges must increase from zero and be finite")
    length = axes[0][-1]
    if any(axis[-1] != length for axis in axes):
        raise ValueError("axis edges must span the same length")
    shape = tuple(len(axis) - 1 for axis in axes)
    if (not isinstance(mask, np.ndarray) or mask.dtype != bool
            or mask.shape != shape or not np.any(mask)):
        raise ValueError("mask must be a nonempty boolean array matching edges")
    if np.any(indices >= np.array(shape)):
        raise ValueError("mode must be resolved on every axis")
    averages = []
    for axis, index in zip(axes, indices):
        widths = np.diff(axis)
        centers = axis[:-1] + widths / 2
        averages.append(np.cos(index * np.pi * centers / length)
                        * np.sinc(index * widths / (2 * length)))
    field = (averages[0][:, None, None] * averages[1][None, :, None]
             * averages[2][None, None, :])
    return field[mask]


def _rms(values, volumes):
    return float(np.sqrt(volumes @ np.square(values) / volumes.sum()))


def diffusion_check(transport, edges, mask, diffusivity=.02, duration=1.):
    """Compare semidiscrete diffusion with exact continuum cell averages.

    The matrix exponential isolates spatial error from explicit time-step error.
    The same continuous initial field is averaged afresh at every resolution.
    """
    diffusivity, duration = _positive("diffusivity", diffusivity), _positive("duration", duration)
    n = mask.shape[0] if isinstance(mask, np.ndarray) and mask.ndim == 3 else 0
    _resolution(n)
    if not np.array_equal(mask, l_prism_mask(n)):
        raise ValueError("diffusion reference requires the L-prism mask")
    spatial = [cosine_average(edges, mode, mask) for mode in MODES]
    length = float(edges[0][-1])
    if any(len(axis) != n + 1 or axis[n // 2] != length / 2 for axis in edges):
        raise ValueError("the midpoint boundary must be resolved on each axis")
    if transport.volumes.shape != spatial[0].shape:
        raise ValueError("transport does not match the active compartments")
    eigenvalues = np.array([np.pi**2 * np.dot(mode, mode) / length**2 for mode in MODES])
    initial = 1 + sum(amplitude * vector for amplitude, vector in zip(AMPLITUDES, spatial))
    exact = 1 + sum(amplitude * vector * np.exp(-diffusivity * duration * eigenvalue)
                    for amplitude, vector, eigenvalue in zip(AMPLITUDES, spatial, eigenvalues))
    numerical = expm_multiply(diffusivity * duration * transport.delta, initial)
    error = numerical - exact
    mass, final_mass = float(transport.volumes @ initial), float(transport.volumes @ numerical)
    report = {
        "l2_error": _rms(error, transport.volumes),
        "linf_error": float(np.max(np.abs(error))),
        "relative_mass_drift": abs(final_mass - mass) / abs(mass),
        "initial_mass": mass, "final_mass": final_mass,
        "exact_mass": float(transport.volumes @ exact),
        "minimum_concentration": float(numerical.min()),
        "mode_consistency_rms": [
            _rms(transport.delta @ vector + eigenvalue * vector, transport.volumes)
            for vector, eigenvalue in zip(spatial, eigenvalues)],
    }
    return report, {"initial": initial, "exact": exact, "numerical": numerical, "error": error}


def arm_pulse_check(n=8, length=1., grading=.35, diffusivity=.02):
    """Test exchange around the notch, which the smooth even modes do not probe.

    A unit concentration in the upper-left arm must travel through the lower-
    left junction into the lower-right arm. This is a topology/equilibration
    diagnostic on one modest mesh, not an independent convergence reference.
    Times are scaled by L²/D, so the test is unchanged by units.
    """
    edges = graded_edges(n, length, grading)
    diffusivity = _positive("diffusivity", diffusivity)
    mask = l_prism_mask(n)
    transport = masked_cartesian_transport(edges=edges, mask=mask)
    source = transport.centers[:, 1] > length / 2
    destination = transport.centers[:, 0] > length / 2
    initial = source.astype(float)
    volumes = transport.volumes
    mass = float(volumes @ initial)
    mean = mass / float(volumes.sum())
    initial_rms = _rms(initial - mean, volumes)
    dimensionless_times = (0., .05, .5, 5.)
    traces, fields = [], {}
    for diffusion_time in dimensionless_times:
        field = (initial.copy() if diffusion_time == 0 else
                 expm_multiply(diffusion_time * length**2 * transport.delta, initial))
        fields[f"pulse_{diffusion_time:g}"] = field
        traces.append({"diffusion_time": diffusion_time,
                       "physical_time": diffusion_time * length**2 / diffusivity,
                       "destination_mass": float(volumes[destination] @ field[destination]),
                       "relative_mass_drift": abs(float(volumes @ field) - mass) / mass,
                       "relative_equilibrium_rms": _rms(field - mean, volumes) / initial_rms,
                       "minimum_concentration": float(field.min())})
    report = {"n": n, "compartments": len(volumes),
              "connected_components": int(connected_components(transport.conductance, directed=False, return_labels=False)),
              "initial_mass": mass, "equilibrium_concentration": mean,
              "destination_equilibrium_mass": float(volumes[destination].sum() * mean),
              "trace": traces}
    return report, {"edges": edges, "mask": mask, **fields}


def observed_orders(errors, resolutions, error_floor=1e-13):
    """Report unresolved convergence rates as None, rather than log(0) or noise."""
    orders = []
    for coarse, fine, n0, n1 in zip(errors, errors[1:], resolutions, resolutions[1:]):
        if min(coarse, fine) <= error_floor or not np.isfinite([coarse, fine]).all():
            orders.append(None)
        else:
            orders.append(float((np.log(coarse) - np.log(fine)) / np.log(n1 / n0)))
    return orders


def run_benchmark(output, resolutions=(4, 8, 16, 32), length=1., grading=.35,
                  diffusivity=.02, duration=1.):
    output = Path(output)
    if output.exists():
        raise FileExistsError("output directory already exists; choose a new path")
    if not isinstance(resolutions, (tuple, list)) or len(resolutions) < 2:
        raise ValueError("supply at least two increasing resolutions")
    for n in resolutions:
        _resolution(n)
    if list(resolutions) != sorted(set(resolutions)):
        raise ValueError("resolutions must be distinct and increasing")
    length = _positive("length", length)
    diffusivity, duration = _positive("diffusivity", diffusivity), _positive("duration", duration)
    graded_edges(resolutions[0], length, grading)  # Validate before creating output.
    grading = float(grading)
    output.mkdir(parents=True)
    started = time.perf_counter()
    report = {
        "schema": 1,
        "parameters": {"resolutions": list(resolutions), "length": length,
                       "grading": grading, "diffusivity": diffusivity, "duration": duration,
                       "modes": [list(mode) for mode in MODES], "amplitudes": list(AMPLITUDES)},
        "domain": "[0,L]^3 excluding x>L/2 AND y>L/2, with homogeneous Neumann boundaries",
        "interpretation": "Smooth pure-diffusion verification on a fixed nonconvex domain with "
                          "nonuniform orthogonal rectangular compartments. Prescribed geometry, not emergent "
                          "embryo shape; no moving cells, nonlinear GM patterns, or complete spectrum claim.",
        "acceptance": {"relative_mass_drift": 1e-10, "scaled_operator_residual": 1e-12,
                       "relative_volume_error": 1e-12, "last_observed_order_min": 1.7,
                       "pulse_relative_equilibrium_rms": 1e-6, "order_error_floor": 1e-13},
        "families": {},
    }
    archive = {}
    for name, coefficient in (("uniform_l_prism", 0.), ("graded_l_prism", grading)):
        rows = []
        for n in resolutions:
            edges, mask = graded_edges(n, length, coefficient), l_prism_mask(n)
            transport = masked_cartesian_transport(edges=edges, mask=mask)
            print(f"{name}: n={n}, {len(transport.volumes):,} compartments", flush=True)
            diffusion, fields = diffusion_check(transport, edges, mask, diffusivity, duration)
            volume = transport.volumes
            widths = np.concatenate([np.diff(axis) for axis in edges])
            constant_residual = float(np.max(np.abs(transport.delta @ np.ones(len(volume)))))
            mass_residual = float(np.max(np.abs(volume @ transport.delta)))
            max_exit = float(np.max(-transport.delta.diagonal()))
            degree = np.asarray(transport.conductance.sum(axis=1)).ravel()
            row = {"n": n, "compartments": len(volume), "domain_volume": float(volume.sum()),
                   "relative_volume_error": abs(float(volume.sum()) / (.75 * length**3) - 1),
                   "min_volume": float(volume.min()), "max_volume": float(volume.max()),
                   "volume_ratio": float(volume.max() / volume.min()),
                   "min_spacing": float(widths.min()), "max_spacing": float(widths.max()),
                   "operator_nonzeros": transport.delta.nnz,
                   "connected_components": int(connected_components(transport.conductance, directed=False, return_labels=False)),
                   "constant_residual": constant_residual,
                   "mass_residual": mass_residual,
                   "scaled_constant_residual": constant_residual / max_exit,
                   "scaled_mass_residual": mass_residual / float(degree.max()),
                   **diffusion}
            rows.append(row)
            prefix = f"{name}_n{n}"
            archive.update({f"{prefix}_{key}": value for key, value in fields.items()})
            archive[f"{prefix}_mask"] = mask
            for axis, values in zip("xyz", edges):
                archive[f"{prefix}_{axis}_edges"] = values
            print(f'  volume-RMS error={row["l2_error"]:.6g}; mass drift={row["relative_mass_drift"]:.3g}', flush=True)
        errors = np.array([row["l2_error"] for row in rows])
        # The grading law is fixed. n is a refinement parameter; report h_max
        # as well, since coarse-grid n ratios and h_max ratios differ slightly.
        orders = observed_orders(errors, resolutions)
        report["families"][name] = {"grading": coefficient, "rows": rows,
                                    "observed_orders": orders}
    pulse, pulse_fields = arm_pulse_check(length=length, grading=grading, diffusivity=diffusivity)
    report["arm_pulse"] = pulse
    for key, value in pulse_fields.items():
        if key == "edges":
            for axis, axis_edges in zip("xyz", value):
                archive[f"pulse_{axis}_edges"] = axis_edges
        else:
            archive[f"pulse_{key}" if key == "mask" else key] = value
    all_rows = [row for family in report["families"].values() for row in family["rows"]]
    report["checks"] = {
        "domain_volume_correct": all(row["relative_volume_error"] < 1e-12 for row in all_rows),
        "domain_connected": all(row["connected_components"] == 1 for row in all_rows) and pulse["connected_components"] == 1,
        "constants_and_amount_preserved": all(max(row["scaled_constant_residual"], row["scaled_mass_residual"]) < 1e-12 for row in all_rows),
        "diffusion_mass_drift_below_1e_10": all(row["relative_mass_drift"] < 1e-10 for row in all_rows),
        "diffusion_positive": all(row["minimum_concentration"] > 0 for row in all_rows),
        "errors_decrease_in_both_families": all(all(a["l2_error"] > b["l2_error"] for a, b in zip(family["rows"], family["rows"][1:])) for family in report["families"].values()),
        "last_order_above_1_7_in_both_families": all(family["observed_orders"][-1] is not None and family["observed_orders"][-1] > 1.7 for family in report["families"].values()),
        "pulse_reaches_other_arm": pulse["trace"][1]["destination_mass"] > 0,
        "pulse_equilibrates": pulse["trace"][-1]["relative_equilibrium_rms"] < 1e-6,
        "pulse_mass_preserved_and_nonnegative": all(item["relative_mass_drift"] < 1e-10 and item["minimum_concentration"] >= -1e-13 for item in pulse["trace"]),
    }
    report["elapsed_seconds"] = time.perf_counter() - started
    (output / "analysis.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    np.savez_compressed(output / "fields.npz", **archive)
    plot_report(report, archive, output)
    return report


def _slice(field, mask):
    values = np.full(mask.shape, np.nan)
    values[mask] = field
    return np.ma.masked_invalid(values[:, :, mask.shape[2] // 2].T)


def plot_report(report, archive, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    colors = ("#237969", "#bb6d41")
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for (name, family), color in zip(report["families"].items(), colors):
        rows = family["rows"]
        ns = np.array([row["n"] for row in rows])
        errors = np.maximum([row["l2_error"] for row in rows], 1e-16)
        label = name.replace("_", " ")
        axes[0].loglog(ns, errors, "o-", color=color, label=label)
        axes[1].semilogy(ns, [max(row["relative_mass_drift"], 1e-16) for row in rows], "o-", color=color, label=label)
    axes[0].loglog(ns, errors[0] * (ns[0] / ns)**2, "k:", label="Second-order reference")
    axes[0].set(title="Same smooth continuum field", xlabel="Compartments per bounding-box axis", ylabel="Volume-weighted RMS error")
    axes[0].legend(fontsize=8)
    axes[1].axhline(1e-10, color="gray", linestyle=":", label="Acceptance threshold")
    axes[1].set(title="Pure diffusion conserves amount", xlabel="Compartments per bounding-box axis", ylabel="Relative mass drift (zeros shown at 1e-16)", ylim=(3e-17, 1e-9))
    axes[1].legend(fontsize=8)
    fig.suptitle("Diffusion convergence on a 3D L-shaped domain")
    fig.savefig(output / "convergence.png", dpi=160)
    plt.close(fig)

    n = report["parameters"]["resolutions"][-1]
    prefix = f"graded_l_prism_n{n}"
    mask = archive[f"{prefix}_mask"]
    x, y = archive[f"{prefix}_x_edges"], archive[f"{prefix}_y_edges"]
    keys = ("initial", "exact", "numerical")
    low = min(float(archive[f"{prefix}_{key}"].min()) for key in keys)
    high = max(float(archive[f"{prefix}_{key}"].max()) for key in keys)
    fig, axes = plt.subplots(1, 4, figsize=(13, 3.8), constrained_layout=True)
    for ax, key, title in zip(axes[:3], keys, ("Initial averages", "Exact continuum averages", "Conservative diffusion")):
        picture = ax.pcolormesh(x, y, _slice(archive[f"{prefix}_{key}"], mask), shading="flat", cmap="viridis", vmin=low, vmax=high, rasterized=True)
        ax.set(title=title, xlabel="x", ylabel="y", aspect="equal")
    fig.colorbar(picture, ax=axes[:3], shrink=.7, label="Concentration (shared scale)")
    error = archive[f"{prefix}_error"]
    bound = max(float(np.max(np.abs(error))), 1e-16)
    picture = axes[3].pcolormesh(x, y, _slice(error, mask), shading="flat", cmap="RdBu_r", vmin=-bound, vmax=bound, rasterized=True)
    axes[3].set(title="Numerical minus exact", xlabel="x", ylabel="y", aspect="equal")
    fig.colorbar(picture, ax=axes[3], shrink=.7, label="Concentration error")
    z = archive[f"{prefix}_z_edges"]
    fig.suptitle(f'Graded L-prism, n={n}; z-slab [{z[n//2]:.3f}, {z[n//2+1]:.3f}]; t={report["parameters"]["duration"]:g}')
    fig.savefig(output / "fields.png", dpi=160)
    plt.close(fig)

    # Boundary faces on the modest pulse mesh make the actual 3D geometry and
    # grading visible. The other panels exercise transport between its arms.
    edges = [archive[f"pulse_{axis}_edges"] for axis in "xyz"]
    mask = archive["pulse_mask"]
    faces = []
    for index in np.argwhere(mask):
        for axis in range(3):
            others = [value for value in range(3) if value != axis]
            for side in (0, 1):
                neighbor = index.copy()
                neighbor[axis] += 2 * side - 1
                if np.all(neighbor >= 0) and np.all(neighbor < np.array(mask.shape)) and mask[tuple(neighbor)]:
                    continue
                face = []
                for a, b in ((0, 0), (1, 0), (1, 1), (0, 1)):
                    vertex = np.zeros(3)
                    vertex[axis] = edges[axis][index[axis] + side]
                    vertex[others[0]] = edges[others[0]][index[others[0]] + a]
                    vertex[others[1]] = edges[others[1]][index[others[1]] + b]
                    face.append(vertex)
                faces.append(face)
    fig = plt.figure(figsize=(13, 4.5))
    grid = fig.add_gridspec(1, 3, left=.065, right=.9, bottom=.12, top=.84, wspace=.32)
    ax = fig.add_subplot(grid[0, 0], projection="3d")
    artist = Poly3DCollection(faces, facecolors="#83b9aa", edgecolors="#365c56", linewidths=.35, alpha=1.)
    ax.add_collection3d(artist)
    length = report["parameters"]["length"]
    ax.set(xlim=(0, length), ylim=(0, length), zlim=(0, length), xlabel="x", ylabel="y", zlabel="z", title="Graded 3D no-flux domain")
    ax.set_box_aspect((1, 1, 1), zoom=.86)
    ax.tick_params(labelsize=8, pad=0)
    for set_ticks, set_labels in ((ax.set_xticks, ax.set_xticklabels),
                                 (ax.set_yticks, ax.set_yticklabels),
                                 (ax.set_zticks, ax.set_zticklabels)):
        set_ticks([0, length / 2, length])
        set_labels(["0", f"{length / 2:g}", f"{length:g}"])
    ax.xaxis.labelpad = ax.yaxis.labelpad = ax.zaxis.labelpad = 0
    ax.title.set_fontsize(11)
    ax.view_init(elev=24, azim=40)
    panels = []
    for position, key, title in ((2, "pulse_0", "Pulse starts in one arm"), (3, "pulse_0.5", "Pulse after Dt/L² = 0.5")):
        ax = fig.add_subplot(grid[0, position - 1])
        picture = ax.pcolormesh(edges[0], edges[1], _slice(archive[key], mask), shading="flat", cmap="viridis", vmin=0, vmax=1)
        ax.set(title=title, xlabel="x", ylabel="y", aspect="equal")
        ax.title.set_fontsize(11)
        panels.append(ax)
    bar_axis = fig.add_axes([.93, .23, .015, .48])
    fig.colorbar(picture, cax=bar_axis, label="Concentration (shared scale)")
    fig.suptitle("Prescribed verification geometry and an independent connectivity probe")
    fig.savefig(output / "mesh.png", dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--resolutions", nargs="+", type=int, default=[4, 8, 16, 32])
    parser.add_argument("--length", type=float, default=1.)
    parser.add_argument("--grading", type=float, default=.35)
    parser.add_argument("--diffusivity", type=float, default=.02)
    parser.add_argument("--duration", type=float, default=1.)
    args = parser.parse_args()
    try:
        result = run_benchmark(**vars(args))
    except (ValueError, FileExistsError) as error:
        parser.error(str(error))
    print(json.dumps(result["checks"], indent=2), flush=True)
    if not all(result["checks"].values()):
        raise SystemExit("Some benchmark criteria were not met; inspect analysis.json.")


if __name__ == "__main__":
    main()
