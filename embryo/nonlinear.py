"""Long nonlinear GM patterns on a fixed L-prism: space, time, and persistence.

A prescribed continuous perturbation is exactly averaged onto every mesh.
Diffusion is implicit (SBDF2), reactions extrapolated explicitly. This method
is not unconditionally positivity-preserving: invalid states abort, never clip.
"""

import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu

from .irregular import graded_edges, l_prism_mask, _positive
from .irregular_signaling import prism_spectrum, stability_report
from .signaling import gm_reaction
from .transport import masked_cartesian_transport


class PrismResolvent:
    """Exact product-domain solve of (alpha I - dt D Delta)c = rhs.

    Transform only the z coordinate to its capacity-weighted eigenbasis, solve
    one sparse 2D shifted system per z mode, and transform back. This is an
    algebraic factorization of the full 3D solve, not dimensional time splitting.
    """

    def __init__(self, spectrum, dt, diffusivity, alpha):
        edges, mask = spectrum.edges, spectrum.mask
        xy = masked_cartesian_transport(edges=(edges[0], edges[1], np.array([0., 1.])),
                                        mask=mask[:, :, :1])
        self.nz = len(spectrum.z_values)
        self.root_volume = np.sqrt(spectrum.transport.volumes).reshape(-1, self.nz)
        self.vectors = spectrum.z_vectors
        identity = sparse.eye(len(xy.volumes), format="csc")
        self.factors = [splu(((alpha + dt * diffusivity * value) * identity
                             + dt * diffusivity * xy.symmetric).tocsc())
                        for value in spectrum.z_values]

    def solve(self, rhs):
        transformed = (self.root_volume * rhs.reshape(-1, self.nz)) @ self.vectors
        for j, factor in enumerate(self.factors):
            transformed[:, j] = factor.solve(transformed[:, j])
        return ((transformed @ self.vectors.T) / self.root_volume).ravel()


class IMEXGM:
    """Fixed-step SBDF2 with a backward-Euler/explicit-reaction startup.

    One first-order startup has O(dt^2) error, consistent with global order two.
    The caller must validate time accuracy and positivity for the intended run.
    Optional zero/custom reactions support independent integrator tests.
    """

    def __init__(self, transport, a, h, dt, beta=2., da=.02, dh=.4,
                 spectrum=None, reaction=None):
        self.dt = _positive("dt", dt)
        self.beta = _positive("beta", beta)
        self.da, self.dh = _positive("da", da), _positive("dh", dh)
        if spectrum is not None and spectrum.transport is not transport:
            raise ValueError("spectrum must belong to the same transport object")
        if np.iscomplexobj(a) or np.iscomplexobj(h):
            raise ValueError("concentrations must be real")
        self.a, self.h = np.asarray(a, dtype=float).copy(), np.asarray(h, dtype=float).copy()
        if self.a.shape != transport.volumes.shape or self.h.shape != self.a.shape:
            raise ValueError("concentrations must match transport size")
        self._check(self.a, self.h)
        self.transport, self.spectrum = transport, spectrum
        self.reaction = reaction if reaction is not None else lambda a, h: gm_reaction(a, h, self.beta)
        self.previous = self.previous_reaction = None
        self.steps = 0
        self.max_amount_balance_residual = 0.
        self.minimum_a, self.minimum_h = float(self.a.min()), float(self.h.min())
        self.solvers = None

    @staticmethod
    def _check(a, h):
        if not np.isfinite(a).all() or not np.isfinite(h).all() or np.any(a <= 0) or np.any(h <= 0):
            raise FloatingPointError("SBDF2 concentrations lost positivity or finiteness; reduce dt")

    def _solvers(self, alpha):
        if self.spectrum is not None:
            return [PrismResolvent(self.spectrum, self.dt, d, alpha) for d in (self.da, self.dh)]
        identity = sparse.eye(len(self.a), format="csc")
        return [splu((alpha * identity - self.dt * d * self.transport.delta).tocsc())
                for d in (self.da, self.dh)]

    def step(self):
        with np.errstate(over="raise", divide="raise", invalid="raise"):
            reaction = self.reaction(self.a, self.h)
            if self.previous is None:
                solvers = self._solvers(1.)
                base = (self.a, self.h)
                source = reaction
                alpha = 1.
            else:
                if self.solvers is None:
                    self.solvers = self._solvers(1.5)
                solvers = self.solvers
                base = tuple(2 * c - .5 * p for c, p in zip((self.a, self.h), self.previous))
                source = tuple(2 * r - p for r, p in zip(reaction, self.previous_reaction))
                alpha = 1.5
            result = tuple(solver.solve(c + self.dt * r) for solver, c, r in zip(solvers, base, source))
            self._check(*result)
            for c, b, r in zip(result, base, source):
                # Reaction changes amount; only diffusion cancels in this balance.
                residual = abs(float(self.transport.volumes @ (alpha * c - b - self.dt * r)))
                scale = max(float(self.transport.volumes @ np.abs(c)), 1e-30)
                self.max_amount_balance_residual = max(self.max_amount_balance_residual, residual / scale)
            self.previous, self.previous_reaction = (self.a, self.h), reaction
            self.a, self.h = result
            self.steps += 1
            self.minimum_a = min(self.minimum_a, float(self.a.min()))
            self.minimum_h = min(self.minimum_h, float(self.h.min()))
        return self.a, self.h


def initial_fields(edges, mask, amplitude=.01):
    """Exact averages of ONE prescribed continuous perturbation on every mesh.

    q = cos(pi*x/L)cos(pi*z/L) + .6cos(pi*y/L)
        + .3cos(2pi*x/L)cos(pi*y/L) - .4/pi.
    Its L-domain mean is zero analytically. No per-mesh normalization or random
    draws are used. q need not satisfy the initial Neumann compatibility at the
    reentrant walls; the finite-volume evolution imposes no flux for t>0.
    """
    amplitude = _positive("amplitude", amplitude)
    if amplitude >= .1:
        raise ValueError("initial amplitude must be below 0.1")
    length = float(edges[0][-1])
    q = np.full(mask.shape, -.4 / np.pi)
    for coefficient, mode in ((1., (1, 0, 1)), (.6, (0, 1, 0)), (.3, (2, 1, 0))):
        averages = [np.cos(np.pi * m * (axis[:-1] + axis[1:]) / (2 * length))
                    * np.sinc(m * np.diff(axis) / (2 * length)) for axis, m in zip(edges, mode)]
        q += coefficient * averages[0][:, None, None] * averages[1][None, :, None] * averages[2][None, None, :]
    q = q[mask]
    return 1 + amplitude * q, 1 - .7 * amplitude * q


def restrict_field(field, fine_volumes, fine_n, coarse_n):
    """Volume-average a nested graded L-prism field, preserving its amount.

    Both meshes must use the same length and grading, with fine_n/coarse_n an
    integer. This is restriction of whole rectangular compartments, not point
    interpolation or a general remesher.
    """
    if type(fine_n) is not int or type(coarse_n) is not int or fine_n < coarse_n or fine_n % coarse_n:
        raise ValueError("restriction requires integer-factor nested resolutions")
    fine_mask, coarse_mask = l_prism_mask(fine_n), l_prism_mask(coarse_n)
    field, volumes = np.asarray(field), np.asarray(fine_volumes)
    if field.shape != (int(fine_mask.sum()),) or volumes.shape != field.shape or np.any(volumes <= 0):
        raise ValueError("field and positive volumes must match the fine mesh")
    if not np.isfinite(field).all() or not np.isfinite(volumes).all():
        raise ValueError("field and volumes must be finite")
    ratio = fine_n // coarse_n
    mass_cube, volume_cube = np.zeros(fine_mask.shape), np.zeros(fine_mask.shape)
    mass_cube[fine_mask], volume_cube[fine_mask] = field * volumes, volumes
    shape = (coarse_n, ratio, coarse_n, ratio, coarse_n, ratio)
    mass = mass_cube.reshape(shape).sum(axis=(1, 3, 5))[coarse_mask]
    capacity = volume_cube.reshape(shape).sum(axis=(1, 3, 5))[coarse_mask]
    return mass / capacity


def field_difference(a, h, reference_a, reference_h, volumes):
    weights = volumes / volumes.sum()
    ca, ch = reference_a - weights @ reference_a, reference_h - weights @ reference_h
    denominator = float(np.sqrt(weights @ (ca**2 + ch**2)))
    absolute = float(np.sqrt(weights @ ((a - reference_a)**2 + (h - reference_h)**2)))
    centered_a = a - weights @ a
    norm = np.sqrt((weights @ centered_a**2) * (weights @ ca**2))
    return {"absolute_joint_rms": absolute,
            "relative_joint_rms": absolute / denominator if denominator > 1e-12 else None,
            "activator_correlation": float(np.clip((weights @ (centered_a * ca)) / norm, -1, 1)) if norm > 1e-24 else None}


def _steps(value, dt, name):
    count = int(round(value / dt))
    if count < 1 or not np.isclose(count * dt, value, rtol=1e-10, atol=1e-12):
        raise ValueError(f"{name} must be a positive integer multiple of dt")
    return count


def run_trajectory(spectrum, duration=100., dt=.05, sample_interval=1., amplitude=.01,
                   beta=2., da=.02, dh=.4):
    duration, dt = _positive("duration", duration), _positive("dt", dt)
    steps = _steps(duration, dt, "duration")
    save_every = _steps(_positive("sample_interval", sample_interval), dt, "sample interval")
    a, h = initial_fields(spectrum.edges, spectrum.mask, amplitude)
    engine = IMEXGM(spectrum.transport, a, h, dt, beta, da, dh, spectrum=spectrum)
    weights = spectrum.transport.volumes / spectrum.transport.volumes.sum()
    frame_steps = {int(round(f * steps)) for f in (0, .2, .4, .6, .75, .8, .9, 1.)}
    late_start = int(np.ceil(.75 * steps))
    anchor = None
    records, times, snapshots_a, snapshots_h = [], [], [], []
    maximum_late_change = None
    started = time.perf_counter()
    for step in range(steps + 1):
        if step == late_start:
            anchor = (engine.a.copy(), engine.h.copy())
        if step % save_every == 0 or step in (late_start, steps):
            mean_a, mean_h = float(weights @ engine.a), float(weights @ engine.h)
            sd = float(np.sqrt(weights @ (engine.a - mean_a)**2))
            record = {"time": step * dt, "mean_a": mean_a, "mean_h": mean_h,
                      "std_a": sd, "std_h": float(np.sqrt(weights @ (engine.h - mean_h)**2)),
                      "max_a": float(engine.a.max()), "min_a": float(engine.a.min()),
                      "min_h": float(engine.h.min())}
            if anchor is not None:
                difference = field_difference(engine.a, engine.h, *anchor, spectrum.transport.volumes)
                record["relative_change_from_late_anchor"] = difference["relative_joint_rms"]
                if difference["relative_joint_rms"] is not None:
                    maximum_late_change = max(maximum_late_change or 0., difference["relative_joint_rms"])
            records.append(record)
        if step in frame_steps:
            times.append(step * dt)
            snapshots_a.append(engine.a.copy())
            snapshots_h.append(engine.h.copy())
        if step < steps:
            engine.step()
    late = [r for r in records if r["time"] >= late_start * dt - 1e-10]
    crossing = next((i for i, row in enumerate(records) if row["std_a"] >= .1), None)
    onset = None
    if crossing is not None:
        if crossing == 0:
            onset = records[0]["time"]
        else:
            lo, hi = records[crossing - 1], records[crossing]
            onset = lo["time"] + (.1 - lo["std_a"]) / (hi["std_a"] - lo["std_a"]) * (hi["time"] - lo["time"])
    report = {"n": spectrum.mask.shape[0], "dt": dt, "steps": steps, "duration": duration,
              "beta": beta, "da": da, "dh": dh, "sample_interval": sample_interval,
              "late_window_start": late_start * dt, "late_samples": len(late),
              "onset_std_0_1": onset, "final_std_a": records[-1]["std_a"],
              "final_mean_a": records[-1]["mean_a"], "minimum_late_std_a": min(r["std_a"] for r in late),
              "maximum_late_relative_field_change": maximum_late_change,
              "minimum_a_all_steps": engine.minimum_a, "minimum_h_all_steps": engine.minimum_h,
              "max_discrete_amount_balance_residual": engine.max_amount_balance_residual,
              "elapsed_seconds": time.perf_counter() - started, "history": records}
    fields = {"time": np.array(times), "a": np.array(snapshots_a), "h": np.array(snapshots_h)}
    return report, fields


def compare_trajectories(coarse, fine, coarse_volumes, fine_volumes, coarse_n, fine_n):
    comparisons = []
    for index, t in enumerate(coarse["time"]):
        matches = np.flatnonzero(np.isclose(fine["time"], t, rtol=0, atol=1e-10))
        if not len(matches):
            continue
        j = int(matches[0])
        fa, fh = fine["a"][j], fine["h"][j]
        if fine_n != coarse_n:
            fa = restrict_field(fa, fine_volumes, fine_n, coarse_n)
            fh = restrict_field(fh, fine_volumes, fine_n, coarse_n)
        comparisons.append({"time": float(t), **field_difference(coarse["a"][index], coarse["h"][index], fa, fh, coarse_volumes)})
    if not comparisons or not np.isclose(comparisons[-1]["time"], coarse["time"][-1]):
        raise ValueError("trajectories must share a final comparison time")
    errors = [r["relative_joint_rms"] for r in comparisons if r["relative_joint_rms"] is not None]
    return {"coarse_n": coarse_n, "fine_n": fine_n, "snapshots": comparisons,
            "maximum_relative_joint_rms": max(errors) if errors else None, "final": comparisons[-1]}


def run_benchmark(output, resolutions=(8, 16, 32), duration=100., dt=.05,
                  length=1., grading=.35, amplitude=.01, sample_interval=1.):
    output = Path(output)
    if output.exists():
        raise FileExistsError("output directory already exists; choose a new path")
    if not isinstance(resolutions, (list, tuple)) or len(resolutions) < 3 or any(
            type(n) is not int or n < 4 or n > 32 or n % 2 for n in resolutions):
        raise ValueError("supply at least three even resolutions from 4 to 32")
    if any(b <= a or b % a for a, b in zip(resolutions, resolutions[1:])):
        raise ValueError("resolutions must increase by integer factors for conservative restriction")
    length, duration, dt = _positive("length", length), _positive("duration", duration), _positive("dt", dt)
    graded_edges(resolutions[0], length, grading)
    initial_fields(graded_edges(resolutions[0], length, grading), l_prism_mask(resolutions[0]), amplitude)
    _steps(duration, dt, "duration")
    _steps(_positive("sample interval", sample_interval), dt, "sample interval")
    if duration < 4 * sample_interval:
        raise ValueError("duration must span at least four sample intervals")
    output.mkdir(parents=True)
    report = {"schema": 1, "parameters": {"resolutions": list(resolutions), "duration": duration,
              "dt": dt, "temporal_steps": [dt, dt / 2, dt / 4], "length": length, "grading": grading,
              "amplitude": amplitude, "sample_interval": sample_interval, "beta": 2., "da": .02, "dh": .4},
              "scope": "One prescribed physical perturbation on a fixed graded L-prism. Finite-window "
                       "nonlinear persistence and numerical refinement, not fate/shape or ensemble validation.",
              "acceptance": {"late_fraction": .25, "persistent_std_min": .1, "late_field_change_max": .01,
                             "fine_spatial_final_relative_error_max": .05, "fine_spatial_correlation_min": .99,
                             "fine_temporal_max_relative_error_max": .005, "stable_control_late_std_max": 1e-6,
                             "amount_balance_residual_max": 1e-9},
              "preflight": [], "spatial_runs": [], "temporal_runs": []}
    archive, spatial_fields, spectra = {}, {}, {}

    def save_fields(prefix, fields):
        archive.update({f"{prefix}_{key}": val for key, val in fields.items()})

    def save_report():
        (output / "analysis.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")

    # All actual spectra are checked before any long nonlinear trajectories.
    for n in resolutions:
        spectra[n] = prism_spectrum(n, length, grading)
        spectrum = spectra[n]
        linear = stability_report(spectrum)
        report["preflight"].append({"n": n, "unstable_mode_count": linear["unstable_mode_count"],
                                    "maximum_spatial_growth": linear["maximum_spatial_growth"],
                                    "local_stable": linear["local_stable"]})
        for axis, edges in zip("xyz", spectrum.edges):
            archive[f"n{n}_{axis}_edges"] = edges
        archive[f"n{n}_mask"], archive[f"n{n}_volumes"] = spectrum.mask, spectrum.transport.volumes
    save_report()
    try:
        for n in resolutions:
            print(f"Nonlinear spatial run n={n}, dt={dt}, t={duration}", flush=True)
            row, fields = run_trajectory(spectra[n], duration, dt, sample_interval, amplitude)
            report["spatial_runs"].append(row)
            spatial_fields[n] = fields
            save_fields(f"n{n}", fields)
            change = row["maximum_late_relative_field_change"]
            change_text = f"{change:.3g}" if change is not None else "undefined (uniform field)"
            print(f'  final SD={row["final_std_a"]:.6g}; late relative change={change_text}', flush=True)
            save_report()
        finest = resolutions[-1]
        temporal_fields = [spatial_fields[finest]]
        report["temporal_runs"].append(report["spatial_runs"][-1])
        for factor in (2, 4):
            print(f"Temporal control n={finest}, dt={dt/factor}", flush=True)
            row, fields = run_trajectory(spectra[finest], duration, dt / factor, sample_interval, amplitude)
            report["temporal_runs"].append(row)
            temporal_fields.append(fields)
            save_fields(f"temporal_{factor}", fields)
            save_report()
        control_n = resolutions[-2]
        report["stable_control_preflight"] = {"n": control_n, "unstable_mode_count": stability_report(spectra[control_n], dh=.02)["unstable_mode_count"]}
        print(f"Equal-diffusivity stable control n={control_n}", flush=True)
        control, fields = run_trajectory(spectra[control_n], duration, dt, sample_interval, amplitude, dh=.02)
        report["stable_control"] = control
        save_fields("stable_control", fields)
    except FloatingPointError as error:
        report["numerical_failure"] = str(error)
        save_report()
        np.savez_compressed(output / "fields.npz", **archive)
        raise
    report["spatial_comparisons"] = [compare_trajectories(spatial_fields[a], spatial_fields[b],
            spectra[a].transport.volumes, spectra[b].transport.volumes, a, b) for a, b in zip(resolutions, resolutions[1:])]
    volume = spectra[finest].transport.volumes
    report["temporal_comparisons"] = [compare_trajectories(a, b, volume, volume, finest, finest)
                                       for a, b in zip(temporal_fields, temporal_fields[1:])]
    spatial = report["spatial_comparisons"]
    temporal = report["temporal_comparisons"]
    e0, e1 = (item["maximum_relative_joint_rms"] for item in temporal)
    report["observed_temporal_order"] = float(np.log2(e0 / e1)) if e0 is not None and e1 is not None and min(e0, e1) > 1e-10 else None
    # Restriction defines like-for-like physical averages, with no shifting,
    # rotations, sign flips, or best-fit registration to improve agreement.
    finest_runs = [report["spatial_runs"][-1], *report["temporal_runs"][1:]]
    all_runs = [*report["spatial_runs"], *report["temporal_runs"][1:], control]
    spatial_error = spatial[-1]["final"]["relative_joint_rms"]
    previous_spatial_error = spatial[-2]["final"]["relative_joint_rms"]
    correlation = spatial[-1]["final"]["activator_correlation"]
    report["checks"] = {
        "initial_fields_agree_under_restriction": all(row["snapshots"][0]["absolute_joint_rms"] < 1e-12 for row in spatial),
        "positive_concentrations_all_steps": all(min(row["minimum_a_all_steps"], row["minimum_h_all_steps"]) > 0 for row in all_runs),
        "discrete_reaction_amount_balance": all(row["max_discrete_amount_balance_residual"] < 1e-9 for row in all_runs),
        "persistent_fine_mesh_contrast": all(row["minimum_late_std_a"] >= .1 for row in finest_runs),
        "fine_mesh_late_pattern_stationary": all(row["maximum_late_relative_field_change"] is not None and row["maximum_late_relative_field_change"] <= .01 for row in finest_runs),
        "spatial_final_difference_decreases": spatial_error is not None and previous_spatial_error is not None and spatial_error < previous_spatial_error,
        "fine_spatial_final_error_below_5_percent": spatial_error is not None and spatial_error <= .05,
        "fine_spatial_final_correlation_above_0_99": correlation is not None and correlation >= .99,
        "temporal_error_decreases": e0 is not None and e1 is not None and (e1 < e0 or max(e0, e1) < 1e-10),
        "fine_temporal_error_below_0_5_percent": e1 is not None and e1 <= .005,
        "equal_diffusivity_control_loses_contrast": max(row["std_a"] for row in control["history"] if row["time"] >= control["late_window_start"]) < 1e-6,
    }
    save_report()
    np.savez_compressed(output / "fields.npz", **archive)
    plot_report(report, archive, output)
    return report


def plot_report(report, archive, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    for row in report["spatial_runs"]:
        axes[0].plot([r["time"] for r in row["history"]], [r["std_a"] for r in row["history"]], label=f'n={row["n"]}')
    control = report["stable_control"]
    axes[0].plot([r["time"] for r in control["history"]], [r["std_a"] for r in control["history"]], "k--", label="Equal diffusion control")
    axes[0].axvspan(.75 * report["parameters"]["duration"], report["parameters"]["duration"], color="gray", alpha=.12)
    axes[0].set(title="Growth and finite-window persistence", xlabel="Time", ylabel="Volume-weighted activator SD")
    axes[0].legend(fontsize=8)
    for item in report["spatial_comparisons"]:
        axes[1].plot([r["time"] for r in item["snapshots"]], [r["relative_joint_rms"] for r in item["snapshots"]], "o-", label=f'n={item["coarse_n"]} vs {item["fine_n"]}')
    axes[1].set(title="Conservative spatial comparison", xlabel="Time", ylabel="Relative two-species field error")
    axes[1].legend(fontsize=8)
    for i, item in enumerate(report["temporal_comparisons"]):
        axes[2].plot([r["time"] for r in item["snapshots"]], [r["relative_joint_rms"] for r in item["snapshots"]], "o-", label=f'dt={report["parameters"]["dt"]/2**i:g} vs {report["parameters"]["dt"]/2**(i+1):g}')
    axes[2].set(title="Time refinement on the finest mesh", xlabel="Time", ylabel="Relative two-species field error")
    axes[2].legend(fontsize=8)
    fig.suptitle("Nonlinear Gierer-Meinhardt signaling on a fixed L-shaped domain")
    fig.savefig(output / "convergence.png", dpi=160)
    plt.close(fig)
    finest = report["parameters"]["resolutions"][-1]
    mask = archive[f"n{finest}_mask"]
    a = archive[f"n{finest}_a"]
    times = archive[f"n{finest}_time"]
    chosen = sorted(set([0, min(2, len(times) - 1), len(times) - 1]))
    x, y, z = (archive[f"n{finest}_{axis}_edges"] for axis in "xyz")
    # Use a fixed near-boundary plane, not a separately selected peak at each time.
    slab = 0
    fig, axes = plt.subplots(1, len(chosen), figsize=(4.2 * len(chosen), 4.4), constrained_layout=True, squeeze=False)
    low, high = float(a.min()), float(a.max())
    for ax, index in zip(axes[0], chosen):
        cube = np.full(mask.shape, np.nan)
        cube[mask] = a[index]
        picture = ax.pcolormesh(x, y, np.ma.masked_invalid(cube[:, :, slab].T), cmap="viridis", vmin=low, vmax=high, shading="flat")
        ax.set(aspect="equal", xlabel="x", ylabel="y", title=f't={times[index]:g}')
    fig.colorbar(picture, ax=axes[0], shrink=.75, label="Activator concentration (shared scale)")
    fig.suptitle(f'Same physical perturbation; n={finest}, z-slab [{z[slab]:.3f}, {z[slab+1]:.3f}]')
    fig.savefig(output / "patterns.png", dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--resolutions", nargs="+", type=int, default=[8, 16, 32])
    for name, default in (("duration", 100.), ("dt", .05), ("length", 1.), ("grading", .35), ("amplitude", .01), ("sample-interval", 1.)):
        parser.add_argument("--" + name, type=float, default=default)
    args = parser.parse_args()
    try:
        report = run_benchmark(**vars(args))
    except (ValueError, FileExistsError) as error:
        parser.error(str(error))
    print(json.dumps(report["checks"], indent=2), flush=True)
    if not all(report["checks"].values()):
        raise SystemExit("Some nonlinear criteria were not met; inspect analysis.json.")


if __name__ == "__main__":
    main()
