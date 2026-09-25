"""Reproducible frozen-graph, graph-replay, and coupled cleavage-timescale study.

Run capture first to record every mechanical step, then analyze without rerunning
3D mechanics. Graph replay is a prescribed-geometry control, not a coupled run.
"""

import argparse
import csv
from dataclasses import asdict, replace
import json
from pathlib import Path
import time

import numpy as np
from scipy.linalg import expm

from .model import Config, Simulation
from .signaling import normalized_graph, stability, integrate, gm_jacobian
from .graph_analysis import cycle, verify_mode


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def first_crossing(times, amplitudes, threshold):
    """First sampled crossing, not an interpolated or persistent-onset estimate."""
    indices = np.flatnonzero(np.asarray(amplitudes) >= threshold)
    return float(times[indices[0]]) if len(indices) else None


def capture(config, output):
    """Record pre-step contacts and states at every step, including all cleavages."""
    config.validate()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    write_json(output / "config.json", asdict(config))
    write_json(output / "graph_preflight.json", {
        f"cycle_{n}": stability(cycle(n), config.signal_beta, config.signal_da, config.signal_dh)
        for n in (4, 8, 16)
    })
    sim = Simulation(config)
    shape = (config.steps + 1, config.max_cells)
    weights = np.zeros((*shape, config.max_cells))
    counts = np.zeros(shape[0], dtype=int)
    ids = np.full(shape, -1, dtype=int)
    activator, inhibitor = np.zeros(shape), np.zeros(shape)
    metrics = []
    started = time.monotonic()
    for step in range(config.steps + 1):
        n = len(sim.phi)
        weights[step, :n, :n] = sim.contacts()[0]
        counts[step] = n
        ids[step, :n] = sim.ids
        activator[step, :n], inhibitor[step, :n] = sim.activator, sim.inhibitor
        if step % config.save_every == 0 or step == config.steps:
            row = sim.metrics()
            graph = normalized_graph(weights[step, :n, :n], config.graph_contact_cutoff)
            row["maximum_spatial_growth"] = stability(
                graph, config.signal_beta, config.signal_da, config.signal_dh)["maximum_spatial_growth"]
            metrics.append(row)
        if step % 200 == 0 or step == config.steps:
            print(f"{output.name}: t={sim.time:.3f}, cells={n}, "
                  f"std(a)={sim.activator.std():.6g}, elapsed={time.monotonic()-started:.1f}s", flush=True)
        if step < config.steps:
            sim.step()
    np.savez_compressed(output / "trace.npz", time=np.arange(config.steps + 1) * config.dt,
                        counts=counts, ids=ids, weights=weights,
                        activator=activator, inhibitor=inhibitor)
    write_json(output / "cleavage_spectra.json", sim.graph_events)
    write_json(output / "metrics.json", metrics)
    with (output / "metrics.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=metrics[0])
        writer.writeheader()
        writer.writerows(metrics)
    sim.checkpoint(output / "final_state.npz")
    write_json(output / "run.json", {"elapsed_seconds": time.monotonic() - started,
                                    "final": metrics[-1], "trace_schema": 1})


def load_capture(directory):
    directory = Path(directory)
    config = Config(**json.loads((directory / "config.json").read_text()))
    with np.load(directory / "trace.npz", allow_pickle=False) as data:
        trace = {key: data[key] for key in data.files}
    events = json.loads((directory / "cleavage_spectra.json").read_text())
    return config, trace, events


def activity_std(trace, name="activator"):
    return np.array([row[:n].std() for row, n in zip(trace[name], trace["counts"])])


def apply_inheritance(values, event, species):
    """Reuse fractional partition draws, not additive jumps from another state."""
    p = np.asarray(event["prolongation"])
    original = p @ np.asarray(event[f"{species}_before"])
    jump = np.asarray(event[f"{species}_partition_jump"])
    if np.any((original == 0) & (jump != 0)):
        raise ValueError("nonzero partition jump from a zero parent activity")
    fractional = np.divide(jump, original, out=np.zeros_like(jump), where=original != 0)
    return (p @ values) * (1 + fractional)


def replay(config, trace, events, precleavage_factor=1.0):
    """Stretch all graph dwell times through last cleavage; later times unchanged.

    The same contacts, ID order, lineage and fractional partition draws are used.
    No mechanical response to altered signals is allowed. Holding factor=1 must
    reconstruct recorded signaling to numerical roundoff.
    """
    if not np.isfinite(precleavage_factor) or precleavage_factor <= 0:
        raise ValueError("precleavage_factor must be finite and positive")
    by_step = {}
    for event in events:
        step = int(round(event["time"] / config.dt))
        by_step.setdefault(step, []).append(event)
    last_division = max(by_step, default=0)
    a = trace["activator"][0, :trace["counts"][0]].copy()
    b = trace["inhibitor"][0, :len(a)].copy()
    result = {"time": np.zeros(len(trace["time"])), "counts": trace["counts"].copy(),
              "activator": np.zeros_like(trace["activator"]),
              "inhibitor": np.zeros_like(trace["inhibitor"])}
    result["activator"][0, :len(a)] = a
    result["inhibitor"][0, :len(b)] = b
    for step in range(len(trace["time"]) - 1):
        n = int(trace["counts"][step])
        if len(a) != n:
            raise ValueError("trace and cleavage events have inconsistent cell counts")
        graph = normalized_graph(trace["weights"][step, :n, :n], config.graph_contact_cutoff)
        dt = trace["time"][step + 1] - trace["time"][step]
        # Use the configured step to exactly match the coupled integrator.
        if not np.isclose(dt, config.dt):
            raise ValueError("replay requires a trace recorded at every solver step")
        dt = config.dt * (precleavage_factor if step < last_division else 1)
        a, b = integrate(a, b, graph, dt, config.signal_beta, config.signal_da, config.signal_dh)
        for event in by_step.get(step + 1, []):
            a = apply_inheritance(a, event, "activator")
            b = apply_inheritance(b, event, "inhibitor")
        if len(a) != trace["counts"][step + 1]:
            raise ValueError("missing or inconsistent cleavage event")
        result["time"][step + 1] = result["time"][step] + dt
        result["activator"][step + 1, :len(a)] = a
        result["inhibitor"][step + 1, :len(b)] = b
    return result


def frozen(config, weights, activator, inhibitor, duration=60., dt=.05):
    """Nonlinear signals and exact linearized evolution from the same real state."""
    if min(duration, dt) <= 0 or not np.isfinite([duration, dt]).all():
        raise ValueError("duration and dt must be finite and positive")
    graph = normalized_graph(weights, config.graph_contact_cutoff)
    count = len(graph.degree)
    steps = int(np.ceil(duration / dt))
    dt = duration / steps
    j, identity = gm_jacobian(config.signal_beta), np.eye(count)
    operator = np.block([[j[0, 0] * identity + config.signal_da * graph.delta, j[0, 1] * identity],
                         [j[1, 0] * identity, j[1, 1] * identity + config.signal_dh * graph.delta]])
    transition = expm(operator * dt)
    linear = np.r_[np.asarray(activator) - 1, np.asarray(inhibitor) - 1]
    a, b = np.asarray(activator).copy(), np.asarray(inhibitor).copy()
    result = {"time": np.arange(steps + 1) * dt,
              "activator": np.zeros((steps + 1, count)),
              "inhibitor": np.zeros((steps + 1, count)),
              "linear_activator_deviation": np.zeros((steps + 1, count)),
              "linear_inhibitor_deviation": np.zeros((steps + 1, count))}
    for step in range(steps + 1):
        result["activator"][step], result["inhibitor"][step] = a, b
        result["linear_activator_deviation"][step] = linear[:count]
        result["linear_inhibitor_deviation"][step] = linear[count:]
        if step < steps:
            a, b = integrate(a, b, graph, dt, config.signal_beta, config.signal_da, config.signal_dh)
            linear = transition @ linear
    return graph, result


def trace_summary(trace):
    amplitudes = activity_std(trace)
    stages = {}
    for n in (4, 8, 16):
        indices = np.flatnonzero(trace["counts"] == n)
        if len(indices):
            index = indices[0]
            changed = np.flatnonzero(trace["counts"][index:] != n)
            stages[str(n)] = {"time": float(trace["time"][index]),
                              "activator_std": float(amplitudes[index]),
                              "time_until_next_count": float(trace["time"][index + changed[0]] - trace["time"][index])
                                  if len(changed) else None}
    result = {"stages": stages, "final_time": float(trace["time"][-1]),
            "final_activator_std": float(amplitudes[-1]),
            "first_std_0.01": first_crossing(trace["time"], amplitudes, .01),
            "first_std_0.1": first_crossing(trace["time"], amplitudes, .1)}
    if "16" in stages:
        for threshold in ("0.01", "0.1"):
            onset = result[f"first_std_{threshold}"]
            result[f"first_std_{threshold}_relative_to_16_cells"] = (
                onset - stages["16"]["time"] if onset is not None else None)
    return result


def matched_followup(traces):
    """Compare amplitude at equal elapsed time after first reaching 16 cells."""
    origins = {name: float(trace["time"][np.flatnonzero(trace["counts"] == 16)[0]])
               for name, trace in traces.items()}
    age = min(trace["time"][-1] - origins[name] for name, trace in traces.items())
    return {"time_since_16_cells": float(age), "activator_std": {
        name: float(np.interp(origins[name] + age, trace["time"], activity_std(trace)))
        for name, trace in traces.items()}}


def analyze(baseline, slow, output, frozen_duration=120., replay_factor=2., frozen_dt=.05):
    output = Path(output)
    if output.exists():
        raise ValueError("analysis output directory already exists")
    config, trace, events = load_capture(baseline)
    slow_config, slow_trace, _ = load_capture(slow)
    differing = {k for k, v in asdict(config).items() if asdict(slow_config)[k] != v}
    if differing - {"division_interval"}:
        raise ValueError(f"coupled comparison changes more than division interval: {differing}")
    if not np.any(trace["counts"] == 16) or not np.any(slow_trace["counts"] == 16):
        raise ValueError("both coupled runs must reach 16 cells")
    output.mkdir(parents=True)
    replayed = replay(config, trace, events)
    errors = {key: float(np.max(np.abs(replayed[key] - trace[key])))
              for key in ("activator", "inhibitor")}
    if max(errors.values()) > 1e-8:
        raise ValueError(f"native graph replay does not reconstruct coupled signals: {errors}")
    stretched = replay(config, trace, events, replay_factor)
    np.savez_compressed(output / "replay_native.npz", **replayed)
    np.savez_compressed(output / "replay_stretched.npz", **stretched)
    summary = {"baseline_source": str(Path(baseline).resolve()),
               "slow_source": str(Path(slow).resolve()),
               "baseline_config": asdict(config), "slow_config": asdict(slow_config),
               "parameters": {"frozen_duration": frozen_duration, "frozen_dt": frozen_dt,
                              "replay_factor": replay_factor,
                              "slow_cycle_factor": slow_config.division_interval / config.division_interval,
                              "std_thresholds": [.01, .1]},
               "replay_max_absolute_error": errors, "frozen": {},
               "coupled_baseline": trace_summary(trace), "coupled_slow": trace_summary(slow_trace),
               "replay_native": trace_summary(replayed), "replay_stretched": trace_summary(stretched),
               "matched_coupled_followup": matched_followup({"baseline": trace, "slow": slow_trace}),
               "matched_replay_followup": matched_followup({"native": replayed, "stretched": stretched}),
               "interpretation": "Single-seed, coarse-grid timescale screen. Frozen graphs start at exact "
                   "first attainment of each count; replay prescribes geometry and reuses fractional noise. "
                   "Thresholds measure signal amplitude, not fate commitment. Linear continuation outside "
                   "the small-perturbation regime is a mathematical comparison, not a nonlinear forecast."}
    curves = {}
    for n in (4, 8, 16):
        index = int(np.flatnonzero(trace["counts"] == n)[0])
        weights = trace["weights"][index, :n, :n]
        graph, result = frozen(config, weights, trace["activator"][index, :n],
                               trace["inhibitor"][index, :n], frozen_duration, frozen_dt)
        np.savez_compressed(output / f"frozen_{n}.npz", weights=weights, **result)
        report = stability(graph, config.signal_beta, config.signal_da, config.signal_dh)
        rate = report["maximum_spatial_growth"]
        amplitude = result["activator"].std(axis=1)
        summary["frozen"][str(n)] = {
            "source_time": float(trace["time"][index]), "stability": report,
            "linear_mode_verification": verify_mode(graph, config.signal_beta, config.signal_da, config.signal_dh),
            "efold_time": 1 / rate if rate > 0 else None,
            "initial_activator_std": float(amplitude[0]), "final_activator_std": float(amplitude[-1]),
            "first_std_0.01": first_crossing(result["time"], amplitude, .01),
            "first_std_0.1": first_crossing(result["time"], amplitude, .1)}
        curves[n] = result
    for name, directory in (("coupled_baseline", baseline), ("coupled_slow", slow)):
        metrics = json.loads((Path(directory) / "metrics.json").read_text())
        summary[name]["final_metrics"] = metrics[-1]
        summary[name]["maximum_sampled_cell_volume_error"] = max(m["max_cell_volume_error"] for m in metrics)
        summary[name]["maximum_sampled_boundary_occupancy"] = max(m["boundary_occupancy"] for m in metrics)
        summary[name]["maximum_sampled_overdue_divisions"] = max(m["overdue_divisions"] for m in metrics)
    write_json(output / "analysis.json", summary)
    plot_results(curves, trace, slow_trace, replayed, stretched, summary, output / "timescales.png")
    print(json.dumps({k: v for k, v in summary.items() if k in
                     ("frozen", "replay_max_absolute_error", "replay_stretched")}, indent=2), flush=True)
    return summary


def plot_results(curves, baseline, slow, native, stretched, summary, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    for n, curve in curves.items():
        line, = axes[0, 0].semilogy(curve["time"], np.maximum(curve["activator"].std(axis=1), 1e-15),
                                  label=f"{n} cells, nonlinear")
        axes[0, 0].semilogy(curve["time"], np.maximum(curve["linear_activator_deviation"].std(axis=1), 1e-15),
                            "--", color=line.get_color(), alpha=.7)
    axes[0, 0].set(title="Frozen actual graphs: solid nonlinear, dashed linear",
                   xlabel="Time since graph was frozen", ylabel="Activator standard deviation", ylim=(1e-8, 2))
    axes[0, 0].legend(fontsize=8)
    for n, report in summary["frozen"].items():
        check = report["linear_mode_verification"]
        if "fitted" in check:
            axes[0, 1].scatter(check["predicted"], check["fitted"], s=50)
            axes[0, 1].annotate(f"  {n} cells", (check["predicted"], check["fitted"]))
    axes[0, 1].plot([-.3, .3], [-.3, .3], "k--", alpha=.4)
    axes[0, 1].axhline(0, color="gray", linewidth=.5)
    axes[0, 1].axvline(0, color="gray", linewidth=.5)
    axes[0, 1].set(title="Independent small-eigenmode growth check",
                   xlabel="Predicted rate", ylabel="Measured early exponential rate")
    for ax, cases in [(axes[1, 0], [(baseline, "Coupled: normal cycle"), (slow, "Coupled: doubled cycle interval")]),
                      (axes[1, 1], [(native, "Recorded geometry: native timing"),
                                    (stretched, "Same geometry: stretched before final cleavage")])]:
        for trace, label in cases:
            ax.semilogy(trace["time"], np.maximum(activity_std(trace), 1e-15), label=label)
        ax.axhline(.01, color="gray", linestyle=":", linewidth=1)
        ax.axhline(.1, color="gray", linestyle=":", linewidth=1)
        ax.set(xlabel="Time since zygote", ylabel="Activator standard deviation", ylim=(1e-5, 2))
        ax.legend(fontsize=8)
    axes[1, 0].set_title("Fully coupled 3D: longer time and slower cleavage")
    axes[1, 1].set_title("Prescribed-geometry control: identical fractional noise draws")
    fig.suptitle(f"Signaling growth versus cleavage timescale · seed {summary['baseline_config']['seed']} · dimensionless units")
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    record = commands.add_parser("capture", help="Run and record coupled 3D mechanics")
    record.add_argument("--output", type=Path, required=True)
    record.add_argument("--duration", type=float, default=60.)
    record.add_argument("--cycle-factor", type=float, default=1.)
    record.add_argument("--grid", type=int, default=40)
    record.add_argument("--dt", type=float, default=.015)
    record.add_argument("--seed", type=int, default=7)
    analysis = commands.add_parser("analyze", help="Frozen graphs and replay; no 3D stepping")
    analysis.add_argument("--baseline", type=Path, required=True)
    analysis.add_argument("--slow", type=Path, required=True)
    analysis.add_argument("--output", type=Path, required=True)
    analysis.add_argument("--frozen-duration", type=float, default=120.)
    analysis.add_argument("--frozen-dt", type=float, default=.05)
    analysis.add_argument("--replay-factor", type=float, default=2.)
    args = parser.parse_args()
    try:
        if args.command == "capture":
            if min(args.duration, args.dt, args.cycle_factor) <= 0 or not np.isfinite(
                    [args.duration, args.dt, args.cycle_factor]).all():
                raise ValueError("duration, dt, and cycle factor must be finite and positive")
            count = round(args.duration / args.dt)
            if not np.isclose(count * args.dt, args.duration):
                raise ValueError("duration must be a whole number of time steps")
            config = replace(Config(), seed=args.seed, grid=args.grid, dt=args.dt, steps=count,
                             division_interval=Config().division_interval * args.cycle_factor)
            capture(config, args.output)
        else:
            analyze(args.baseline, args.slow, args.output, args.frozen_duration, args.replay_factor, args.frozen_dt)
    except (ValueError, FileExistsError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
