"""Run an experiment and export an offline viewer, metrics, and full checkpoint."""

import argparse
import csv
from dataclasses import asdict
import json
from pathlib import Path
import time

from .model import Config, Simulation
from .signaling import normalized_graph, stability
from .graph_analysis import cycle


def export_plot(frames, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    fig = plt.figure(figsize=(13, 7), facecolor="#101827")
    counts = [1, 4, 8, frames[-1]["metrics"]["cells"]]
    selected = []
    for count in counts[:-1]:
        selected.append(next((f for f in frames if f["metrics"]["cells"] >= count), frames[-1]))
    selected.append(frames[-1])
    cmap = plt.get_cmap("coolwarm")
    for k, frame in enumerate(selected):
        ax = fig.add_subplot(2, 4, k + 1, projection="3d", facecolor="#101827")
        for cell in frame["cells"]:
            p = np.asarray(cell["points"])
            if len(p):
                color = cmap(np.clip(0.5 + cell["fate"] / 2.4, 0, 1))
                ax.scatter(*p.T, color=color, s=2, alpha=0.7, linewidths=0)
        ax.set(xlim=(-1.2, 1.2), ylim=(-1.2, 1.2), zlim=(-1.2, 1.2))
        ax.set_box_aspect((1, 1, 1))
        ax.set_axis_off()
        count = frame["metrics"]["cells"]
        ax.set_title(f't = {frame["metrics"]["time"]:.2f} · {count} cell' + ("s" if count != 1 else ""), color="white")
    ts = [f["metrics"]["time"] for f in frames]
    ax = fig.add_subplot(2, 2, 3, facecolor="#101827")
    for key, color, label in [("fate_a", "#ed8d77", "A tendency"), ("fate_b", "#7ab9ef", "B tendency"), ("uncommitted", "#d3d7df", "Uncommitted")]:
        ax.plot(ts, [f["metrics"][key] for f in frames], color=color, label=label)
    ax.set_ylabel("Cells", color="white")
    ax.legend(facecolor="#101827", labelcolor="white", frameon=False)
    ax2 = fig.add_subplot(2, 2, 4, facecolor="#101827")
    ax2.plot(ts, [f["metrics"]["axis_ratio"] for f in frames], color="#a7e1c0")
    ax2.set_ylabel("Principal axis ratio", color="white")
    for plot in (ax, ax2):
        plot.tick_params(colors="#b5c1d4")
        plot.set_xlabel("Dimensionless time", color="#b5c1d4")
        for spine in plot.spines.values():
            spine.set_color("#354258")
    fig.suptitle("Embryogenesis · deformable cells and two fate tendencies", color="white", fontsize=18)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=150, facecolor=fig.get_facecolor())
    plt.close(fig)


def run(config, output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    sim = Simulation(config)
    # Analyze finite reference graphs before advancing any 3D state.
    import numpy as np
    preflight = {}
    for n in (4, 8, 16):
        preflight[f"cycle_{n}"] = stability(cycle(n), config.signal_beta, config.signal_da, config.signal_dh)
        preflight[f"complete_{n}"] = stability(normalized_graph(np.ones((n, n)) - np.eye(n)),
                                               config.signal_beta, config.signal_da, config.signal_dh)
    (output / "graph_preflight.json").write_text(json.dumps(preflight, indent=2) + "\n")
    frames = []
    started = time.monotonic()
    (output / "config.json").write_text(json.dumps(asdict(config), indent=2) + "\n")
    for step in range(config.steps + 1):
        if step % config.save_every == 0 or step == config.steps:
            metrics = sim.metrics()
            frames.append({"metrics": metrics, "cells": sim.surfaces(), "graph": sim.graph_snapshot()})
            if step % (config.save_every * 5) == 0 or step == config.steps:
                print(f'{step:5d}/{config.steps}  t={sim.time:6.2f}  cells={len(sim.phi):2d} '
                      f'dividing={len(sim.divisions)} '
                      f'A/B/U={metrics["fate_a"]}/{metrics["fate_b"]}/{metrics["uncommitted"]} '
                      f'volume error={metrics["relative_volume_error"]:+.2%}', flush=True)
        if step < config.steps:
            sim.step()
    sim.checkpoint(output / "final_state.npz")
    (output / "lineage.json").write_text(json.dumps(sim.lineage, indent=2) + "\n")
    (output / "graph_history.json").write_text(json.dumps([f["graph"] for f in frames], indent=2) + "\n")
    (output / "cleavage_spectra.json").write_text(json.dumps(sim.graph_events, indent=2) + "\n")
    with (output / "metrics.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=frames[0]["metrics"].keys())
        writer.writeheader()
        writer.writerows(frame["metrics"] for frame in frames)
    payload = {"config": asdict(config), "frames": frames}
    text = json.dumps(payload, separators=(",", ":"), allow_nan=False)
    (output / "trajectory.json").write_text(text)
    template = Path(__file__).with_name("viewer.html").read_text()
    (output / "viewer.html").write_text(template.replace("__SIMULATION_DATA__", text))
    export_plot(frames, output / "summary.png")
    diagnostics = {
        "elapsed_seconds": time.monotonic() - started,
        "maximum_absolute_volume_error": max(abs(f["metrics"]["relative_volume_error"]) for f in frames),
        "maximum_boundary_occupancy": max(f["metrics"]["boundary_occupancy"] for f in frames),
        "minimum_radius_grid_cells": min(f["metrics"]["min_radius_grid_cells"] for f in frames),
        "maximum_overdue_divisions": max(f["metrics"]["overdue_divisions"] for f in frames),
        "final": frames[-1]["metrics"],
        "interpretation": "Thresholded fate tendencies; commitment and biological validation are not established.",
    }
    (output / "diagnostics.json").write_text(json.dumps(diagnostics, indent=2) + "\n")
    print(f'Wrote {output / "viewer.html"} ({diagnostics["elapsed_seconds"]:.1f}s)', flush=True)
    return diagnostics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, help="JSON file with Config overrides")
    parser.add_argument("--output", type=Path, default=Path("outputs/demo"), help="New output directory (never overwritten)")
    parser.add_argument("--seed", type=int)
    parser.add_argument("--grid", type=int)
    parser.add_argument("--steps", type=int)
    parser.add_argument("--max-cells", type=int)
    parser.add_argument("--no-feedback", action="store_true", help="Disable fate-dependent mechanics")
    parser.add_argument("--mechanics-only", action="store_true", help="Disable differentiation and its feedback")
    parser.add_argument("--no-signaling", action="store_true", help="Disable Gierer–Meinhardt evolution and partition noise")
    parser.add_argument("--no-polarity", action="store_true", help="Disable apical–basal polarity and its tension contrast")
    args = parser.parse_args()
    values = json.loads(args.config.read_text()) if args.config else {}
    for name in ("seed", "grid", "steps", "max_cells"):
        if getattr(args, name) is not None:
            values[name] = getattr(args, name)
    if args.no_feedback or args.mechanics_only:
        values["feedback"] = False
    if args.mechanics_only:
        values["differentiation"] = False
    if args.no_signaling or args.mechanics_only:
        values["signaling"] = False
    if args.no_polarity or args.mechanics_only:
        values["polarity_enabled"] = False
    try:
        config = Config(**values)
        config.validate()
        if args.output.exists():
            parser.error(f"output directory exists: {args.output}; choose a new path")
        run(config, args.output)
    except (ValueError, TypeError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
