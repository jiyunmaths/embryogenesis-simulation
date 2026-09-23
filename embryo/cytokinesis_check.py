"""Reproducible single-cleavage volume/contact and resolution checks."""

import argparse
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np

from .model import Config, Simulation


def measure(name, model, config):
    sim = model(config)

    def sample():
        return {"time": sim.time, "cells": len(sim.phi),
                "volume_error": float(np.max(abs(sim.volumes() / sim.target - 1))),
                "contact": float(sim.contacts()[0].sum() / 2)}

    trace = [sample()]
    sim.divide(0, [1, 2, 3])  # Same oblique axis in every case; isolate cytokinesis.
    trace.append(sample())
    for _ in range(round(3 / config.dt)):
        sim.step()
        trace.append(sample())
    return {"model": name, "grid": config.grid, "dt": config.dt,
            "abscission": sim.lineage[0]["division"],
            "max_cell_volume_error": max(t["volume_error"] for t in trace),
            "final_volume_error": trace[-1]["volume_error"],
            "final_axis_ratio": sim.metrics()["axis_ratio"],
            "max_contact_increment": max(abs(b["contact"] - a["contact"]) for a, b in zip(trace, trace[1:])),
            "trace": trace}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="New directory")
    parser.add_argument("--baseline", type=Path, help="Optional archived local model.py for instantaneous-division comparison")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output already exists; choose a new directory")
    cases = []
    if args.baseline:
        spec = importlib.util.spec_from_file_location("archived_embryo_model", args.baseline)
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        cases.append(("instantaneous", module.Simulation,
                      module.Config(max_cells=2, differentiation=False, feedback=False)))
    for name, grid, dt in [("progressive", 40, .015), ("progressive_half_dt", 40, .0075),
                           ("progressive_grid48", 48, .015)]:
        cases.append((name, Simulation, Config(grid=grid, dt=dt, max_cells=2,
                                               differentiation=False, feedback=False,
                                               signaling=False, polarity_enabled=False)))
    results = []
    for case in cases:
        result = measure(*case)
        results.append(result)
        print({k: v for k, v in result.items() if k != "trace"}, flush=True)
    args.output.mkdir(parents=True)
    (args.output / "comparison.json").write_text(json.dumps(results, indent=2) + "\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    for result in results:
        t = [row["time"] for row in result["trace"]]
        axes[0].plot(t, [100 * row["volume_error"] for row in result["trace"]], label=result["model"])
        axes[1].plot(t, [row["contact"] for row in result["trace"]])
    axes[0].set_ylabel("Maximum individual volume error (%)")
    axes[1].set_ylabel("Diffuse contact proxy (not physical area)")
    for ax in axes:
        ax.set_xlabel("Time since cytokinesis initiation")
        ax.grid(alpha=.2)
    axes[0].legend(fontsize=8)
    fig.suptitle("Single cleavage · identical oblique axis · regulatory feedback disabled")
    fig.savefig(args.output / "comparison.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
