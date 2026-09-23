"""Small deterministic four-cell checks; not a full convergence certification."""

import argparse
from dataclasses import asdict
import json
from pathlib import Path

from .model import Config, Simulation


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output already exists; choose a new path")
    cases = [(32, .015), (40, .015), (48, .015), (40, .0075)]
    rows = []
    for grid, dt in cases:
        c = Config(grid=grid, dt=dt, steps=round(6 / dt), max_cells=4,
                   division_interval=1.5, cycle_jitter=0, differentiation=False,
                   feedback=False, signaling=False, polarity_enabled=False)
        sim = Simulation(c)
        for _ in range(c.steps):
            sim.step()
        row = {"config": asdict(c), "final": sim.metrics()}
        rows.append(row)
        print(f'grid={grid} dt={dt}: axis={row["final"]["axis_ratio"]:.5f}, '
              f'volume error={row["final"]["relative_volume_error"]:+.3%}', flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(rows, indent=2) + "\n")


if __name__ == "__main__":
    main()
