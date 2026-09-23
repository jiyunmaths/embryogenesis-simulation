"""Paired-seed screening of mechanics-only, one-way, and two-way feedback."""

import argparse
import csv
from dataclasses import asdict
import json
from pathlib import Path

from .model import Config, Simulation


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", nargs="+", type=int, default=[7, 11, 23])
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--grid", type=int, default=40)
    parser.add_argument("--steps", type=int, default=1000)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output directory already exists; choose a new path")
    Config(grid=args.grid, steps=args.steps).validate()
    args.output.mkdir(parents=True)
    rows, configs = [], []
    for seed in args.seeds:
        for mode, differentiation, feedback in [("mechanics", False, False), ("one_way", True, False), ("two_way", True, True)]:
            c = Config(seed=seed, grid=args.grid, steps=args.steps,
                       differentiation=differentiation, feedback=feedback,
                       signaling=differentiation, polarity_enabled=differentiation)
            sim = Simulation(c)
            history = [sim.metrics()]
            for step in range(c.steps):
                sim.step()
                if (step + 1) % c.save_every == 0 or step + 1 == c.steps:
                    history.append(sim.metrics())
            row = {"seed": seed, "mode": mode, **history[-1],
                   "max_observed_volume_error": max(abs(m["relative_volume_error"]) for m in history)}
            rows.append(row)
            configs.append({"mode": mode, **asdict(c)})
            (args.output / f"{seed}_{mode}.json").write_text(json.dumps(history, indent=2) + "\n")
            print(f'{seed} {mode}: A/B/U={row["fate_a"]}/{row["fate_b"]}/{row["uncommitted"]}, '
                  f'axis={row["axis_ratio"]:.3f}, volume error={row["relative_volume_error"]:+.2%}', flush=True)
    with (args.output / "comparison.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
    (args.output / "configs.json").write_text(json.dumps(configs, indent=2) + "\n")


if __name__ == "__main__":
    main()
