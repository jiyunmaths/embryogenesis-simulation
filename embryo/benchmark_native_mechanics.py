"""Idle-machine backend comparison; no scientific studies are resumed.

Run with OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1. Compilation, restore,
checkpoint IO and comparison diagnostics are outside the timed step loops.
"""
import argparse
import json
import os
from pathlib import Path
import platform
import time
import numpy as np
from .attribute_development import AttributeSimulation
from .fast_mechanics import FastAttributeSimulation
from .native_mechanics import NativeSimulation, kernel
from .feedback_long import digest
from .resolution import write_json


def run(checkpoint, output, steps=128, thread_steps=24):
    checkpoint = Path(checkpoint).resolve()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    kernel()
    simulations = {
        'reference': AttributeSimulation.restore(checkpoint),
        'accepted_fast': FastAttributeSimulation.restore(checkpoint),
        'native_4': NativeSimulation.restore(checkpoint),
    }
    simulations['native_4'].native_threads = 4
    elapsed = {key: [] for key in simulations}
    errors = {key: dict(phi_abs=0., chemical_log_abs=0., polarity_abs=0., volume_relative=0.)
              for key in simulations if key != 'reference'}
    for sim in simulations.values():
        sim.step()
    for step in range(steps):
        order = list(simulations)
        order = order[step % len(order):] + order[:step % len(order)]
        for name in order:
            start = time.perf_counter()
            simulations[name].step()
            elapsed[name].append(time.perf_counter()-start)
        ref = simulations['reference']
        for name, maxima in errors.items():
            sim = simulations[name]
            current = dict(phi_abs=float(abs(ref.phi-sim.phi).max()),
                chemical_log_abs=float(abs(np.log(np.array([ref.activator, ref.inhibitor])/
                    np.array([sim.activator, sim.inhibitor]))).max()),
                polarity_abs=float(abs(ref.polarity-sim.polarity).max()),
                volume_relative=float(abs(ref.volumes()/sim.volumes()-1).max()))
            for key, value in current.items():
                maxima[key] = max(maxima[key], value)
        if (step+1) % 32 == 0:
            print(json.dumps(dict(completed_steps=step+1, maximum_errors=errors)), flush=True)
    timings = {name: dict(mean=float(np.mean(rows)), median=float(np.median(rows)))
               for name, rows in elapsed.items()}
    report = dict(checkpoint=str(checkpoint), checkpoint_sha256=digest(checkpoint),
        steps=steps, dt=ref.config.dt, grid=ref.config.grid, cells=len(ref.ids),
        seconds_per_step=timings, maximum_errors=errors,
        exact_agreement=all(value == 0 for row in errors.values() for value in row.values()),
        native_speedup_vs_reference=timings['reference']['mean']/timings['native_4']['mean'],
        native_speedup_vs_accepted_fast=timings['accepted_fast']['mean']/timings['native_4']['mean'],
        scope='Short mature trajectory comparison, interleaved idle-machine timings. Not a full developmental or full-horizon response validation.',
        platform=platform.platform(), thread_environment={key: os.environ.get(key) for key in
            ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS')},
        source_sha256={str(p.resolve()): digest(p) for p in [Path(__file__), *[
            Path(__file__).with_name(f) for f in ('native_mechanics.py', 'native_mechanics.cpp',
                'fast_mechanics.py', 'fast_mechanics.c', 'attribute_development.py', 'model.py', 'polarity.py')]]})
    write_json(output/'comparison.json', report)
    print(json.dumps(report, indent=2), flush=True)
    del simulations
    rows = []
    # Alternate order over two rounds to reduce warmup/order effects.
    for trial, order in enumerate(((1, 2, 4, 6), (6, 4, 2, 1))):
        for threads in order:
            sim = NativeSimulation.restore(checkpoint)
            sim.native_threads = threads
            sim.step()
            start = time.perf_counter()
            for _ in range(thread_steps):
                sim.step()
            seconds = time.perf_counter()-start
            rows.append(dict(trial=trial, threads=threads, steps=thread_steps,
                             seconds_per_step=seconds/thread_steps))
            write_json(output/'thread_scaling.json', dict(trials=rows))
            print(json.dumps(rows[-1]), flush=True)
    if not report['exact_agreement']:
        raise RuntimeError('Expected exact agreement failed; inspect report before adoption')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checkpoint', type=Path, default=Path('outputs/cell-exchange-response-moving/unexchanged/unexchanged_control/latest_state.npz'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--steps', type=int, default=128)
    parser.add_argument('--thread-steps', type=int, default=24)
    args = parser.parse_args()
    if min(args.steps, args.thread_steps) < 1:
        parser.error('step counts must be positive')
    run(args.checkpoint, args.output, args.steps, args.thread_steps)
