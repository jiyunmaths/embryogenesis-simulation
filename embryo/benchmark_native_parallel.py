"""Compare worker/thread allocations without resuming biological studies."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
from pathlib import Path
import time
import numpy as np
from .native_mechanics import NativeSimulation, kernel
from .resolution import write_json
from .feedback_long import digest

_BARRIER = None


def initialize(barrier):
    global _BARRIER
    _BARRIER = barrier


def trial(args):
    checkpoint, threads, steps = args
    sim = NativeSimulation.restore(checkpoint)
    sim.native_threads = threads
    kernel()
    sim.step()
    _BARRIER.wait(timeout=120)
    start = time.perf_counter()
    for _ in range(steps):
        sim.step()
    end = time.perf_counter()
    return dict(start=start, end=end, seconds_per_step=(end-start)/steps,
                chemical_state=np.array([sim.activator, sim.inhibitor]).tolist())


def run(checkpoint, output, steps):
    context = mp.get_context('spawn')  # Never fork an initialized OpenMP runtime.
    output = Path(output)
    if output.exists():
        raise FileExistsError(output)
    rows = []
    for workers, threads in ((1, 4), (2, 2), (3, 2), (4, 1), (2, 1)):
        barrier = context.Barrier(workers)
        with ProcessPoolExecutor(max_workers=workers, mp_context=context,
                                 initializer=initialize, initargs=(barrier,)) as pool:
            values = list(pool.map(trial, [(str(checkpoint), threads, steps)]*workers))
        seconds = max(v['end'] for v in values)-min(v['start'] for v in values)
        errors = max(float(abs(np.asarray(v['chemical_state'])-values[0]['chemical_state']).max()) for v in values)
        row = dict(workers=workers, threads_per_worker=threads, steps_per_worker=steps,
                   aggregate_steps_per_second=steps*workers/seconds,
                   wall_seconds=seconds, individual_seconds_per_step=[v['seconds_per_step'] for v in values],
                   between_worker_chemical_error=errors)
        rows.append(row)
        write_json(output, dict(checkpoint=str(Path(checkpoint).resolve()),
            checkpoint_sha256=digest(checkpoint), trials=rows,
            scope='Idle machine, synchronized worker start, startup/restore/warmup excluded. Short throughput screen; identical checkpoint per worker.'))
        print(row, flush=True)
    return rows


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--checkpoint', type=Path, default=Path('outputs/cell-exchange-response-moving/unexchanged/unexchanged_control/latest_state.npz'))
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--steps', type=int, default=32)
    a = p.parse_args()
    if a.steps < 1:
        p.error('steps must be positive')
    run(a.checkpoint, a.output, a.steps)
