"""Matched coupled CPU/GPU timing after the full scientific agreement gate."""
import argparse
import json
from pathlib import Path
import time

import numpy as np
import torch

from .cell_response_moving import observe
from .feedback_long import digest
from .gpu_backend import GpuSimulation
from .gpu_response_runner import require_validation
from .native_mechanics import NativeSimulation
from .resolution import write_json
from .validate_gpu_backend import discrepancies


def cpu_audit(sim):
    volumes = sim.volumes()
    if (np.max(abs(volumes/sim.target-1)) >= .05 or sim.clipped_fraction or
            np.min((3*volumes/(4*np.pi))**(1/3)/sim.dx) < 4 or
            not all(np.isfinite(x).all() for x in (volumes, sim.activator, sim.inhibitor, sim.polarity)) or
            np.any(sim.activator <= 0) or np.any(sim.inhibitor <= 0)):
        raise RuntimeError('CPU benchmark quality failure')


def run(output, validation, steps=120, threads=4):
    if steps < 40 or steps % 40 or not 1 <= threads <= 6:
        raise ValueError('Use a positive multiple of 40 steps and 1–6 CPU threads')
    output = Path(output)
    if output.exists():
        raise FileExistsError(output)
    require_validation(validation)
    p = json.loads((Path(validation)/'protocol.json').read_text())
    checkpoint = Path(p['baseline'])/'unexchanged/source.npz'
    torch.set_num_threads(1)
    cpu = NativeSimulation.restore(checkpoint)
    cpu.native_threads = threads
    gpu = GpuSimulation.restore(checkpoint)
    # Equal warmup advances; timed trajectories cover the same physical interval.
    for _ in range(5):
        cpu.step(); gpu.step()
    torch.cuda.synchronize()
    output.mkdir(parents=True)
    records, endpoints = {}, {}
    for backend, sim in (('native', cpu), ('gpu', gpu)):
        started = time.perf_counter()
        samples = []
        for step in range(steps):
            sim.step()
            if backend == 'native': cpu_audit(sim)
            else: sim.audit()
            if (step+1) % 40 == 0:
                row = observe(sim, (step+1)*sim.config.dt) if backend == 'native' else sim.observe((step+1)*sim.config.dt)
                if row['boundary_occupancy'] >= .01:
                    raise RuntimeError('Boundary screen failed')
                samples.append(row)
        torch.cuda.synchronize()
        wall = time.perf_counter()-started
        started = time.perf_counter()
        sim.checkpoint(output/f'{backend}.npz')
        checkpoint_wall = time.perf_counter()-started
        endpoints[backend] = samples[-1]
        records[backend] = dict(steps=steps, wall_seconds=wall, seconds_per_step=wall/steps,
                               checkpoint_seconds=checkpoint_wall)
        write_json(output/f'{backend}-history.json', samples)
    errors = discrepancies(endpoints['gpu'], endpoints['native'])
    phi_error = float(abs(gpu.phi.cpu().numpy()-cpu.phi).max())
    passed = all(v <= p['criteria'][k] for k, v in errors.items()) and phi_error <= p['criteria']['final_phi_abs_max']
    report = dict(passed=bool(passed), timings=records,
        coupled_step_speedup=records['native']['seconds_per_step']/records['gpu']['seconds_per_step'],
        endpoint_errors=errors, phi_abs_max=phi_error, native_threads=threads,
        device=torch.cuda.get_device_name(0), dt=cpu.config.dt, steps=steps,
        source_sha256={str(Path(__file__).resolve()): digest(__file__)},
        validation_protocol_sha256=digest(Path(validation)/'protocol.json'),
        scope='Sequential matched short coupled continuations after equal five-step warmup. Whole step includes every-step quality audits and 0.15-spaced scientific diagnostics. Standard checkpoint time is reported separately; no startup/compile or full-assay I/O speedup claim. One GPU versus one native worker with the reported OpenMP thread count, BLAS one thread. Not a multi-worker or multi-GPU benchmark.')
    write_json(output/'comparison.json', report)
    print(json.dumps(report, indent=2), flush=True)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('outputs/resident-gpu-coupled-benchmark'))
    parser.add_argument('--validation', type=Path, default=Path('outputs/gpu-backend-validation'))
    parser.add_argument('--steps', type=int, default=120)
    parser.add_argument('--threads', type=int, default=4)
    args = parser.parse_args()
    run(args.output, args.validation, args.steps, args.threads)
