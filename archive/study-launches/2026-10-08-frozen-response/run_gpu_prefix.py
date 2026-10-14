"""Qualify every actual frozen path on CUDA before starting full production."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import sys
import time

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from embryo import phase_carry_frozen_response as study
from embryo.parameter_robustness import verify
from embryo.transport import integrate_gm


OUTPUT = ROOT / 'output/studies/phase-carry-frozen-response'
REPORT = Path(__file__).with_name('gpu-prefix.json')


def main():
    if REPORT.exists():
        raise FileExistsError(REPORT)
    p = study.read(OUTPUT / 'protocol.json')
    verify(p)
    device = study.runtime_check(p)
    if device.type != 'cuda':
        raise ValueError('This qualification requires the intended CUDA device')
    torch.set_num_threads(1)
    started = datetime.now(timezone.utc).isoformat()
    initial, matrix, mapping, count = study.batch_layout(p)
    cases = []
    contexts = {}
    for trial in p['trials']:
        _, delta, volumes, values, ids = study.load_context(p, trial)
        contexts[trial['key']] = (study.fixed_transport(delta, volumes), values, ids)
    for level in p['levels']:
        dt = level['dt']
        total, every = study._steps(.6, dt), study._steps(p['interval'], dt)
        a, b, tensor = [torch.as_tensor(x.copy(), dtype=torch.float64, device=device)
                        for x in (initial[0], initial[1], matrix)]
        substeps = study.substeps(matrix, dt, p['beta'], p['da'], p['db'])
        gpu = [study.unpack(initial, mapping, count)]
        torch.cuda.synchronize(device)
        clock = time.perf_counter()
        for step in range(1, total + 1):
            a, b = study.checked_torch_step(a, b, tensor, dt, p, substeps)
            if step % every == 0:
                values = np.stack([a.detach().cpu().numpy(), b.detach().cpu().numpy()])
                gpu.append(study.unpack(values, mapping, count))
        torch.cuda.synchronize(device)
        gpu_seconds = time.perf_counter() - clock
        gpu = np.array(gpu)
        cpu = []
        clock = time.perf_counter()
        for trial in p['trials']:
            transport, values, _ = contexts[trial['key']]
            aa, bb = values.copy()
            path = [values]
            for step in range(1, total + 1):
                aa, bb = integrate_gm(aa, bb, transport, dt, p['beta'], p['da'], p['db'])
                if step % every == 0:
                    path.append(np.stack([aa, bb]))
            cpu.append(path)
        cpu_seconds = time.perf_counter() - clock
        cpu = np.array(cpu).transpose(1, 0, 2, 3)
        raw_error = float(abs(np.log(gpu / cpu)).max())
        pulse_errors = []
        for i, trial in enumerate(p['trials']):
            if trial['target'] is None:
                continue
            control = next(j for j, t in enumerate(p['trials'])
                           if t['context'] == trial['context'] and t['target'] is None)
            amplitude = abs(np.log(trial['factor']))
            error = float(abs(np.log(gpu[:, i] / gpu[:, control])
                              - np.log(cpu[:, i] / cpu[:, control])).max() / amplitude)
            pulse_errors.append(dict(trial=trial['key'], initial_log_amplitude=amplitude,
                                     normalized_response_error=error))
        normalized = max(row['normalized_response_error'] for row in pulse_errors)
        passed = (raw_error <= p['criteria']['cpu_torch_log_max']
                  and normalized <= p['criteria']['normalized_response_max'])
        cases.append(dict(level=level, horizon=.6, logical_paths=33, samples=5,
                          substeps=substeps, raw_log_error=raw_error,
                          normalized_response_error=normalized, pulse_errors=pulse_errors,
                          gpu_wall_seconds=gpu_seconds, serial_cpu_wall_seconds=cpu_seconds,
                          estimated_full_gpu_seconds=100*gpu_seconds, passed=bool(passed)))
        print(level['level'], 'GPU prefix seconds', gpu_seconds,
              'raw/normalized errors', raw_error, normalized, flush=True)
    verify(p)
    report = dict(passed=bool(all(row['passed'] for row in cases)),
                  started_utc=started, finished_utc=datetime.now(timezone.utc).isoformat(),
                  protocol_sha256=study.digest(OUTPUT / 'protocol.json'),
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  runtime=p['runtime'], device=p['device'],
                  gpu_name=torch.cuda.get_device_name(device), levels=cases,
                  estimated_full_gpu_seconds=sum(row['estimated_full_gpu_seconds'] for row in cases),
                  scope='All 33 logical paths at three timesteps through .6 units. '
                        'Unchanged float64 PyTorch stages checked on CUDA and compared '
                        'with NumPy/SciPy; no full-window scientific outcome or mechanics stepping.')
    study.write_json(REPORT, report)
    if not report['passed']:
        raise RuntimeError('Actual CUDA frozen-prefix qualification failed')
    print('Actual GPU prefixes qualify:', REPORT, flush=True)


if __name__ == '__main__':
    main()
