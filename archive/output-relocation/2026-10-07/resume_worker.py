"""Resume the unchanged scientific study after verified output relocation."""
from datetime import datetime, timezone
from zoneinfo import ZoneInfo
from pathlib import Path
import hashlib
import json
import os
import subprocess

ROOT = Path('/home/sji/Documents/Programming/Python/embryogenesis-simulation')
ARCHIVE = ROOT/'archive/output-relocation/2026-10-07'
STUDY = ROOT/'outputs/phase-carry-redistribution-response'
PYTHON = '/home/sji/miniforge3/envs/deeplearning/bin/python'
GPU_UUID = 'GPU-2991ab3d-dcbd-41f3-eac1-133be278f56e'
read = lambda path: json.loads(Path(path).read_text())


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


proof = read(ARCHIVE/'compatibility-verification.json')
assert proof['passed'] and proof['manifest_sha256'] == sha(ARCHIVE/'manifest.json')
assert read(ARCHIVE/'retirement.json')['passed']
assert not (ROOT/'outputs.relocation-original-2026-10-07').exists()
assert not (ARCHIVE/'resume.json').exists(), 'A restart is already recorded'
old = read(ARCHIVE/'interruption.json')
p = read(STUDY/'protocol.json')
assert sha(STUDY/'protocol.json') == old['protocol_sha256']
assert sha(old['checkpoint']['path']) == old['checkpoint']['sha256']
for path, digest in p['source_sha256'].items():
    assert sha(path) == digest, path
processes = subprocess.check_output(['ps', '-C', 'python', '-o', 'pid=,args='], text=True)
assert not any('-m embryo.' in line for line in processes.splitlines()), 'Another simulation worker is active'
inventory = subprocess.check_output(['nvidia-smi', '--query-gpu=index,name,uuid', '--format=csv,noheader'], text=True)
assert '1, NVIDIA GeForce GTX 1080 Ti, '+GPU_UUID in inventory
apps = subprocess.check_output(['nvidia-smi', '--query-compute-apps=gpu_uuid,pid,process_name', '--format=csv,noheader'], text=True)
assert GPU_UUID not in apps, 'Qualified GPU has an active compute process'
command = [PYTHON, str(ROOT/'scripts/with_output_path_alias.py'), '--record', str(ARCHIVE/'manifest.json'),
    '--record-sha256', sha(ARCHIVE/'manifest.json'), '-m', 'embryo.phase_carry_redistribution_response',
    'run', '--output', str(STUDY)]
environment = os.environ.copy()
environment.update(CUDA_DEVICE_ORDER='PCI_BUS_ID', CUDA_VISIBLE_DEVICES='1', PYTHONPATH='.',
    OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', PYTHONUNBUFFERED='1')
now = datetime.now(timezone.utc)
with (STUDY/'run.log').open('ab', buffering=0) as stream:
    stream.write(('\nStorage relocation resume '+now.isoformat()+'; unchanged protocol, committed elapsed=3 checkpoint.\n').encode())
    worker = subprocess.Popen(command, cwd=ROOT, env=environment, stdin=subprocess.DEVNULL,
        stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
record = dict(passed=True, pid=worker.pid, previous_pid=old['old_pid'], command=command,
    launched_utc=now.isoformat(), launched_local=now.astimezone(ZoneInfo('America/New_York')).isoformat(),
    manifest_sha256=sha(ARCHIVE/'manifest.json'), protocol_sha256=sha(STUDY/'protocol.json'),
    resumed_checkpoint=old['checkpoint'], launcher_sha256=sha(ROOT/'scripts/with_output_path_alias.py'),
    resume_script_sha256=sha(__file__), device_uuid=GPU_UUID, physical_gpu=1,
    scope='Resume only from the last committed unchanged state; completed paths are validated in place. No new protocol or scientific outcome.')
(ARCHIVE/'resume.json').write_text(json.dumps(record, indent=2)+'\n')
print(json.dumps({k:v for k,v in record.items() if k != 'resumed_checkpoint'}, indent=2), flush=True)
