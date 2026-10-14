import hashlib
import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path
import sys
import traceback

import numpy as np
import torch

ROOT = Path('/home/sji/Documents/Programming/Python/embryogenesis-simulation')
sys.path.insert(0, str(ROOT))
ARCHIVE = ROOT/'archive/output-relocation/2026-10-07'
RECORD = ARCHIVE/'manifest.json'
RECORD_SHA = 'edcf6b433ed7d3db70adb89c8ca701ec49b0d73a5b981380c61a966bf89a3c02'
PROOF = ARCHIVE/'scientific-compatibility.json'
DEST = Path('/mnt/data/Documents/Python/embryogenesis-simulation/outputs')
CURRENT = 'phase-carry-redistribution-response'
study_root = DEST/CURRENT

def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def metadata(path):
    s = Path(path).lstat()
    return dict(dev=s.st_dev, ino=s.st_ino, mode=s.st_mode, uid=s.st_uid,
        gid=s.st_gid, nlink=s.st_nlink, size=s.st_size, mtime_ns=s.st_mtime_ns)

assert not PROOF.exists(), 'Never replace an existing scientific proof'
manifest = json.loads(RECORD.read_text())
assert sha(RECORD) == RECORD_SHA
spec = importlib.util.spec_from_file_location('relocation_path_alias', ROOT/'scripts/with_output_path_alias.py')
alias = importlib.util.module_from_spec(spec)
spec.loader.exec_module(alias)
proof = dict(passed=False, manifest_sha256=RECORD_SHA, destination_root=str(study_root),
    started_utc=datetime.now(timezone.utc).isoformat(), gpu_disabled=True,
    launcher_sha256=sha(ROOT/'scripts/with_output_path_alias.py'),
    verifier_sha256=sha(__file__))

try:
    torch.set_num_threads(1)
    assert not torch.cuda.is_available(), 'GPU must remain disabled'
    p = json.loads((study_root/'protocol.json').read_text())
    before_source = {path:sha(path) for path in p['source_sha256']}
    assert len(before_source) == 53 and before_source == p['source_sha256']
    all_files = {name:entry for name,entry in manifest['files'].items() if name.startswith(CURRENT+'/')}
    before_metadata = {name:metadata(DEST/name) for name in all_files}
    assert all(before_metadata[name] == entry['destination_stat'] for name,entry in all_files.items())
    protected = [name for name in all_files if Path(name).name in (
        'source.npz', 'latest_state.npz', 'protocol.json', 'preflight.json',
        'pilot-refinement.json', 'full-refinement.json', 'pilot-quality.json')]
    print(f'Hashing {len(protected)} prepared checkpoints and reports before recomputation', flush=True)
    before_hashes = {name:sha(DEST/name) for name in protected}
    assert all(before_hashes[name] == all_files[name]['sha256'] for name in protected)
    interruption = json.loads((ARCHIVE/'interruption.json').read_text())
    cp = Path(interruption['checkpoint']['path'])
    with np.load(cp, allow_pickle=False) as data:
        cp_meta = json.loads(str(data['metadata']))
        saved = json.loads(str(data['long_experiment']))
    cp_sha = sha(cp)
    assert cp_sha == interruption['checkpoint']['sha256']
    assert cp_meta['time'] == 513.0 and saved['history'][-1]['elapsed'] == 3.0
    assert saved['job'] == interruption['checkpoint']['job']
    assert saved['protocol_hash'] == interruption['protocol_sha256'] == sha(study_root/'protocol.json')
    reports = []
    with alias.validated_alias(RECORD, RECORD_SHA):
        from embryo import phase_carry_redistribution_response as study
        for label,horizon in [('pilot',p['pilot_duration']), ('full',p['duration'])]:
            file = study_root/'pairs/seed-7'/f'{label}-refinement.json'
            assert file.is_file()
            original = json.loads(file.read_text())
            old_sha = sha(file)
            print(f'Recomputing existing seed7 {label} decision from copied destination', flush=True)
            recomputed = study.pair_report(study_root,p,7,horizon)
            assert recomputed == original and recomputed['passed'] is True
            assert list(recomputed['evidence_sha256']) == list(original['evidence_sha256'])
            assert sha(file) == old_sha
            reports.append(dict(label=label, horizon=horizon, passed=True,
                report_sha256=old_sha, exact_decision_identity=True,
                exact_evidence_path_identity=True,
                evidence_sha256=recomputed['evidence_sha256']))
    print('Rechecking checkpoint/report hashes and all current-study file metadata', flush=True)
    after_hashes = {name:sha(DEST/name) for name in protected}
    after_metadata = {name:metadata(DEST/name) for name in all_files}
    assert after_hashes == before_hashes and after_metadata == before_metadata
    assert {path:sha(path) for path in p['source_sha256']} == before_source
    assert sha(cp) == cp_sha
    with np.load(cp, allow_pickle=False) as data:
        assert json.loads(str(data['metadata'])) == cp_meta
        assert json.loads(str(data['long_experiment'])) == saved
    proof.update(passed=True, reports=reports, protocol_sha256=sha(study_root/'protocol.json'),
        source_files_verified=53, source_sha256=before_source,
        current_study_files_metadata_unchanged=len(all_files),
        protected_checkpoint_and_report_sha256_unchanged=before_hashes,
        committed_checkpoint=dict(path=str(cp), sha256=cp_sha, time=cp_meta['time'],
            step_number=cp_meta['step_number'], elapsed=saved['history'][-1]['elapsed'],
            history_rows=len(saved['history']), metadata=cp_meta,
            filesystem_metadata=metadata(cp), exact_metadata_identity=True),
        scope='CPU-only paired report recomputation under the pinned output path alias; existing output bytes and metadata unchanged. No numerical advancement or scientific source modification.')
except BaseException as error:
    proof.update(error=str(error), traceback=traceback.format_exc())
    raise
finally:
    proof['finished_utc'] = datetime.now(timezone.utc).isoformat()
    PROOF.write_text(json.dumps(proof, indent=2, allow_nan=False)+'\n')
    print('Scientific compatibility passed:', proof['passed'], flush=True)
