"""Independent initial-state and evidence audit; no embryo imports."""
from pathlib import Path
import hashlib
import json
import numpy as np

ROOT=Path('output/studies/development-time-diagnosis')
HERE=Path(__file__).parent


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def main():
    p=json.loads((ROOT/'protocol.json').read_text());protocol_hash=sha(ROOT/'protocol.json')
    status=json.loads((ROOT/'status.json').read_text());assert status['state']=='prepared' and status['protocol_sha256']==protocol_hash
    for file,expected in {**p['source_sha256'],**p['input_sha256']}.items():assert sha(file)==expected,file
    parent=Path(p['parent_study']);assert sha(parent/'protocol.json')==p['parent_protocol_sha256']
    baseline=json.loads((parent/'protocol.json').read_text())
    cfg=next(j['config'] for j in baseline['jobs'] if j['key']=='off_coarse')
    assert not cfg['feedback'] and cfg['fate_tension']==cfg['fate_adhesion']==cfg['polarity_tension']==0
    dts=[.00375,.001875,.0009375,.00046875]
    with np.load(parent/'off_coarse/source.npz',allow_pickle=False) as old:
        oldmeta=json.loads(str(old['metadata']));rows=[]
        for job,dt in zip(p['jobs'],dts):
            changed={k:v for k,v in job['config'].items() if v!=cfg[k]}
            assert set(changed).issubset({'dt','steps','save_every'})
            assert job['config']['dt']==dt and job['config']['steps']==round(18/dt)
            assert job['config']['save_every']==round(.15/dt)
            with np.load(job['source'],allow_pickle=False) as new:
                assert set(old.files)==set(new.files)
                for key in old.files:
                    if key!='metadata':np.testing.assert_array_equal(old[key],new[key])
                meta=json.loads(str(new['metadata']));assert meta['config']==job['config']
                assert {k:v for k,v in meta.items() if k!='config'}=={k:v for k,v in oldmeta.items() if k!='config'}
                assert np.all(new['phase_carry']==0) and np.all(new['chemical_correction']==0)
                assert np.all(new['polarity']==0) and np.all(new['fate']==0)
                rows.append(dict(case=job['key'],dt=dt,initial_arrays_and_rng_match=True))
    gate=json.loads((ROOT/'instrumentation-gate.json').read_text())
    assert gate['passed'] and gate['protocol_sha256']==protocol_hash
    assert [(r['stage'],r['steps'],r['bitwise_equal']) for r in gate['rows']]==[(2,8,True),(4,8,True),(8,8,True)]
    for file,expected in gate['evidence_sha256'].items():assert sha(file)==expected,file
    preserved=json.loads(Path(p['preserved_manifest']).read_text())
    for file,expected in preserved['preserved_sha256'].items():assert sha(file)==expected,file
    assert '38 passed' in (HERE/'tests.log').read_text()
    result=dict(passed=True,protocol_sha256=protocol_hash,source_count=len(p['source_sha256']),input_count=len(p['input_sha256']),
        preserved_count=len(preserved['preserved_sha256']),initial_states=rows,instrumentation_gate_sha256=sha(ROOT/'instrumentation-gate.json'),
        test_log_sha256=sha(HERE/'tests.log'),script_sha256=sha(__file__),new_scientific_histories=0)
    (HERE/'prelaunch-audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
