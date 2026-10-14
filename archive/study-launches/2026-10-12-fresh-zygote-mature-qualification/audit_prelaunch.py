"""Independent mature-clock audit: NumPy and hashes, no embryo imports."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json

import numpy as np

ROOT=Path('output/studies/fresh-zygote-mature-qualification')
PARENT=Path('output/studies/fresh-zygote-control-completion')
ORIGINAL=Path('output/studies/fresh-zygote-pilot')
HERE=Path(__file__).resolve().parent
PARENT_SHA='44af9b2bfdd77176891e015e5024db2aeca6acca413f7d8dd63e5b8a86ffd926'
HISTORIES=(1101,1102,1103)
COUPLINGS=dict(off=(False,0.,0.,0.),tension_adhesion=(True,.25,.35,0.),polarity=(True,0.,0.,.35),full=(True,.25,.35,.35))
LEVELS=(('coarse',.00375),('fine',.001875))


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def read(path):return json.loads(Path(path).read_text())


def main():
    output=HERE/'prelaunch-audit.json';assert not output.exists()
    p=read(ROOT/'protocol.json');ph=sha(ROOT/'protocol.json');parent=read(PARENT/'protocol.json')
    assert sha(PARENT/'protocol.json')==p['parent_protocol_sha256']==PARENT_SHA
    assert read(ROOT/'status.json')==dict(state='prepared',protocol_sha256=ph)
    for name,expected in {**p['source_sha256'],**p['input_sha256']}.items():assert sha(name)==expected,name
    for mapping in ('source_sha256','input_sha256'):
        for name,expected in parent[mapping].items():assert p[mapping][name]==expected,name
    for key in ('physical_limits','time_limits','dense_time_spectrum_limit','prefix_steps','prefix_limits',
                'method','gpu_uuid','runtime_binary_sha256'):assert p[key]==parent[key],key
    assert p['time_limits']['field_max']==p['time_limits']['relative_spectrum_max']==p['dense_time_spectrum_limit']==.01
    assert p['histories']==list(HISTORIES) and p['branches']==list(COUPLINGS)
    assert p['new_scientific_histories']==0 and p['new_numerical_paths']==24 and len(p['jobs'])==24
    assert p['start_time']==18. and p['duration']==6. and p['interval']==.15 and p['diagnostic_interval']==.015
    assert p['snapshots']==[0,3,6] and p['restart_elapsed']==3. and p['restart_steps']==8
    assert p['homologous_limits']==dict(chemical_log_max=.01,polarity_abs_max=.01,relative_volume_max=.005,
        relative_transport_max=.01,unmatched_observations_max=0)
    assert p['levels']==[dict(level=l,dt=d) for l,d in LEVELS]
    a=read('docs/fresh_zygote_control_completion_assessment.json');record=a['independent_review'];review=read(record['path'])
    assert a['review_completed'] and a['scoped_combined_developmental_preparations_accepted'] and not a['original_pilot_accepted']
    assert sha(record['path'])==record['sha256'] and review['scoped_combined_developmental_preparations_accepted']
    expected=[];array_count=0
    for h in HISTORIES:
        for branch,coupling in COUPLINGS.items():
            endpoint=next(e for e in a['combined_primary_endpoints'] if e['history']==h and e['branch']==branch)
            source=(PARENT if h==1103 else ORIGINAL)/endpoint['case']/'gpu-t18.npz'
            with np.load(source,allow_pickle=False) as z:arrays={k:z[k].copy() for k in z.files}
            meta=json.loads(str(arrays['metadata']));base=meta.pop('config');old_step=meta.pop('step_number')
            assert base['grid']==120 and base['extent']==2.24 and base['dt']==.000234375 and base['max_cells']==16
            assert old_step==base['steps']==76800 and meta['time']==18. and not meta['divisions']
            assert arrays['phi'].shape==(16,120,120,120) and arrays['phi'].dtype==np.float32
            assert arrays['phase_carry'].shape==arrays['phi'].shape and arrays['phase_carry'].dtype==np.float64
            assert arrays['chemical_correction'].shape==(2,16) and arrays['chemical_correction'].dtype==np.float64
            assert len(arrays['ids'])==16 and len(json.loads(str(arrays['event_audits'])))==15
            assert tuple(base[k] for k in ('feedback','fate_tension','fate_adhesion','polarity_tension'))==coupling
            for level,dt in LEVELS:
                key=f'h{h}_{branch}_{level}';expected.append(key);job=next(j for j in p['jobs'] if j['key']==key)
                folder=ROOT/key;cfg=job['config']
                assert not job['reused'] and job['seed']==h and job['branch']==branch and job['level']==level
                assert job['source_case']==endpoint['case'] and Path(job['source_checkpoint']).resolve()==source.resolve()
                assert Path(job['source']).resolve()==(folder/'source.npz').resolve() and not (folder/'gpu-result.json').exists()
                assert cfg==dict(base,dt=dt,steps=round(24/dt),save_every=round(.15/dt))
                with np.load(job['source'],allow_pickle=False) as z:
                    assert set(z.files)==set(arrays)
                    current=json.loads(str(z['metadata']));assert current.pop('config')==cfg
                    assert current.pop('step_number')==round(18/dt) and current==meta
                    for name in set(arrays)-{'metadata'}:
                        assert z[name].dtype==arrays[name].dtype and z[name].shape==arrays[name].shape
                        assert z[name].tobytes()==arrays[name].tobytes(),(key,name)
                        array_count+=1
    assert [j['key'] for j in p['jobs']]==expected
    assert p['time_pairs']==[dict(seed=h,branch=b,cases=[f'h{h}_{b}_{l}' for l,_ in LEVELS]) for h in HISTORIES for b in COUPLINGS]
    assert read(ORIGINAL/'status.json')['state']=='failed' and not read(ORIGINAL/'time-h1103-off.json')['passed']
    runtime=read(PARENT/'runtime.json');assert sha(runtime['binary'])==p['runtime_binary_sha256']==runtime['binary_sha256']
    tests=read(HERE/'test-verification.json');assert tests['exit_code']==0 and sha(HERE/'tests.log')==tests['log_sha256']
    migration=Path('archive/output-relocation/2026-10-10-output/studies-manifest.json')
    assert sha(migration)=='8fb2545d8bae35a1a673ec3a43d3cc191148f03cac58892f29bd368182a8aff6'
    result=dict(passed=True,utc=datetime.now(timezone.utc).isoformat(),protocol_sha256=ph,
        source_count=len(p['source_sha256']),input_count=len(p['input_sha256']),
        losslessly_checked_start_arrays=array_count,all24_complete_starts_preserved=True,
        model_couplings_and_compensation_unchanged=True,original_failure_retained=True,
        histories=list(HISTORIES),new_scientific_histories=0,new_numerical_paths=24,long_science_launched=False,
        test_log_sha256=sha(HERE/'tests.log'),test_verification_sha256=sha(HERE/'test-verification.json'),
        migration_record_sha256=sha(migration),runtime_binary_sha256=p['runtime_binary_sha256'],audit_script_sha256=sha(__file__))
    output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
