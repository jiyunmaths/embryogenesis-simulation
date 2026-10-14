"""Independently check the diagnostic design and saved starts; no embryo imports."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json

import numpy as np

ROOT = Path('output/studies/fresh-zygote-time-refinement')
PARENT = Path('output/studies/fresh-zygote-pilot')
HERE = Path(__file__).resolve().parent
PARENT_SHA = '728f584676ea4eeeee530656015aec61931e9296ff095eff8b7c94f9832a13ee'


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def read(path):return json.loads(Path(path).read_text())


def main():
    target=HERE/'prelaunch-audit.json'
    if target.exists():raise FileExistsError('Never overwrite a prelaunch audit')
    p=read(ROOT/'protocol.json');ph=sha(ROOT/'protocol.json')
    assert read(ROOT/'status.json')==dict(state='prepared',protocol_sha256=ph)
    assert sha(PARENT/'protocol.json')==p['parent_protocol_sha256']==PARENT_SHA
    parent=read(PARENT/'protocol.json')
    for file,h in {**p['source_sha256'],**p['input_sha256']}.items():assert sha(file)==h,file
    for mapping in ('source_sha256','input_sha256'):
        for file,h in parent[mapping].items():assert p[mapping][file]==h,file
    inherited=('duration','interval','checkpoint_interval','diagnostic_interval','snapshots',
        'time_limits','physical_limits','restart_steps','restart_stages','minimum_mature_window',
        'method','gpu_uuid','dense_time_spectrum_limit','instrumentation_gate',
        'instrumentation_protocol_sha256','runtime_binary_sha256','prefix_steps','prefix_limits')
    for key in inherited:assert p[key]==parent[key],key
    assert p['time_limits']['field_max']==p['time_limits']['relative_spectrum_max']==.01
    assert p['dense_time_spectrum_limit']==.01
    assert p['new_scientific_histories']==0 and p['new_numerical_paths']==1
    assert p['reused_numerical_paths']==1 and p['histories']==[1103] and p['branches']==['off']
    a,b=p['jobs'];assert a['key']=='h1103_off_development_half' and a['reused']
    assert b['key']=='h1103_off_development_quarter' and not b['reused']
    assert a['seed']==b['seed']==1103 and a['branch']==b['branch']=='off'
    original=next(j for j in parent['jobs'] if j['key']==a['key'])
    assert a==dict(original,reused=True)
    cfg=b['config'];base=a['config']
    assert {k for k,v in cfg.items() if v!=base[k]}=={'dt','steps','save_every'}
    assert base['dt']==.000234375 and cfg['dt']==.0001171875
    assert cfg['steps']==153600 and cfg['save_every']==1280
    assert cfg['grid']==120 and cfg['graph_contact_cutoff']==.02
    assert cfg['feedback'] is False and cfg['fate_tension']==cfg['fate_adhesion']==cfg['polarity_tension']==0
    assert p['time_pair']==dict(seed=1103,branch='off',cases=[a['key'],b['key']])
    ref=ROOT/a['key'];new=ROOT/b['key']
    assert ref.is_symlink() and ref.resolve()==(PARENT/a['key']).resolve()
    assert not new.is_symlink() and not (new/'gpu-result.json').exists()
    report=read(ref/'gpu-result.json');old=p['reuse_provenance'][a['key']]
    assert report['passed'] and report['job']==old['job']==original
    assert report['protocol_sha256']==old['protocol_sha256']==PARENT_SHA
    assert sha(ref/'gpu-result.json')==old['report_sha256']
    for file,h in report['evidence_sha256'].items():assert sha(file)==h,file
    failed=read(PARENT/'time-h1103-off.json');status=read(PARENT/'status.json')
    assert not failed['passed'] and status['state']=='failed' and status['completed']==18
    assert failed['dense']['peak']['relative_spectrum']>.01
    assert sha(PARENT/'status.json')==p['original_failed_status_sha256']
    assert sha(PARENT/'time-h1103-off.json')==p['original_failed_pair_sha256']
    assessment=read(p['parent_assessment'])
    assert assessment['review_completed'] and not assessment['qualification_accepted']
    review_record=assessment['independent_saved_evidence_review']
    assert sha(review_record['path'])==review_record['sha256']
    review=read(review_record['path']);assert review['review_completed'] and review['numerical_stop_confirmed']
    # Compare complete arrays and non-config metadata directly, not production signatures.
    with np.load(a['source'],allow_pickle=False) as x,np.load(b['source'],allow_pickle=False) as y:
        assert set(x.files)==set(y.files)
        mx,my=[json.loads(str(z['metadata'])) for z in (x,y)]
        assert mx.pop('config')==base and my.pop('config')==cfg and mx==my
        assert mx['time']==0 and mx['step_number']==0 and not mx['divisions']
        assert mx['next_id']==1 and not json.loads(str(y['event_audits']))
        for name in set(x.files)-{'metadata'}:np.testing.assert_array_equal(x[name],y[name],err_msg=name)
        assert y['phi'].shape==(1,120,120,120) and y['phi'].dtype==np.float32
        assert y['phase_carry'].dtype==y['chemical_correction'].dtype==np.float64
        np.testing.assert_array_equal(y['activator'],[1.03]);np.testing.assert_array_equal(y['inhibitor'],[.98])
        np.testing.assert_array_equal(y['ids'],[0]);np.testing.assert_array_equal(y['parents'],[-1])
        for name in ('phase_carry','chemical_correction','polarity','fate','rounding'):assert not y[name].any()
        assert str(y['attribute_mode'])=='no_feedback' and str(y['chemical_method'])==p['method']
        # Independently regenerate the original history's streams and cycle schedule.
        rngs=[np.random.default_rng(s) for s in np.random.SeedSequence(1103).spawn(3)]
        due=cfg['division_interval']*rngs[0].uniform(1-cfg['cycle_jitter'],1+cfg['cycle_jitter'])
        np.testing.assert_array_equal(y['due'],[due])
        for name,rng in zip(('rng','fate_rng','signal_rng'),rngs):assert my[name]==rng.bit_generator.state
    gate=read(p['instrumentation_gate'])
    assert gate['passed'] and gate['protocol_sha256']==p['instrumentation_protocol_sha256']
    for file,h in gate['evidence_sha256'].items():assert sha(file)==h,file
    runtime=read(PARENT/'runtime.json')
    assert sha(runtime['binary'])==runtime['binary_sha256']==p['runtime_binary_sha256']
    tests=read(HERE/'test-verification.json')
    assert tests['exit_code']==0 and sha(HERE/'tests.log')==tests['log_sha256']
    assert 'passed' in (HERE/'tests.log').read_text() and 'failed' not in (HERE/'tests.log').read_text()
    migration=Path('archive/output-relocation/2026-10-10-output/studies-manifest.json')
    assert sha(migration)=='8fb2545d8bae35a1a673ec3a43d3cc191148f03cac58892f29bd368182a8aff6'
    result=dict(passed=True,utc=datetime.now(timezone.utc).isoformat(),protocol_sha256=ph,
        source_count=len(p['source_sha256']),input_count=len(p['input_sha256']),
        immutable_reference_evidence_count=len(report['evidence_sha256']),
        original_failure_retained=True,unchanged_model_and_limits=True,complete_initial_arrays_match=True,
        initial_rng_independently_reproduced=True,new_scientific_histories=0,new_numerical_paths=1,
        reused_numerical_paths=1,test_log_sha256=sha(HERE/'tests.log'),
        test_verification_sha256=sha(HERE/'test-verification.json'),
        migration_record_sha256=sha(migration),runtime_binary_sha256=p['runtime_binary_sha256'],
        audit_script_sha256=sha(__file__),post_cleavage_not_launched=True)
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
