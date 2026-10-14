"""Independent frozen-design/array audit with no embryo imports or GPU steps."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json

import numpy as np

ROOT = Path('output/studies/fresh-zygote-control-completion')
PARENT = Path('output/studies/fresh-zygote-time-refinement')
ORIGINAL = Path('output/studies/fresh-zygote-pilot')
HERE = Path(__file__).resolve().parent
PARENT_SHA = '54632a1f284affcbcec93cd0f1036d2d9c3788bd93ef5480df5fa1f14336c6c3'
REFERENCES = ('h1103_off_development_half','h1103_off_development_quarter')
COUPLINGS = dict(tension_adhesion=(True,.25,.35,0.), polarity=(True,0.,0.,.35), full=(True,.25,.35,.35))
LEVELS = (('development_half',.000234375),('development_quarter',.0001171875))


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024), b''):h.update(block)
    return h.hexdigest()


def read(path):return json.loads(Path(path).read_text())


def main():
    output = HERE/'prelaunch-audit.json'
    assert not output.exists()
    p=read(ROOT/'protocol.json'); ph=sha(ROOT/'protocol.json'); parent=read(PARENT/'protocol.json')
    assert sha(PARENT/'protocol.json')==p['parent_protocol_sha256']==PARENT_SHA
    assert read(ROOT/'status.json')==dict(state='prepared',protocol_sha256=ph)
    for name,expected in {**p['source_sha256'],**p['input_sha256']}.items():assert sha(name)==expected,name
    for mapping in ('source_sha256','input_sha256'):
        for name,expected in parent[mapping].items():assert p[mapping][name]==expected,name
    inherited=('duration','interval','checkpoint_interval','diagnostic_interval','snapshots',
        'time_limits','physical_limits','restart_steps','restart_stages','minimum_mature_window',
        'method','gpu_uuid','dense_time_spectrum_limit','instrumentation_gate',
        'instrumentation_protocol_sha256','runtime_binary_sha256','prefix_steps','prefix_limits')
    for key in inherited:assert p[key]==parent[key],key
    assert p['time_limits']['field_max']==p['time_limits']['relative_spectrum_max']==p['dense_time_spectrum_limit']==.01
    assert p['histories']==[1103] and list(COUPLINGS)==p['branches']
    assert p['new_scientific_histories']==0 and p['new_numerical_paths']==6 and p['reused_numerical_paths']==2
    assert len(p['jobs'])==8 and [j['key'] for j in p['jobs'][:2]]==list(REFERENCES)
    base=p['jobs'][0]['config']; evidence_count=0
    for job in p['jobs'][:2]:
        old=p['reuse_provenance'][job['key']]; folder=ROOT/job['key']
        assert job['reused'] and folder.is_symlink() and folder.resolve()==Path(old['folder']).resolve()
        report=read(folder/'gpu-result.json')
        assert report['passed'] and report['job']==old['job']
        assert report['protocol_sha256']==old['protocol_sha256'] and sha(folder/'gpu-result.json')==old['report_sha256']
        assert report['job']['config']==job['config']
        for name,expected in report['evidence_sha256'].items():assert sha(name)==expected,name
        evidence_count+=len(report['evidence_sha256'])
    expected=[]
    for branch,coupling in COUPLINGS.items():
        for level,dt in LEVELS:
            key=f'h1103_{branch}_{level}'; expected.append(key)
            job=next(j for j in p['jobs'] if j['key']==key); cfg=job['config']; folder=ROOT/key
            assert not job['reused'] and not folder.is_symlink() and not (folder/'gpu-result.json').exists()
            assert job['seed']==cfg['seed']==1103 and cfg['grid']==120 and cfg['graph_contact_cutoff']==.02
            assert job['branch']==branch and job['level']==level and cfg['dt']==dt
            assert cfg['steps']==round(18/dt) and cfg['save_every']==round(.15/dt)
            assert tuple(cfg[k] for k in ('feedback','fate_tension','fate_adhesion','polarity_tension'))==coupling
            assert {k for k,v in cfg.items() if v!=base[k]} <= {
                'dt','steps','save_every','feedback','fate_tension','fate_adhesion','polarity_tension'}
            assert Path(job['source']).resolve()==(folder/'source.npz').resolve()
    assert [j['key'] for j in p['jobs'][2:]]==expected
    assert p['time_pairs']==[dict(seed=1103,branch=b,cases=[f'h1103_{b}_{l}' for l,_ in LEVELS]) for b in COUPLINGS]
    with np.load(p['jobs'][0]['source'],allow_pickle=False) as z:
        arrays={name:z[name].copy() for name in z.files}; meta=json.loads(str(z['metadata']));meta.pop('config')
    assert meta['time']==0 and meta['step_number']==0 and not meta['divisions']
    for job in p['jobs']:
        with np.load(job['source'],allow_pickle=False) as z:
            assert set(z.files)==set(arrays)
            current=json.loads(str(z['metadata']));assert current.pop('config')==job['config'] and current==meta
            for name in set(arrays)-{'metadata','attribute_mode'}:np.testing.assert_array_equal(z[name],arrays[name],err_msg=(job['key'],name))
            assert str(z['attribute_mode'])==('no_feedback' if job['branch']=='off' else 'direct')
    rngs=[np.random.default_rng(s) for s in np.random.SeedSequence(1103).spawn(3)]
    due=base['division_interval']*rngs[0].uniform(1-base['cycle_jitter'],1+base['cycle_jitter'])
    np.testing.assert_array_equal(arrays['due'],[due])
    for name,rng in zip(('rng','fate_rng','signal_rng'),rngs):assert meta[name]==rng.bit_generator.state
    accepted=read(p['accepted_assessment']); record=accepted['independent_review']; review=read(record['path'])
    assert accepted['review_completed'] and accepted['scoped_refinement_accepted'] and not accepted['original_cohort_accepted']
    assert sha(record['path'])==record['sha256'] and review['review_completed'] and review['scoped_refinement_accepted']
    assert read(ORIGINAL/'status.json')['state']=='failed' and not read(ORIGINAL/'time-h1103-off.json')['passed']
    runtime=read(PARENT/'runtime.json'); assert sha(runtime['binary'])==p['runtime_binary_sha256']==runtime['binary_sha256']
    tests=read(HERE/'test-verification.json'); assert tests['exit_code']==0 and sha(HERE/'tests.log')==tests['log_sha256']
    migration=Path('archive/output-relocation/2026-10-10-output/studies-manifest.json')
    assert sha(migration)=='8fb2545d8bae35a1a673ec3a43d3cc191148f03cac58892f29bd368182a8aff6'
    result=dict(passed=True,utc=datetime.now(timezone.utc).isoformat(),protocol_sha256=ph,
        source_count=len(p['source_sha256']),input_count=len(p['input_sha256']),
        immutable_reference_evidence_count=evidence_count,complete_initial_arrays_match=True,
        initial_rng_independently_reproduced=True,unchanged_model_and_limits=True,
        original_failure_retained=True,new_scientific_histories=0,new_numerical_paths=6,reused_numerical_paths=2,
        test_log_sha256=sha(HERE/'tests.log'),test_verification_sha256=sha(HERE/'test-verification.json'),
        migration_record_sha256=sha(migration),runtime_binary_sha256=p['runtime_binary_sha256'],audit_script_sha256=sha(__file__))
    output.write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps(result,indent=2))


if __name__=='__main__':main()
