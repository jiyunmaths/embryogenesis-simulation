"""Check the fresh-history design and saved starts without embryo imports."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json

import numpy as np
from scipy.integrate import quad
from scipy.special import expit

HERE=Path(__file__).parent
ROOT=Path('output/studies/fresh-zygote-pilot')


def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):digest.update(block)
    return digest.hexdigest()


def read(path):return json.loads(Path(path).read_text())


def main():
    p=read(ROOT/'protocol.json');ph=sha(ROOT/'protocol.json')
    assert read(ROOT/'status.json')==dict(state='prepared',protocol_sha256=ph)
    for file,h in {**p['source_sha256'],**p['input_sha256']}.items():assert sha(file)==h,file
    preserved=read(p['preserved_manifest'])['preserved_sha256']
    for file,h in preserved.items():assert sha(file)==h,file
    parent_root=Path('output/studies/development-control-qualification')
    parent=read(parent_root/'protocol.json');assert sha(parent_root/'protocol.json')==p['parent_protocol_sha256']
    for file,h in parent['source_sha256'].items():assert p['source_sha256'][file]==h,file
    accepted=read(p['accepted_parent_assessment'])
    assert accepted['review_completed'] and accepted['qualification_accepted']
    assert accepted['protocol_sha256']==p['parent_protocol_sha256']
    assert sha(p['accepted_parent_review'])==accepted['independent_review']['sha256']
    for key in ('duration','interval','checkpoint_interval','diagnostic_interval','snapshots',
        'time_limits','physical_limits','restart_steps','restart_stages','minimum_mature_window',
        'method','gpu_uuid','dense_time_spectrum_limit','runtime_binary_sha256'):
        assert p[key]==parent[key],key
    inventory=read(ROOT/'seed-inventory.json')
    assert inventory['proposed_histories']==p['histories']==[1101,1102,1103]
    assert not set(p['histories'])&set(inventory['used_seed_or_history_ids'])
    assert inventory['all_proposed_unused'] and inventory['input_count']==len(inventory['files'])
    for record in inventory['files']:assert sha(record['path'])==record['sha256']
    assert p['independent_histories']==p['new_scientific_histories']==3
    assert p['new_numerical_paths']==len(p['jobs'])==24
    branches=['off','tension_adhesion','polarity','full'];assert p['branches']==branches
    levels=[('development',.00046875),('development_half',.000234375)]
    couplings={'off':(False,0.,0.,0.),'tension_adhesion':(True,.25,.35,0.),
        'polarity':(True,0.,0.,.35),'full':(True,.25,.35,.35)}
    expected=[(s,b,l,dt) for s in p['histories'] for b in branches for l,dt in levels]
    assert p['time_pairs']==[dict(seed=s,branch=b,cases=[f'h{s}_{b}_{l}' for l,_ in levels])
        for s in p['histories'] for b in branches]
    base=next(j for j in parent['jobs'] if j['key']=='off120')['config']
    assert p['baseline_config']==base
    integral=quad(lambda r:4*np.pi*r*r*(lambda x:x*x*(3-2*x))(expit(np.sqrt(2)*(.8-r)/.085)),
        0,.8+30*.085,epsabs=1e-12,epsrel=1e-12)[0]
    matched={};signatures={};rows=[]
    for job,(seed,branch,level,dt) in zip(p['jobs'],expected):
        assert (job['seed'],job['branch'],job['level'])==(seed,branch,level)
        assert job['key']==f'h{seed}_{branch}_{level}' and not job['reused']
        cfg=job['config'];assert cfg['seed']==seed and cfg['grid']==120 and cfg['dt']==dt
        assert cfg['steps']==round(18/dt) and cfg['save_every']==round(.15/dt)
        assert tuple(cfg[k] for k in ('feedback','fate_tension','fate_adhesion','polarity_tension'))==couplings[branch]
        assert {k for k,v in cfg.items() if v!=base[k]} <= {
            'seed','dt','steps','save_every','feedback','fate_tension','fate_adhesion','polarity_tension'}
        assert not (ROOT/job['key']).is_symlink() and not (ROOT/job['key']/'gpu-result.json').exists()
        with np.load(job['source'],allow_pickle=False) as z:
            arrays={k:z[k].copy() for k in z.files};meta=json.loads(str(z['metadata']))
            assert meta['config']==cfg and meta['time']==0 and meta['step_number']==0 and meta['next_id']==1
            assert not meta['divisions'] and meta['lineage']==[dict(id=0,parent=-1,birth=0.,division=None)]
            assert not json.loads(str(z['event_audits']))
            np.testing.assert_array_equal(z['ids'],[0]);np.testing.assert_array_equal(z['parents'],[-1])
            axis=(np.arange(120)+.5)*(2*2.24/120)-2.24
            radius=np.sqrt(axis[:,None,None]**2+axis[None,:,None]**2+axis[None,None,:]**2)
            phase=expit(np.sqrt(2)*(.8-radius)/.085).astype(np.float32)
            np.testing.assert_array_equal(z['phi'][0],phase)
            assert z['phi'].shape==(1,120,120,120) and z['phi'].dtype==np.float32
            assert abs(float(z['target'][0])/integral-1)<2e-14
            assert str(z['chemical_method'])==p['method'] and int(z['development_schema'])==1
            assert str(z['attribute_mode'])==('no_feedback' if branch=='off' else 'direct')
            np.testing.assert_array_equal(z['activator'],[1.03]);np.testing.assert_array_equal(z['inhibitor'],[.98])
            for name in ('phase_carry','chemical_correction','fate','polarity','rounding'):
                assert not z[name].any(),name
            assert z['phase_carry'].dtype==z['chemical_correction'].dtype==np.float64
            assert z['phase_carry'].shape==z['phi'].shape and z['chemical_correction'].shape==(2,1)
            # Independently regenerate all three spawned streams and first cycle.
            rngs=[np.random.default_rng(x) for x in np.random.SeedSequence(seed).spawn(3)]
            due=cfg['division_interval']*rngs[0].uniform(1-cfg['cycle_jitter'],1+cfg['cycle_jitter'])
            np.testing.assert_array_equal(z['due'],[due])
            for name,rng in zip(('rng','fate_rng','signal_rng'),rngs):assert meta[name]==rng.bit_generator.state
            del meta['config'];arrays.pop('metadata');arrays.pop('attribute_mode')
            if seed in matched:
                old_arrays,old_meta=matched[seed];assert meta==old_meta
                assert set(arrays)==set(old_arrays)
                for name in arrays:np.testing.assert_array_equal(arrays[name],old_arrays[name])
            else:matched[seed]=(arrays,meta)
            signatures[seed]=json.dumps({name:meta[name] for name in ('rng','fate_rng','signal_rng')},sort_keys=True)
        rows.append(dict(case=job['key'],seed=seed,branch=branch,dt=dt,
            analytic_phase_exact=True,initial_rng_independently_reproduced=True,matched_within_history=True))
    assert len(set(signatures.values()))==3
    gate=read(p['instrumentation_gate']);assert gate['passed']
    assert gate['protocol_sha256']==p['instrumentation_protocol_sha256']
    for file,h in gate['evidence_sha256'].items():assert sha(file)==h
    runtime=read(parent_root/'runtime.json')
    assert sha(runtime['binary'])==runtime['binary_sha256']==p['runtime_binary_sha256']
    tests=(HERE/'tests.log').read_text();assert '38 passed' in tests and 'failed' not in tests
    smoke=read(HERE/'gpu-smoke.json');assert smoke['passed'] and len(smoke['results'])==8
    assert all(r['passed'] and r['exact_gpu_restart'] for r in smoke['results'])
    migration=Path('archive/output-relocation/2026-10-10-output/studies-manifest.json')
    assert sha(migration)=='8fb2545d8bae35a1a673ec3a43d3cc191148f03cac58892f29bd368182a8aff6'
    result=dict(passed=True,utc=datetime.now(timezone.utc).isoformat(),protocol_sha256=ph,
        independent_histories=3,new_numerical_paths=24,rows=rows,
        source_count=len(p['source_sha256']),input_count=len(p['input_sha256']),
        preserved_count=len(preserved),seed_inventory_count=inventory['input_count'],
        test_log_sha256=sha(HERE/'tests.log'),gpu_smoke_sha256=sha(HERE/'gpu-smoke.json'),
        migration_record_sha256=sha(migration),runtime_binary_sha256=p['runtime_binary_sha256'],
        audit_script_sha256=sha(__file__),post_cleavage_not_launched=True)
    (HERE/'prelaunch-audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))


if __name__=='__main__':main()
