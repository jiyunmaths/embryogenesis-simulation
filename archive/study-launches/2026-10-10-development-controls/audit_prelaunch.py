"""Independent branch/analytic-start/provenance audit, with no embryo imports."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import re
import numpy as np
from scipy.integrate import quad
from scipy.special import expit

ROOT=Path('output/studies/development-control-qualification')
HERE=Path(__file__).parent

def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):digest.update(block)
    return digest.hexdigest()

def main():
    p=json.loads((ROOT/'protocol.json').read_text());ph=sha(ROOT/'protocol.json')
    assert json.loads((ROOT/'status.json').read_text())==dict(state='prepared',protocol_sha256=ph)
    for f,h in {**p['source_sha256'],**p['input_sha256']}.items():assert sha(f)==h,f
    preserved=json.loads(Path(p['preserved_manifest']).read_text())
    for f,h in preserved['preserved_sha256'].items():assert sha(f)==h,f
    parent_root=Path('output/studies/development-space-refinement')
    parent=json.loads((parent_root/'protocol.json').read_text())
    assert p['parent_protocol_sha256']==sha(parent_root/'protocol.json')
    for f,h in parent['source_sha256'].items():assert p['source_sha256'][f]==h,f
    accepted=json.loads(Path('docs/development_space_refinement_assessment.json').read_text())
    assert accepted['review_completed'] and accepted['scoped_screen_accepted']
    assert sha(accepted['independent_review']['path'])==accepted['independent_review']['sha256']
    assert accepted['protocol_sha256']==p['parent_protocol_sha256']
    for key in ('duration','interval','checkpoint_interval','diagnostic_interval','snapshots',
        'time_limits','space_limits','physical_limits','restart_steps','restart_stages',
        'minimum_mature_window','primary_probe','primary_order','probes','orders','calibration_limit',
        'diagnostic_spread_limit','trend_floor','dense_time_spectrum_limit','method','gpu_uuid'):
        assert p[key]==parent[key],key
    assert p['new_numerical_paths']==12 and p['reused_numerical_paths']==1
    assert p['new_scientific_histories']==0 and p['full_restart_replays']==12
    branches=['tension_adhesion','polarity','full']
    levels=[('space88',88,.00046875),('space104',104,.00046875),('space120',120,.00046875),('time120_half',120,.000234375)]
    couplings={'off':(False,0.,0.,0.),'tension_adhesion':(True,.25,.35,0.),
        'polarity':(True,0.,0.,.35),'full':(True,.25,.35,.35)}
    assert p['branches']==branches and len(p['jobs'])==13
    assert p['space_pairs']=={b:[[b+'_space88',b+'_space104'],[b+'_space104',b+'_space120']] for b in branches}
    assert p['time_pairs']=={b:[b+'_space120',b+'_time120_half'] for b in branches}
    assert p['replay_cases']==['off120']+[b+'_space120' for b in branches]
    assert p['restart_stages']==[2,4,8]
    old=p['reuse_provenance']['off120'];base=old['job']['config'];old_folder=Path(old['folder'])
    assert (ROOT/'off120').is_symlink() and (ROOT/'off120').resolve()==old_folder
    assert sha(old_folder/'gpu-result.json')==old['report_sha256']
    old_report=json.loads((old_folder/'gpu-result.json').read_text())
    assert old_report['passed'] and old_report['job']==old['job']
    assert old_report['protocol_sha256']==old['protocol_sha256']==p['parent_protocol_sha256']
    with np.load(old['job']['source'],allow_pickle=False) as z:
        old_arrays={k:z[k].copy() for k in z.files};old_meta=json.loads(str(z['metadata']))
    integral=quad(lambda r:4*np.pi*r*r*(lambda f:f*f*(3-2*f))(expit(np.sqrt(2)*(.8-r)/.085)),
        0,.8+30*.085,epsabs=1e-12,epsrel=1e-12)[0]
    expected_jobs=[(b+'_'+level,b,level,g,dt,False) for b in branches for level,g,dt in levels]+[('off120','off',old['job']['level'],120,.00046875,True)]
    rows=[]
    for job,(key,branch,level,grid,dt,reused) in zip(p['jobs'],expected_jobs):
        assert (job['key'],job['branch'],job['level'],job['reused'])==(key,branch,level,reused)
        cfg=job['config'];assert cfg['grid']==grid and cfg['dt']==dt
        assert cfg['steps']==round(18/dt) and cfg['save_every']==round(.15/dt)
        assert tuple(cfg[k] for k in ('feedback','fate_tension','fate_adhesion','polarity_tension'))==couplings[branch]
        assert {k for k,v in cfg.items() if v!=base[k]} <= {'grid','dt','steps','save_every','feedback','fate_tension','fate_adhesion','polarity_tension'}
        if not reused:assert not (ROOT/key).is_symlink() and not (ROOT/key/'gpu-result.json').exists()
        with np.load(job['source'],allow_pickle=False) as z:
            assert set(z.files)==set(old_arrays)
            meta=json.loads(str(z['metadata']));assert meta['config']==cfg
            assert {k:v for k,v in meta.items() if k!='config'}=={k:v for k,v in old_meta.items() if k!='config'}
            for k in z.files:
                if k not in ('phi','phase_carry','metadata','attribute_mode'):
                    np.testing.assert_array_equal(z[k],old_arrays[k])
            assert str(z['attribute_mode'])==('no_feedback' if reused else 'direct')
            axis=(np.arange(grid)+.5)*(2*cfg['extent']/grid)-cfg['extent']
            radius=np.sqrt(axis[:,None,None]**2+axis[None,:,None]**2+axis[None,None,:]**2)
            expected=expit(np.sqrt(2)*(.8-radius)/cfg['interface_width']).astype(np.float32)
            assert z['phi'].dtype==np.float32 and z['phi'].shape==(1,grid,grid,grid)
            np.testing.assert_array_equal(z['phi'][0],expected)
            assert z['phase_carry'].dtype==np.float64 and not z['phase_carry'].any()
            assert not z['chemical_correction'].any() and not z['polarity'].any() and not z['fate'].any()
            assert abs(float(z['target'][0])/integral-1)<2e-14
        rows.append(dict(case=key,grid=grid,dt=dt,branch=branch,reused=reused,
            analytic_phase_exact=True,matched_arrays_metadata_rng_exact=True))
    calibration=json.loads(Path(p['calibration_reference']).read_text())
    cm=json.loads(Path(p['calibration_manifest']).read_text())
    assert calibration['passed'] and len(calibration['rows'])==48
    assert calibration['protocol_sha256']==cm['protocol_sha256']==p['parent_protocol_sha256']
    assert cm['calibration_sha256']==sha(p['calibration_reference'])
    assert all(r['relative_l2']<=.0025 for r in calibration['rows'])
    combinations={(r['grid'],r['probe'],r['order']) for r in calibration['rows']}
    assert {(g,probe,order) for g in (88,104,120) for probe in (112,128) for order in (3,5)}<=combinations
    gate=json.loads(Path(p['instrumentation_gate']).read_text())
    assert gate['passed'] and gate['protocol_sha256']==p['instrumentation_protocol_sha256']
    for f,h in gate['evidence_sha256'].items():assert sha(f)==h
    runtime=json.loads((parent_root/'runtime.json').read_text())
    assert sha(runtime['binary'])==runtime['binary_sha256']==p['runtime_binary_sha256']
    testlog=(HERE/'tests.log').read_text();assert '40 passed' in testlog and 'failed' not in testlog
    import shutil
    disk=shutil.disk_usage(ROOT);assert disk.free>12*1024**3
    result=dict(passed=True,checked_at_utc=datetime.now(timezone.utc).isoformat(),protocol_sha256=ph,
        source_count=len(p['source_sha256']),input_count=len(p['input_sha256']),
        preserved_count=len(preserved['preserved_sha256']),initial_states=rows,
        inherited_calibration_rows=48,calibration_sha256=sha(p['calibration_reference']),
        inherited_instrumentation_sha256=sha(p['instrumentation_gate']),
        runtime_binary_sha256=p['runtime_binary_sha256'],test_log_sha256=sha(HERE/'tests.log'),
        script_sha256=sha(__file__),new_scientific_histories=0,free_GiB=disk.free/1024**3)
    target=HERE/'prelaunch-audit.json'
    if target.exists():raise FileExistsError(target)
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('passed','source_count','input_count','preserved_count','free_GiB')},indent=2))

if __name__=='__main__':main()
