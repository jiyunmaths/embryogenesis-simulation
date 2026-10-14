"""Independent analytic-start and provenance audit; no embryo imports."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import numpy as np
from scipy.integrate import quad
from scipy.special import expit

ROOT = Path('output/studies/development-space-refinement')
HERE = Path(__file__).parent


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''): h.update(block)
    return h.hexdigest()


def main():
    p = json.loads((ROOT/'protocol.json').read_text()); ph = sha(ROOT/'protocol.json')
    status = json.loads((ROOT/'status.json').read_text())
    assert status['state'] == 'prepared' and status['protocol_sha256'] == ph
    for path,expected in {**p['source_sha256'],**p['input_sha256']}.items(): assert sha(path) == expected,path
    preserved = json.loads(Path(p['preserved_manifest']).read_text())
    for path,expected in preserved['preserved_sha256'].items(): assert sha(path) == expected,path
    prior = json.loads(Path('output/studies/development-time-diagnosis/protocol.json').read_text())
    for path,expected in prior['source_sha256'].items(): assert p['source_sha256'][path] == expected
    assessment = json.loads(Path('docs/development_time_diagnosis_assessment.json').read_text())
    assert assessment['review_completed'] and assessment['temporal_screen_passed']
    assert sha(assessment['independent_review']) == assessment['independent_review_sha256']
    assert p['time_limits'] == prior['time_limits'] and p['physical_limits'] == prior['physical_limits']
    original = json.loads(Path('output/studies/repeated-development-gpu-qualification-v2/protocol.json').read_text())
    assert p['space_limits'] == original['space_limits']
    assert p['probes'] == original['probes'] and p['orders'] == original['orders']
    assert p['primary_probe'] == original['primary_probe'] and p['primary_order'] == original['primary_order']
    for key in ('calibration_limit','diagnostic_spread_limit','trend_floor'):
        assert p[key] == original[key]
    expected_cases = [('off72',72,.00046875,True),('off88',88,.00046875,False),
        ('off104',104,.00046875,False),('off120',120,.00046875,False),('off120_half',120,.000234375,False)]
    assert len(p['jobs']) == 5 and p['space_pairs'] == [['off72','off88'],['off88','off104'],['off104','off120']]
    assert p['time_pair'] == ['off120','off120_half'] and p['dense_time_spectrum_limit'] == .01
    provenance = p['reuse_provenance']; old_folder = Path(provenance['folder'])
    assert (ROOT/'off72').is_symlink() and (ROOT/'off72').resolve() == old_folder
    assert sha(old_folder/'gpu-result.json') == provenance['report_sha256']
    base = provenance['job']['config']; rows = []
    assert not base['feedback'] and base['fate_tension'] == base['fate_adhesion'] == base['polarity_tension'] == 0
    integral = quad(lambda r: 4*np.pi*r*r*(lambda f:f*f*(3-2*f))(expit(np.sqrt(2)*(.8-r)/.085)),
        0,.8+30*.085,epsabs=1e-12,epsrel=1e-12)[0]
    with np.load(provenance['job']['source'],allow_pickle=False) as z:
        old = {k:z[k] for k in z.files}; old_meta = json.loads(str(old['metadata']))
    for job,(key,grid,dt,reused) in zip(p['jobs'],expected_cases):
        assert job['key'] == key and job['reused'] == reused
        cfg = job['config']; assert cfg['grid'] == grid and cfg['dt'] == dt
        assert cfg['steps'] == round(18/dt) and cfg['save_every'] == round(.15/dt)
        changed = {k for k,v in cfg.items() if v != base[k]}
        assert changed <= {'grid','dt','steps','save_every'}
        folder = ROOT/key
        if not reused: assert not folder.is_symlink() and not (folder/'gpu-result.json').exists()
        with np.load(job['source'],allow_pickle=False) as z:
            assert set(z.files) == set(old)
            meta = json.loads(str(z['metadata'])); assert meta['config'] == cfg
            assert {k:v for k,v in meta.items() if k != 'config'} == {k:v for k,v in old_meta.items() if k != 'config'}
            for name in z.files:
                if name not in ('phi','phase_carry','metadata'): np.testing.assert_array_equal(z[name],old[name])
            phi = z['phi']; assert phi.dtype == np.float32 and phi.shape == (1,grid,grid,grid)
            axis = (np.arange(grid)+.5)*(2*cfg['extent']/grid)-cfg['extent']
            radius = np.sqrt(axis[:,None,None]**2+axis[None,:,None]**2+axis[None,None,:]**2)
            expected = expit(np.sqrt(2)*(.8-radius)/cfg['interface_width']).astype(np.float32)
            np.testing.assert_array_equal(phi[0],expected)
            assert z['phase_carry'].dtype == np.float64 and not z['phase_carry'].any()
            assert not z['chemical_correction'].any() and not z['polarity'].any() and not z['fate'].any()
            assert abs(float(z['target'][0])/integral-1) < 2e-14
            h = phi*phi*(np.float32(3)-np.float32(2)*phi)
            volume = float(h.sum(dtype=float))*(2*cfg['extent']/grid)**3
            rows.append(dict(case=key,grid=grid,dt=dt,reused=reused,analytic_phase_exact=True,
                common_arrays_metadata_rng_exact=True,initial_relative_volume_error=abs(volume/integral-1)))
    calibration = json.loads((ROOT/'calibration.json').read_text())
    cm = json.loads((ROOT/'calibration-manifest.json').read_text())
    assert calibration['passed'] and len(calibration['rows']) == 48
    assert calibration['protocol_sha256'] == cm['protocol_sha256'] == ph
    assert cm['calibration_sha256'] == sha(ROOT/'calibration.json')
    assert all(i['relative_l2'] <= .0025 for i in calibration['rows'])
    assert {(r['grid'],r['probe'],r['order']) for r in calibration['rows']} == {
        (g,pr,o) for g in (72,88,104,120) for pr in (112,128) for o in (3,5)}
    gate = json.loads(Path(p['instrumentation_gate']).read_text())
    assert gate['passed'] and gate['protocol_sha256'] == provenance['protocol_sha256']
    assert [(r['stage'],r['steps'],r['bitwise_equal']) for r in gate['rows']] == [(2,8,True),(4,8,True),(8,8,True)]
    for path,expected in gate['evidence_sha256'].items(): assert sha(path) == expected
    assert '29 passed' in (HERE/'tests.log').read_text()
    result = dict(passed=True,checked_at_utc=datetime.now(timezone.utc).isoformat(),protocol_sha256=ph,
        source_count=len(p['source_sha256']),input_count=len(p['input_sha256']),
        preserved_count=len(preserved['preserved_sha256']),initial_states=rows,
        calibration_rows=48,calibration_max=max(i['relative_l2'] for i in calibration['rows']),
        calibration_sha256=sha(ROOT/'calibration.json'),inherited_instrumentation_gate_sha256=sha(p['instrumentation_gate']),
        test_log_sha256=sha(HERE/'tests.log'),script_sha256=sha(__file__),new_scientific_histories=0)
    target = HERE/'prelaunch-audit.json'
    if target.exists(): raise FileExistsError(target)
    target.write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps(result,indent=2))


if __name__ == '__main__': main()
