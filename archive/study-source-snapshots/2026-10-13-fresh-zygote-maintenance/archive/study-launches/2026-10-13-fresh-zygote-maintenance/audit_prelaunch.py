"""Independent NumPy/hash audit of maintenance inheritance; no embryo imports."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = Path('output/studies/fresh-zygote-maintenance')
PARENT = Path('output/studies/fresh-zygote-formation')
EXPECTED = 'e193c41f6d018cf4cfc7ed88777e00f7e62e99aa89fb98e783baab0c8e199d06'
COUPLINGS = {'off': (False,0.,0.,0.), 'tension_adhesion': (True,.25,.35,0.),
             'polarity': (True,0.,0.,.35), 'full': (True,.25,.35,.35)}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def main():
    output = HERE / 'prelaunch-audit.json'
    assert not output.exists()
    p, parent = read(ROOT / 'protocol.json'), read(PARENT / 'protocol.json')
    ph = sha(ROOT / 'protocol.json')
    assert sha(PARENT / 'protocol.json') == p['parent_protocol_sha256'] == EXPECTED
    assert read(ROOT / 'status.json') == dict(state='prepared', protocol_sha256=ph)
    print('Verifying pinned numerical sources and evidence', flush=True)
    for name, expected in {**p['source_sha256'], **p['input_sha256']}.items():
        assert sha(name) == expected, name
    for mapping in ('source_sha256', 'input_sha256'):
        assert all(p[mapping][name] == value for name, value in parent[mapping].items())
    inherited = ('histories','branches','levels','time_pairs','physical_limits','time_limits',
        'homologous_limits','dense_time_spectrum_limit','prefix_steps','prefix_limits','method',
        'gpu_uuid','runtime_binary_sha256','interval','diagnostic_interval')
    assert all(p[key] == parent[key] for key in inherited)
    assert p['histories'] == [1101,1102,1103] and p['branches'] == list(COUPLINGS)
    assert p['new_numerical_paths'] == 24 and p['new_scientific_histories'] == 0
    assert p['start_time'] == 78. and p['duration'] == 180.
    assert p['formation_origin'] == 18. and p['total_elapsed_at_end'] == 240.
    assert p['snapshots'] == [0,36,72,108,144,180] and p['restart_elapsed'] == [36,108] and p['restart_steps'] == 8
    assert p['thresholds'] == [.05,.1,.2] and p['sustained_window'] == p['late_window'] == 24.
    assert p['transition_time_limit'] == .3 and p['growth_abs_limit'] == .001
    a = read('docs/fresh_zygote_formation_assessment.json')
    rec = a['independent_review']
    assert a['review_completed'] and a['scoped_formation_timestep_accepted']
    assert sha(rec['path']) == rec['sha256']
    review = read(rec['path'])
    assert review['scoped_formation_timestep_accepted'] and review['protocol_sha256'] == a['protocol_sha256'] == EXPECTED
    count = 0
    assert len(p['jobs']) == len(parent['jobs']) == 24
    for job, old in zip(p['jobs'], parent['jobs']):
        endpoint = (PARENT / old['key'] / 'gpu-e60.npz').absolute()
        report_path = (PARENT / old['key'] / 'gpu-result.json').absolute()
        for path in (endpoint, report_path):
            assert sha(path) == review['original_files_sha256'][str(path)]
        report = read(report_path)
        assert report['passed'] and report['final_time'] == 78. and report['final_cells'] == 16
        assert report['job'] == old and report['protocol_sha256'] == EXPECTED
        assert all(x['sustained_formation'] and x['late_persistent'] for x in report['formation'])
        cfg = job['config']
        assert cfg == dict(old['config'], steps=round(258 / old['config']['dt']))
        assert cfg['dt'] in (.00375,.001875) and cfg['grid'] == 120 and cfg['extent'] == 2.24
        assert cfg['signal_transport'] == 'conservative' and cfg['max_cells'] == 16
        assert tuple(cfg[k] for k in ('feedback','fate_tension','fate_adhesion','polarity_tension')) == COUPLINGS[job['branch']]
        assert job['source_case'] == old['key'] and Path(job['source_checkpoint']).resolve() == (PARENT / old['key'] / 'gpu-e60.npz').resolve()
        assert Path(job['source']).resolve() == (ROOT / job['key'] / 'source.npz').resolve()
        assert not (ROOT / job['key'] / 'gpu-result.json').exists()
        with np.load(job['source_checkpoint'], allow_pickle=False) as z:
            arrays = {name: z[name] for name in z.files}
        with np.load(job['source'], allow_pickle=False) as z:
            current = {name: z[name] for name in z.files}
        assert set(arrays) == set(current)
        ma, mb = [json.loads(str(x['metadata'])) for x in (arrays, current)]
        assert ma.pop('config') == old['config'] and mb.pop('config') == cfg and ma == mb
        assert ma['time'] == 78. and ma['step_number'] * cfg['dt'] == 78. and not ma['divisions']
        assert current['phi'].shape == (16,120,120,120) and current['phi'].dtype == np.float32
        assert current['phase_carry'].shape == current['phi'].shape and current['phase_carry'].dtype == np.float64
        assert current['chemical_correction'].shape == (2,16) and current['chemical_correction'].dtype == np.float64
        assert np.isfinite(current['activator']).all() and current['activator'].min() > 0
        assert np.isfinite(current['inhibitor']).all() and current['inhibitor'].min() > 0
        assert float(np.log(current['activator']).std()) > max(p['thresholds'])
        assert len(json.loads(str(current['event_audits']))) == 15
        for name in set(arrays) - {'metadata'}:
            assert current[name].dtype == arrays[name].dtype and current[name].shape == arrays[name].shape
            assert current[name].tobytes() == arrays[name].tobytes(), (job['key'], name)
            count += 1
        del arrays, current
        print('Independently checked complete initial arrays:', job['key'], flush=True)
    assert count == 432 and read(ROOT / 'transfer-verification.json')['checked_arrays'] == count
    assert read('output/studies/fresh-zygote-mature-qualification/status.json')['state'] == 'failed'
    assert read('output/studies/fresh-zygote-pilot/status.json')['state'] == 'failed'
    assert not read('output/studies/fresh-zygote-pilot/time-h1103-off.json')['passed']
    tests = read(HERE / 'test-verification.json')
    assert tests['exit_code'] == 0 and sha(HERE / 'tests.log') == tests['log_sha256']
    migration = Path('archive/output-relocation/2026-10-10-output/studies-manifest.json')
    assert sha(migration) == '8fb2545d8bae35a1a673ec3a43d3cc191148f03cac58892f29bd368182a8aff6'
    runtime = read(PARENT / 'runtime.json')
    assert sha(runtime['binary']) == p['runtime_binary_sha256'] == runtime['binary_sha256']
    snapshot = Path('archive/study-document-snapshots/2026-10-13-before-fresh-zygote-maintenance/manifest.json')
    for name, expected in read(snapshot)['preserved_sha256'].items():
        # The notebook is intentionally updated after successful launch.
        if name != str(Path('docs/lab_notebook.md').absolute()):
            assert sha(name) == expected
    result = dict(passed=True, utc=datetime.now(timezone.utc).isoformat(), protocol_sha256=ph,
        source_count=len(p['source_sha256']), input_count=len(p['input_sha256']),
        losslessly_checked_start_arrays=count, all24_complete_starts_preserved=True,
        only_metadata_change='config.steps', model_couplings_and_compensation_unchanged=True,
        original_failures_retained=True, histories=p['histories'], independent_histories=3,
        inherited_start_time=78., final_physical_time=258., total_elapsed_since_t18=240.,
        readme_science_manuscript_and_ledger_unchanged=True,
        new_scientific_histories=0, new_numerical_paths=24, maintenance_extension_started=False, downstream_assays_started=False,
        test_log_sha256=sha(HERE / 'tests.log'), test_verification_sha256=sha(HERE / 'test-verification.json'),
        migration_record_sha256=sha(migration), runtime_binary_sha256=p['runtime_binary_sha256'],
        audit_script_sha256=sha(__file__))
    output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
