"""Read-only NumPy/SciPy audit of the stopped fresh-zygote pilot.

Reuse independently written saved-array checks, never the embryo comparison
helpers or GPU stepping. A verified numerical failure remains a failure.
"""
from datetime import datetime, timezone
from pathlib import Path
import importlib.util
import json
import re
import time

import numpy as np
from scipy.linalg import eigh

ROOT = Path('output/studies/fresh-zygote-pilot')
HERE = Path(__file__).parent
HELPER = Path('archive/study-reviews/2026-10-11-development-controls/review_saved_arrays.py')
EXPECTED = '728f584676ea4eeeee530656015aec61931e9296ff095eff8b7c94f9832a13ee'
spec = importlib.util.spec_from_file_location('independent_controls_review', HELPER)
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)
c.ROOT = ROOT
s, t, r = c.s, c.t, c.r


def spectrum(g, v):
    return eigh(np.diag(g.sum(1)) - g, np.diag(v), eigvals_only=True)


def unthresholded(row, width):
    w = np.array(row['contact']); centers = np.array(row['centers'])
    distances = np.linalg.norm(centers[:, None] - centers[None, :], axis=-1)
    g = np.divide(6 * np.sqrt(2) * w / width, distances,
                  out=np.zeros_like(w), where=distances > 1e-12)
    np.fill_diagonal(g, 0.)
    return g


def cutoff_diagnostic(a, b, cfg):
    records = []; unequal = []
    for stamp in sorted(a):
        x, y = a[stamp], b[stamp]
        if len(x['ids']) != len(y['ids']):
            unequal.append(stamp)
            continue
        sx, sy = r.eig(x)[0], r.eig(y)[0]
        ga, gb = unthresholded(x, cfg['interface_width']), unthresholded(y, cfg['interface_width'])
        # Check the geometry conversion without calling production transport.
        for row, raw in ((x, ga), (y, gb)):
            mask = np.array(row['contact']) > row['cutoff']
            np.fill_diagonal(mask, False)
            np.testing.assert_allclose(raw * mask, row['conductance'], rtol=2e-12, atol=1e-13)
        ua, ub = spectrum(ga, x['volumes']), spectrum(gb, y['volumes'])
        records.append(dict(time=stamp, cells=len(sx), edges=[x['edges'], y['edges']],
            relative_spectrum=float(np.linalg.norm(sx-sy)/(np.linalg.norm(sy) or 1.)),
            posthoc_no_cutoff_relative_spectrum=float(np.linalg.norm(ua-ub)/(np.linalg.norm(ub) or 1.))))
    peak = max(records, key=lambda row: row['relative_spectrum'])
    stamp = peak['time']; x, y = a[stamp], b[stamp]
    assert x['ids'] == y['ids'] and x['parents'] == y['parents']
    ga, gb = np.array(x['conductance']), np.array(y['conductance'])
    ua, ub = unthresholded(x, cfg['interface_width']), unthresholded(y, cfg['interface_width'])
    differences = []
    for i, j in zip(*np.triu_indices(len(ga), 1)):
        if (ga[i, j] > 0) != (gb[i, j] > 0):
            differences.append(dict(cells=[x['ids'][i], x['ids'][j]],
                contact=[x['contact'][i][j], y['contact'][i][j]],
                cutoff=[x['cutoff'], y['cutoff']],
                cutoff_margins=[x['contact'][i][j]-x['cutoff'], y['contact'][i][j]-y['cutoff']],
                conductance=[float(ga[i, j]), float(gb[i, j])]))
    masks = {'intersection': (ga > 0) & (gb > 0), 'union': (ga > 0) | (gb > 0)}
    matched = {}
    for key, mask in masks.items():
        ma, mb = spectrum(ua * mask, x['volumes']), spectrum(ub * mask, y['volumes'])
        matched[key] = float(np.linalg.norm(ma-mb)/(np.linalg.norm(mb) or 1.))
    return dict(original_peak=peak, above_limit_times=[row['time'] for row in records if row['relative_spectrum'] > .01],
        posthoc_no_cutoff_peak=max(records, key=lambda row: row['posthoc_no_cutoff_relative_spectrum']),
        original_unequal_count_times=unequal, same_cell_ids_and_parents_at_peak=True,
        peak_ids=x['ids'], differing_edges=differences, matched_support_peak_errors=matched,
        peak_relative_volume_error=float(abs(np.array(x['volumes'])/y['volumes']-1).max()),
        peak_active_furrows=[len(x['active_divisions']), len(y['active_divisions'])],
        furrow_start_time_differences=[dict(parent=ev['parent'], difference=ev['start']-other['start'])
            for ev in x['active_divisions'] for other in y['active_divisions'] if ev['parent']==other['parent']],
        time_series=records,
        scope='Post hoc algebraic reassembly of saved geometries only. No chemistry or mechanics rerun; no threshold, gate or recorded outcome changed. Does not establish a new qualified solver.')


def main():
    started = time.perf_counter(); output = HERE/'independent_review.json'
    assert not output.exists(), 'Never overwrite a recorded review'
    for name in (Path(__file__), HELPER, c.HELPER, s.HELPER, t.HELPER):
        r.pin(name)
    p = r.read(ROOT/'protocol.json'); ph = r.pin(ROOT/'protocol.json', EXPECTED)
    status = r.read(ROOT/'status.json')
    assert status['state']=='failed' and status['completed']==18 and status['total']==24
    assert status['protocol_sha256']==ph and not status['post_cleavage_started']
    assert 'Developmental timestep gate failed' in status['failure']['message']
    assert not (ROOT/'summary.json').exists()
    for name, expected in {**p['source_sha256'], **p['input_sha256']}.items():
        r.pin(name, expected)
    runtime = r.read(ROOT/'runtime.json'); r.pin(runtime['binary'], runtime['binary_sha256'])
    assert runtime['binary_sha256']==p['runtime_binary_sha256']
    recorded = {row['case']: row for row in status['paths']}
    paths = []; reports = {}; histories = {}; networks = {}; jobs = {}; prefixes = []
    for job in p['jobs']:
        key = job['key']
        if key not in recorded:
            assert not (ROOT/key/'gpu-result.json').exists()
            assert not (ROOT/key/'actual-start-prefix.json').exists()
            continue
        check, report, history, network = c.review_path(job, p, ph)
        assert report['passed']==recorded[key]['passed']
        r.close(report['seconds'], recorded[key]['seconds'], ('coordinator runtime', key))
        prefix = r.read(ROOT/key/'actual-start-prefix.json')
        assert prefix['passed'] and prefix['failure'] is None and prefix['steps']==p['prefix_steps']
        assert prefix['protocol_sha256']==ph and prefix['exact_gpu_restart']
        for name, value in prefix['errors'].items(): assert value <= p['prefix_limits'][name]
        for name, expected in prefix['evidence_sha256'].items(): r.pin(name, expected)
        prefixes.append(dict(case=key, recorded_reference_passed=True, recorded_errors=prefix['errors']))
        paths.append(check); reports[key]=report; histories[key]=history; networks[key]=network; jobs[key]=job
    assert len(paths)==len(prefixes)==18
    pairs=[]
    for original in status['pairs']:
        file=ROOT/f"time-h{original['seed']}-{original['branch']}.json"
        saved=r.read(file); assert saved==original
        for name, expected in original['evidence_sha256'].items(): r.pin(name, expected)
        row=s.compare_histories(original['temporal'],histories,reports,networks)
        dense=row['dense']; old_dense=original['dense']
        r.close(dense['peak']['relative_spectrum'],old_dense['peak']['relative_spectrum'],'dense spectrum')
        assert dense['peak']['time']==old_dense['peak']['time']
        assert dense['unmatched_times']==old_dense['unmatched_times']
        assert dense['topology_mismatch_times']==old_dense['topology_mismatch_times']
        failures={key:value for key,value in row['errors'].items() if value>p['time_limits'][key]}
        assert not failures and row['maturity_time_error']<=p['time_limits']['maturity_time_difference_max']
        if dense['peak']['relative_spectrum']>p['dense_time_spectrum_limit']:
            failures['dense_relative_spectrum']=dense['peak']['relative_spectrum']
        assert bool(failures)==(not original['passed'])
        pairs.append(dict(history=original['seed'],branch=original['branch'],cases=original['cases'],
            independently_recomputed_regular_errors=row['errors'],maturity_time_error=row['maturity_time_error'],
            independently_recomputed_dense=dense,recorded_field_max=original['temporal']['field_max'],
            density_fields_recomputed=False,failed_gates=failures,passed=not failures))
    assert len(pairs)==9 and sum(row['passed'] for row in pairs)==8
    # Validate all 18 actual endpoints; recompute all seven checkpoint fields for the failed pair.
    checkpoints=[]; endpoints=[]; failed_field_rows=[]
    bad=('h1103_off_development','h1103_off_development_half')
    for stamp in (0,3,6,9,12,15,18):
        fields={}; ranges={}
        keys=list(jobs) if stamp==18 else list(bad)
        for key in keys:
            z=r.load(ROOT/key/f'gpu-t{stamp:02d}.npz')
            cfg=jobs[key]['config']
            check=r.checkpoint(z,cfg,histories[key][round(stamp/.15)],stamp==18)
            assert check['volume_error']<p['physical_limits']['volume_error_max']
            assert check['minimum_radius']>=p['physical_limits']['radius_min']
            assert z['phase_carry'].dtype==np.float64 and z['phase_carry'].shape==z['phi'].shape
            assert z['chemical_correction'].dtype==np.float64 and z['chemical_correction'].shape==(2,len(z['ids']))
            assert np.isfinite(z['phase_carry']).all() and np.isfinite(z['chemical_correction']).all()
            assert json.loads(str(z['event_audits']))==[e for e in reports[key]['events'] if e['time']<=stamp+1e-12]
            checkpoints.append(dict(case=key,**check,carry_and_chemical_corrections_present=True))
            if key in bad: fields[key],ranges[key]=s.density(z,None)
            if stamp==18:
                h=histories[key][-1]; lam,_=r.eig(h)
                beta,da,db=cfg['signal_beta'],cfg['signal_da'],cfg['signal_dh']
                jac=np.array([[1.,-1.],[2*beta,-beta]])
                growth=np.array([np.linalg.eigvals(jac-x*np.diag([da,db])).real.max() for x in lam])
                endpoints.append(dict(case=key,history=jobs[key]['seed'],branch=jobs[key]['branch'],level=jobs[key]['level'],
                    time=h['time'],cells=len(h['ids']),axis_ratio=h['axis_ratio'],log_activator_sd=h['log_activator_sd'],
                    frozen_uniform_unstable_modes=int(np.sum((lam>1e-10)&(growth>1e-10))),
                    max_frozen_uniform_growth=float(growth[lam>1e-10].max()),
                    active_furrows=len(h['dividing_ids']),chemistry_range=[[min(a),max(a)] for a in h['chemistry']]))
            del z
        values=r.field_error(fields[bad[0]],fields[bad[1]])
        expected=next(row for row in status['pairs'][-1]['temporal']['fields'] if row['time']==stamp)
        for key,value in values.items(): r.close(value,expected['errors'][key],('failed-pair field',stamp,key))
        np.testing.assert_allclose([ranges[key] for key in bad],expected['phase_ranges'],rtol=1e-6,atol=1e-9)
        failed_field_rows.append(dict(time=stamp,errors=values,phase_ranges=[ranges[key] for key in bad]))
        del fields
        print('Reviewed failed-pair fields and saved arrays:',stamp,flush=True)
    maximum=max(value for row in failed_field_rows for value in row['errors'].values())
    r.close(maximum,status['pairs'][-1]['temporal']['field_max'],'failed-pair field max')
    assert maximum<=p['time_limits']['field_max']
    pairs[-1]['density_fields_recomputed']=True
    diagnostic=cutoff_diagnostic(networks[bad[0]],networks[bad[1]],jobs[bad[0]]['config'])
    assert diagnostic['above_limit_times']==[11.43]
    # Check the authorized editorial exceptions separately from the numerical failure.
    record_path=Path('archive/study-document-snapshots/2026-10-11-authorized-comment-resolution/record.json')
    record=r.read(record_path); allowed={row['original_path']:row for row in record['non_input_document_relocations']}
    manifest=r.read(p['preserved_manifest']); unchanged=0; reconciled=[]; annotated=[]
    for name,expected in manifest['preserved_sha256'].items():
        actual=r.pin(name)
        if actual==expected: unchanged+=1; continue
        entry=allowed[name]
        assert expected==entry['expected_sha256']
        r.pin(entry['preserved_bytes'],expected)
        assert name not in p['source_sha256'] and name not in p['input_sha256']
        if actual!=entry['edited_sha256']:
            assert Path(name)==Path('manuscript/sections.md').resolve()
            text=Path(name).read_text(); comments=[]
            for match in re.finditer(r'\[([^\]\n]+)\]',text):
                comment=match.group(1)
                if comment.startswith(('what does the related','is there a citation','what does \'leave','what is the animal','what are the frozen')):
                    comments.append(comment)
            assert len(comments)==5
            r.pin(HERE/'sections-observed-with-new-comments.md',actual)
            annotated.append(dict(path=name,current_sha256=actual,
                last_resolved_sha256=entry['edited_sha256'],comments=comments,
                source_not_edited_by_this_review=True,new_comments_unaddressed=True,
                current_draft_matches_last_editorial_record=False,
                full_current_editorial_delta_reconciliation_pending=True,
                archived_original_preservation_hash_verified=True))
        reconciled.append(dict(entry,current_observed_sha256=actual))
    assert unchanged==1769 and len(reconciled)==5
    for name,expected in list(r.WATCHED.items()): assert r.sha(name)==expected,name
    result=dict(review_completed=True,qualification_accepted=False,developmental_design_completed=False,
        numerical_stop_confirmed=True,reviewed_utc=datetime.now(timezone.utc).isoformat(),seconds=time.perf_counter()-started,
        protocol_sha256=ph,registered_histories=[1101,1102,1103],full_matched_histories=[1101,1102],
        completed_paths=18,total_paths=24,completed_pairs=9,passed_pairs=8,total_pairs=12,
        source_count=len(p['source_sha256']),input_count=len(p['input_sha256']),paths=paths,prefixes=prefixes,
        pairs=pairs,checkpoints=checkpoints,endpoints=endpoints,failed_pair_field_rows=failed_field_rows,
        cutoff_diagnostic=diagnostic,unexecuted_cases=[job['key'] for job in p['jobs'] if job['key'] not in jobs],
        preserved_history=dict(unchanged_count=unchanged,authorized_editorial_count=len(reconciled),
            authorization_record=str(record_path),original_archived_bytes_verified=True,exceptions=reconciled,
            additional_user_annotation_deltas=annotated,
            all_current_drafts_match_latest_editorial_record=False,
            caused_observed_stop=False,original_manifest_unchanged=True),
        scope='Independent NumPy/SciPy saved-evidence audit of a stopped experiment. All completed histories, graph/pressure/event records and hashes checked; all 18 endpoints and all seven failed-pair density checkpoints checked. Other eight pairs\u2019 full density histories, full GPU/CPU mechanics and long restart replays are not rerun. No production solver imports or scientific simulations.',
        biological_identity_qualified=False,general_continuum_qualified=False,post_cleavage_started=False,
        original_files_sha256=r.WATCHED)
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:result[key] for key in ('review_completed','qualification_accepted','completed_paths','passed_pairs','numerical_stop_confirmed','seconds')},indent=2),flush=True)


if __name__=='__main__': main()
