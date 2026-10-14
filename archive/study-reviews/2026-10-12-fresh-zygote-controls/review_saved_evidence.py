"""Independent NumPy/SciPy audit of history-1103's finer matched controls.

Rebuild comparisons from saved states, graphs and events; no embryo imports,
production comparison calls, fresh GPU stepping or mutation of original runs.
"""
from datetime import datetime, timezone
from pathlib import Path
import importlib.util
import json
import time

import numpy as np

ROOT = Path('output/studies/fresh-zygote-control-completion')
PARENT = Path('output/studies/fresh-zygote-time-refinement')
ORIGINAL = Path('output/studies/fresh-zygote-pilot')
HERE = Path(__file__).resolve().parent
EXPECTED = '44af9b2bfdd77176891e015e5024db2aeca6acca413f7d8dd63e5b8a86ffd926'
HELPER = Path('archive/study-reviews/2026-10-11-development-controls/review_saved_arrays.py')
BRANCHES = ('off', 'tension_adhesion', 'polarity', 'full')
COUPLINGS = {'off': (False,0.,0.,0.), 'tension_adhesion': (True,.25,.35,0.),
    'polarity': (True,0.,0.,.35), 'full': (True,.25,.35,.35)}
LEVELS = (('development_half', .000234375), ('development_quarter', .0001171875))


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec); spec.loader.exec_module(obj)
    return obj


c = module('independent_control_review', HELPER)
c.ROOT = ROOT
s, t, r = c.s, c.t, c.r


def endpoint(history, cfg):
    row = history[-1]; lam, _ = r.eig(row)
    u = np.array(row['chemistry']); beta = cfg['signal_beta']
    jac = np.array([[1., -1.], [2*beta, -beta]])
    growth = np.array([np.linalg.eigvals(jac - ell*np.diag([cfg['signal_da'], cfg['signal_dh']])).real.max() for ell in lam])
    sd = float(np.log(u[0]).std()); r.close(sd, row['log_activator_sd'], 'endpoint contrast')
    return dict(time=18, cells=len(row['ids']), axis_ratio=row['axis_ratio'],
        log_activator_sd=sd, chemistry_ranges=[[float(v.min()),float(v.max())] for v in u],
        frozen_uniform_unstable_modes=int(np.sum((lam>1e-10)&(growth>1e-10))),
        max_frozen_uniform_growth=float(growth[lam>1e-10].max()), active_furrows=len(row['dividing_ids']))


def main():
    started = time.perf_counter(); output = HERE/'independent_review.json'
    assert not output.exists(), 'Never overwrite a completed review'
    for name in (Path(__file__), HELPER, c.HELPER, s.HELPER, t.HELPER):r.pin(name)
    p = r.read(ROOT/'protocol.json'); ph = r.pin(ROOT/'protocol.json', EXPECTED)
    status = r.read(ROOT/'status.json'); summary = r.read(ROOT/'summary.json')
    assert status['state']=='completed_provisional' and status['completed']==status['total']==6
    assert status['pairs_completed']==status['pairs_total']==3
    assert summary['execution_completed'] and summary['developmental_screens_passed']
    assert summary['independent_review_pending'] and not summary['qualification_accepted']
    assert status['protocol_sha256']==summary['protocol_sha256']==ph
    assert summary['new_scientific_histories']==0 and summary['new_numerical_paths']==6
    assert summary['reused_numerical_paths']==2 and not summary['post_cleavage_started']
    for name, expected in {**p['source_sha256'], **p['input_sha256']}.items():r.pin(name, expected)
    parent = r.read(PARENT/'protocol.json')
    assert r.pin(PARENT/'protocol.json',p['parent_protocol_sha256'])=='54632a1f284affcbcec93cd0f1036d2d9c3788bd93ef5480df5fa1f14336c6c3'
    inherited=('duration','interval','checkpoint_interval','diagnostic_interval','snapshots',
        'time_limits','physical_limits','restart_steps','restart_stages','minimum_mature_window',
        'method','gpu_uuid','dense_time_spectrum_limit','instrumentation_gate',
        'instrumentation_protocol_sha256','runtime_binary_sha256','prefix_steps','prefix_limits')
    for key in inherited:assert p[key]==parent[key],key
    assert p['dense_time_spectrum_limit']==p['time_limits']['relative_spectrum_max']==p['time_limits']['field_max']==.01
    assert p['new_scientific_histories']==0 and p['new_numerical_paths']==6 and p['reused_numerical_paths']==2
    accepted = r.read(p['accepted_assessment']); review_record = accepted['independent_review']
    prior_off = r.read(review_record['path']);r.pin(review_record['path'],review_record['sha256'])
    assert accepted['scoped_refinement_accepted'] and prior_off['scoped_refinement_accepted']
    old_assessment = r.read('docs/fresh_zygote_pilot_assessment.json')
    record = old_assessment['independent_saved_evidence_review']
    old_review = r.read(record['path']);r.pin(record['path'],record['sha256'])
    old_status = r.read(ORIGINAL/'status.json');old_pair = r.read(ORIGINAL/'time-h1103-off.json')
    assert old_status['state']=='failed' and old_status['completed']==18 and not old_pair['passed']
    assert old_review['numerical_stop_confirmed'] and not old_review['qualification_accepted']
    assert len(old_review['unexecuted_cases'])==6
    runtime = r.read(ROOT/'runtime.json');r.pin(runtime['binary'],runtime['binary_sha256'])
    assert runtime['binary_sha256']==p['runtime_binary_sha256']
    expected_keys=[f'h1103_{b}_{l}' for b in BRANCHES for l,_ in LEVELS]
    assert [j['key'] for j in p['jobs']]==expected_keys
    jobs={j['key']:j for j in p['jobs']};paths=[];reports={};histories={};networks={};prefixes=[]
    initial=None;base=None
    for job in p['jobs']:
        key=job['key'];cfg=job['config'];assert job['seed']==cfg['seed']==1103 and cfg['grid']==120
        assert job['reused']==(job['branch']=='off') and cfg['graph_contact_cutoff']==.02
        assert cfg['dt']==dict(LEVELS)[job['level']] and cfg['steps']*cfg['dt']==18
        assert cfg['save_every']*cfg['dt']==.15
        assert tuple(cfg[k] for k in ('feedback','fate_tension','fate_adhesion','polarity_tension'))==COUPLINGS[job['branch']]
        if base is None:base=cfg
        assert {k for k,v in cfg.items() if v!=base[k]} <= {'dt','steps','save_every','feedback','fate_tension','fate_adhesion','polarity_tension'}
        z=r.load(job['source']);meta=json.loads(str(z['metadata']));assert meta.pop('config')==cfg
        assert meta['time']==0 and meta['step_number']==0 and not meta['divisions']
        if initial is None:initial=(meta,{k:z[k].copy() for k in z if k not in ('metadata','attribute_mode')})
        assert meta==initial[0] and set(z)-{'metadata','attribute_mode'}==set(initial[1])
        for name,array in initial[1].items():np.testing.assert_array_equal(z[name],array,err_msg=(key,name))
        assert str(z['attribute_mode'])==('no_feedback' if job['branch']=='off' else 'direct');del z
        checked,report,history,network=c.review_path(job,p,ph)
        paths.append(checked);reports[key]=report;histories[key]=history;networks[key]=network
        prefix=r.read(ROOT/key/'actual-start-prefix.json')
        assert prefix['passed'] and prefix['failure'] is None and prefix['steps']==8 and prefix['exact_gpu_restart']
        assert prefix['protocol_sha256']==report['protocol_sha256']
        for name,value in prefix['errors'].items():assert value<=p['prefix_limits'][name]
        for name,expected in prefix['evidence_sha256'].items():r.pin(name,expected)
        if not job['reused']:assert prefix==next(x for x in summary['prefixes'] if x['case']==key)
        prefixes.append(dict(case=key,recorded_passed=True,recorded_errors=prefix['errors']))
    del initial
    pairs=[];samples=[];live=[];endpoints=[]
    for branch in BRANCHES:
        cases=[f'h1103_{branch}_{level}' for level,_ in LEVELS]
        original=r.read((PARENT if branch=='off' else ROOT)/f'time-h1103-{branch}.json')
        assert original['passed'] and original['cases']==cases
        if branch!='off':assert original==next(x for x in summary['pairs'] if x['branch']==branch)
        row=s.compare_histories(original['temporal'],histories,reports,networks)
        dense=row['dense'];saved=original['dense']
        r.close(dense['peak']['relative_spectrum'],saved['peak']['relative_spectrum'],'dense peak')
        assert dense['peak']['time']==saved['peak']['time'] and dense['peak']['edges']==saved['peak']['edges']
        assert dense['unmatched_times']==saved['unmatched_times'] and dense['topology_mismatch_times']==saved['topology_mismatch_times']
        for stamp in (0,3,6,9,12,15,18):
            fields={};ranges={}
            for key in cases:
                z=r.load(ROOT/key/f'gpu-t{stamp:02d}.npz')
                check=r.checkpoint(z,jobs[key]['config'],histories[key][round(stamp/.15)],stamp==18)
                assert check['volume_error']<p['physical_limits']['volume_error_max'] and check['minimum_radius']>=p['physical_limits']['radius_min']
                assert json.loads(str(z['event_audits']))==[e for e in reports[key]['events'] if e['time']<=stamp+1e-12]
                samples.append(dict(case=key,**check,corrections_present=True))
                fields[key],ranges[key]=s.density(z,None)
                if stamp==18:endpoints.append(dict(case=key,history=1103,branch=branch,level=jobs[key]['level'],dt=jobs[key]['config']['dt'],**endpoint(histories[key],jobs[key]['config'])))
                del z
            errors=r.field_error(fields[cases[0]],fields[cases[1]])
            expected=next(x for x in original['temporal']['fields'] if x['time']==stamp)
            for name,value in errors.items():r.close(value,expected['errors'][name],('fields',branch,stamp,name))
            np.testing.assert_allclose([ranges[key] for key in cases],expected['phase_ranges'],rtol=1e-6,atol=1e-9)
            row['field_checks'].append(dict(time=stamp,probe=None,order=3,errors=errors));del fields
            print('Independently reconstructed density fields:',branch,stamp,flush=True)
        for key in cases:
            for stage in (2,4,8):
                z=r.load(ROOT/key/f'restart-stage-{stage}.npz');check=r.checkpoint(z,jobs[key]['config'])
                assert stage<check['cells']<2*stage and check['active']>0
                live.append(dict(case=key,stage=stage,**check));del z
        s.decision(row,original['temporal'],p['time_limits'])
        assert row['independently_passed'] and dense['peak']['relative_spectrum']<=p['dense_time_spectrum_limit']
        pairs.append(dict(history=1103,branch=branch,**row))
    # Earlier accepted pairs retained recorded field errors. Rebuild those fields
    # now before making a combined numerical-preparation statement.
    cohort=[];cohort_endpoints=[]
    old_p=r.read(ORIGINAL/'protocol.json')
    assert old_p['time_limits']==p['time_limits'] and old_p['physical_limits']==p['physical_limits']
    for pair in old_review['pairs']:
        if pair['history'] not in (1101,1102):continue
        assert pair['passed']
        old_jobs={};old_reports={};old_histories={};old_networks={}
        for case in pair['cases']:
            job=next(j for j in old_p['jobs'] if j['key']==case);old_jobs[case]=job
            report=r.read(ORIGINAL/case/'gpu-result.json')
            for name,expected in report['evidence_sha256'].items():r.pin(name,expected)
            old_reports[case]=report;old_histories[case]=r.read(ORIGINAL/case/'gpu-history.json')
            records=t.jsonl(ORIGINAL/case/'network.jsonl')
            old_networks[case]={round(x['time'],9):x for x in records if abs(x['time']/.015-round(x['time']/.015))<1e-7}
        original=r.read(ORIGINAL/f"time-h{pair['history']}-{pair['branch']}.json")
        row=s.compare_histories(original['temporal'],old_histories,old_reports,old_networks)
        r.close(row['dense']['peak']['relative_spectrum'],original['dense']['peak']['relative_spectrum'],'cohort dense')
        assert row['dense']['unmatched_times']==original['dense']['unmatched_times']
        assert row['dense']['topology_mismatch_times']==original['dense']['topology_mismatch_times']
        for stamp in (0,3,6,9,12,15,18):
            fields={};ranges={}
            for case in pair['cases']:
                z=r.load(ORIGINAL/case/f'gpu-t{stamp:02d}.npz')
                check=r.checkpoint(z,old_jobs[case]['config'],old_histories[case][round(stamp/.15)],stamp==18)
                assert check['volume_error']<p['physical_limits']['volume_error_max'] and check['minimum_radius']>=p['physical_limits']['radius_min']
                samples.append(dict(case=case,**check,corrections_present=True,earlier_study=True))
                fields[case],ranges[case]=s.density(z,None);del z
            a,b=pair['cases'];errors=r.field_error(fields[a],fields[b])
            expected=next(x for x in original['temporal']['fields'] if x['time']==stamp)
            for name,value in errors.items():r.close(value,expected['errors'][name],('cohort fields',a,b,stamp,name))
            np.testing.assert_allclose([ranges[a],ranges[b]],expected['phase_ranges'],rtol=1e-6,atol=1e-9)
            row['field_checks'].append(dict(time=stamp,probe=None,order=3,errors=errors));del fields
        s.decision(row,original['temporal'],p['time_limits'])
        assert row['independently_passed'] and row['dense']['peak']['relative_spectrum']<=p['dense_time_spectrum_limit']
        cohort.append(dict(history=pair['history'],branch=pair['branch'],passed=True,
            dts=[.00046875,.000234375],primary_dt=.000234375,
            prior_review=str(record['path']),dense_peak=pair['independently_recomputed_dense']['peak']['relative_spectrum'],
            field_max=row['field_max'],evidence_reused=True,density_fields_recomputed=True,
            independent_comparison=row))
        print('Reconstructed earlier accepted pair:',pair['history'],pair['branch'],flush=True)
    assert len(cohort)==8
    for e in old_review['endpoints']:
        if e['history'] not in (1101,1102) or e['level']!='development_half':continue
        key=e['case'];old_p=r.read(ORIGINAL/'protocol.json');job=next(j for j in old_p['jobs'] if j['key']==key)
        history=r.read(ORIGINAL/key/'gpu-history.json');z=r.load(ORIGINAL/key/'gpu-t18.npz')
        r.checkpoint(z,job['config'],history[-1],True);del z
        measured=endpoint(history,job['config'])
        for name in ('log_activator_sd','axis_ratio','frozen_uniform_unstable_modes'):r.close(measured[name],e[name],('cohort endpoint',key,name))
        cohort_endpoints.append(dict(case=key,history=e['history'],branch=e['branch'],dt=.000234375,**measured))
    assert len(cohort_endpoints)==8
    for pair in pairs:
        cohort.append(dict(history=1103,branch=pair['branch'],passed=True,
            dts=[.000234375,.0001171875],primary_dt=.000234375,
            dense_peak=pair['dense']['peak']['relative_spectrum'],field_max=pair['field_max'],evidence_reused=pair['branch']=='off'))
    cohort_endpoints += [e for e in endpoints if e['level']=='development_half']
    assert {(x['history'],x['branch']) for x in cohort}=={(h,b) for h in (1101,1102,1103) for b in BRANCHES}
    assert len(cohort)==len(cohort_endpoints)==12
    for name,expected in list(r.WATCHED.items()):assert r.sha(name)==expected,name
    result=dict(review_completed=True,scoped_control_completion_accepted=True,
        scoped_combined_developmental_preparations_accepted=True,original_pilot_accepted=False,
        reviewed_utc=datetime.now(timezone.utc).isoformat(),seconds=time.perf_counter()-started,
        protocol_sha256=ph,source_count=len(p['source_sha256']),input_count=len(p['input_sha256']),
        new_scientific_histories=0,new_numerical_paths=6,reused_numerical_paths=2,
        combined_registered_histories=[1101,1102,1103],combined_history_branch_preparations=12,
        paths=paths,prefixes=prefixes,pairs=pairs,checkpoints=samples,live_checkpoints=live,
        endpoints=endpoints,combined_preparations=cohort,combined_primary_endpoints=cohort_endpoints,
        original_failure_retained=True,original_unexecuted_cases=old_review['unexecuted_cases'],
        biological_identity_qualified=False,general_continuum_qualified=False,
        general_developmental_convergence_accepted=False,post_cleavage_started=False,
        scope='Independent NumPy/SciPy saved-array audit: all eight history1103 paths, regular histories, dense graphs/pressure/events, 56 regular checkpoint states/densities and 24 live-cleavage states checked. Actual-start prefixes and exact eight-step restarts verified from saved evidence, not rerun. Earlier accepted history1101/1102 pairs linked through original review and immutable evidence, with regular/dense comparisons and all 112 checkpoint states/density reconstructions now independently repeated. Twelve primary endpoint arrays checked across the combined preparations. Different timestep pairs are declared by history. No fresh mechanics replay, full fine CPU horizon, long-window organization, universal convergence or original failed pilot acceptance.',
        original_files_sha256=dict(r.WATCHED))
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('review_completed','scoped_control_completion_accepted','combined_history_branch_preparations','seconds')},indent=2),flush=True)


if __name__=='__main__':
    alias=module('verified_study_alias',Path('scripts/with_study_output_path_alias.py'))
    with alias.base.validated_alias(Path('archive/output-relocation/2026-10-10-output/studies-manifest.json'),
            '8fb2545d8bae35a1a673ec3a43d3cc191148f03cac58892f29bd368182a8aff6',
            source_root=alias.SOURCE_ROOT,target_root=alias.TARGET_ROOT):main()
