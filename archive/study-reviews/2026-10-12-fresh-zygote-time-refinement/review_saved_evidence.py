"""Independent saved-array review of history-1103's finer developmental pair.

NumPy/SciPy only: no embryo solver/comparison imports and no GPU stepping.
Preserve every original failed report and all production records.
"""
from datetime import datetime, timezone
from pathlib import Path
import importlib.util
import json
import time

import numpy as np

ROOT = Path('output/studies/fresh-zygote-time-refinement')
PARENT = Path('output/studies/fresh-zygote-pilot')
HERE = Path(__file__).resolve().parent
EXPECTED = '54632a1f284affcbcec93cd0f1036d2d9c3788bd93ef5480df5fa1f14336c6c3'
HELPER = Path('archive/study-reviews/2026-10-11-development-controls/review_saved_arrays.py')
PRIOR = Path('archive/study-reviews/2026-10-11-fresh-zygote-stopped/review_saved_evidence.py')


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


c = module('independent_development_controls', HELPER)
c.ROOT = ROOT
s, t, r = c.s, c.t, c.r


def main():
    started = time.perf_counter()
    output = HERE/'independent_review.json'
    assert not output.exists(), 'Never overwrite a recorded review'
    for name in (Path(__file__),HELPER,c.HELPER,s.HELPER,t.HELPER,PRIOR):r.pin(name)
    p = r.read(ROOT/'protocol.json'); ph = r.pin(ROOT/'protocol.json',EXPECTED)
    status = r.read(ROOT/'status.json'); summary = r.read(ROOT/'summary.json')
    original = r.read(ROOT/'time-h1103-off.json')
    assert status['state']=='completed_provisional' and status['completed']==status['total']==1
    assert summary['execution_completed'] and summary['developmental_screens_passed']
    assert not summary['qualification_accepted'] and summary['independent_review_pending']
    assert status['protocol_sha256']==summary['protocol_sha256']==original['protocol_sha256']==ph
    assert summary['new_scientific_histories']==0 and summary['new_numerical_paths']==1
    assert summary['reused_numerical_paths']==1 and not summary['post_cleavage_started']
    assert summary['pair']==original and original['passed']
    for name,expected in {**p['source_sha256'],**p['input_sha256']}.items():r.pin(name,expected)
    for name,expected in original['evidence_sha256'].items():r.pin(name,expected)
    runtime = r.read(ROOT/'runtime.json');r.pin(runtime['binary'],runtime['binary_sha256'])
    assert runtime['binary_sha256']==p['runtime_binary_sha256']
    parent = r.read(PARENT/'protocol.json')
    assert r.pin(PARENT/'protocol.json',p['parent_protocol_sha256'])=='728f584676ea4eeeee530656015aec61931e9296ff095eff8b7c94f9832a13ee'
    r.pin(PARENT/'status.json',p['original_failed_status_sha256'])
    r.pin(PARENT/'time-h1103-off.json',p['original_failed_pair_sha256'])
    old_status=r.read(PARENT/'status.json');old_pair=r.read(PARENT/'time-h1103-off.json')
    assert old_status['state']=='failed' and old_status['completed']==18 and not old_pair['passed']
    assert old_pair['dense']['peak']['relative_spectrum']>.01
    inherited=('duration','interval','checkpoint_interval','diagnostic_interval','snapshots',
        'time_limits','physical_limits','restart_steps','restart_stages','minimum_mature_window',
        'method','gpu_uuid','dense_time_spectrum_limit','instrumentation_gate',
        'instrumentation_protocol_sha256','runtime_binary_sha256','prefix_steps','prefix_limits')
    for key in inherited:assert p[key]==parent[key],key
    assert p['time_limits']['field_max']==p['time_limits']['relative_spectrum_max']==p['dense_time_spectrum_limit']==.01
    paths=[];reports={};histories={};networks={};jobs={}
    for job in p['jobs']:
        check,report,history,network=c.review_path(job,p,ph)
        paths.append(check);reports[job['key']]=report;histories[job['key']]=history
        networks[job['key']]=network;jobs[job['key']]=job
    assert len(paths)==2 and sum(j['reused'] for j in p['jobs'])==1
    a,b=original['cases'];assert [a,b]==[j['key'] for j in p['jobs']]
    assert {k for k,v in jobs[b]['config'].items() if v!=jobs[a]['config'][k]}=={'dt','steps','save_every'}
    assert jobs[a]['config']['dt']==.000234375 and jobs[b]['config']['dt']==.0001171875
    # Exact matching of all arrays and non-config metadata, including random streams.
    x,y=[r.load(jobs[key]['source']) for key in (a,b)]
    assert x.keys()==y.keys()
    mx,my=[json.loads(str(z['metadata'])) for z in (x,y)]
    assert mx.pop('config')==jobs[a]['config'] and my.pop('config')==jobs[b]['config']
    assert mx==my and my['time']==0 and my['step_number']==0 and not my['divisions']
    for name in set(x)-{'metadata'}:np.testing.assert_array_equal(x[name],y[name],err_msg=name)
    del x,y
    prefix=r.read(ROOT/b/'actual-start-prefix.json')
    assert prefix==summary['prefix'] and prefix['passed'] and prefix['failure'] is None
    assert prefix['protocol_sha256']==ph and prefix['steps']==8 and prefix['exact_gpu_restart']
    for key,value in prefix['errors'].items():assert value<=p['prefix_limits'][key]
    for name,expected in prefix['evidence_sha256'].items():r.pin(name,expected)
    row=s.compare_histories(original['temporal'],histories,reports,networks)
    dense=row['dense'];saved=original['dense']
    r.close(dense['peak']['relative_spectrum'],saved['peak']['relative_spectrum'],'dense spectrum peak')
    assert dense['peak']['time']==saved['peak']['time'] and dense['peak']['edges']==saved['peak']['edges']
    assert dense['unmatched_times']==saved['unmatched_times']
    assert dense['topology_mismatch_times']==saved['topology_mismatch_times']
    checkpoints=[];live=[];endpoints=[]
    for stamp in (0,3,6,9,12,15,18):
        fields={};ranges={}
        for key in (a,b):
            z=r.load(ROOT/key/f'gpu-t{stamp:02d}.npz')
            check=r.checkpoint(z,jobs[key]['config'],histories[key][round(stamp/.15)],stamp==18)
            assert check['volume_error']<p['physical_limits']['volume_error_max']
            assert check['minimum_radius']>=p['physical_limits']['radius_min']
            assert json.loads(str(z['event_audits']))==[e for e in reports[key]['events'] if e['time']<=stamp+1e-12]
            checkpoints.append(dict(case=key,**check,corrections_present=True))
            fields[key],ranges[key]=s.density(z,None)
            if stamp==18:
                h=histories[key][-1];lam,_=r.eig(h);cfg=jobs[key]['config']
                jac=np.array([[1.,-1.],[2*cfg['signal_beta'],-cfg['signal_beta']]])
                growth=np.array([np.linalg.eigvals(jac-u*np.diag([cfg['signal_da'],cfg['signal_dh']])).real.max() for u in lam])
                endpoints.append(dict(case=key,dt=cfg['dt'],history=1103,branch='off',time=18,
                    cells=len(h['ids']),axis_ratio=h['axis_ratio'],log_activator_sd=h['log_activator_sd'],
                    chemistry_ranges=[[min(u),max(u)] for u in h['chemistry']],
                    frozen_uniform_unstable_modes=int(np.sum((lam>1e-10)&(growth>1e-10))),
                    max_frozen_uniform_growth=float(growth[lam>1e-10].max()),active_furrows=len(h['dividing_ids'])))
            del z
        errors=r.field_error(fields[a],fields[b])
        expected=next(item for item in original['temporal']['fields'] if item['time']==stamp)
        for key,value in errors.items():r.close(value,expected['errors'][key],('field',stamp,key))
        np.testing.assert_allclose([ranges[a],ranges[b]],expected['phase_ranges'],rtol=1e-6,atol=1e-9)
        row['field_checks'].append(dict(time=stamp,probe=None,order=3,errors=errors))
        del fields
        print('Independently rebuilt both checkpoint fields:',stamp,flush=True)
    for key in (a,b):
        for stage in (2,4,8):
            z=r.load(ROOT/key/f'restart-stage-{stage}.npz')
            check=r.checkpoint(z,jobs[key]['config'])
            assert stage<check['cells']<2*stage and check['active']>0
            live.append(dict(case=key,stage=stage,**check));del z
    s.decision(row,original['temporal'],p['time_limits'])
    assert row['independently_passed'] and dense['peak']['relative_spectrum']<=p['dense_time_spectrum_limit']
    # Recompute the new three-level graph trend, retaining unequal-count clocks.
    coarser='h1103_off_development'
    old_networks={}
    for key in (coarser,a):
        records=t.jsonl(PARENT/key/'network.jsonl')
        old_networks[key]={round(x['time'],9):x for x in records if abs(x['time']/.015-round(x['time']/.015))<1e-7}
    old_dense=s.dense_compare(old_networks[coarser],old_networks[a])
    r.close(old_dense['peak']['relative_spectrum'],old_pair['dense']['peak']['relative_spectrum'],'retained failure')
    assert old_dense['unmatched_times']==old_pair['dense']['unmatched_times']
    spectral_series=[]
    for stamp in sorted(networks[a]):
        old_x,old_y=old_networks[coarser][stamp],old_networks[a][stamp]
        new_x,new_y=networks[a][stamp],networks[b][stamp]
        values=dict(time=stamp,old_relative_spectrum=None,new_relative_spectrum=None,
            old_edges=[old_x['edges'],old_y['edges']],new_edges=[new_x['edges'],new_y['edges']])
        for label,x,y in [('old',old_x,old_y),('new',new_x,new_y)]:
            if len(x['ids'])==len(y['ids']):
                ex,ey=r.eig(x)[0],r.eig(y)[0]
                values[label+'_relative_spectrum']=float(np.linalg.norm(ex-ey)/(np.linalg.norm(ey) or 1.))
        spectral_series.append(values)
    former=next(x for x in spectral_series if x['time']==11.43)
    assert former['new_relative_spectrum'] is not None and former['new_relative_spectrum']<.01
    peak_time=dense['peak']['time'];nx,ny=networks[a][peak_time],networks[b][peak_time]
    peak_edges=[]
    if nx['ids']==ny['ids'] and nx['parents']==ny['parents']:
        ga,gb=np.array(nx['conductance']),np.array(ny['conductance'])
        for i,j in zip(*np.triu_indices(len(ga),1)):
            if (ga[i,j]>0)!=(gb[i,j]>0):peak_edges.append(dict(cells=[nx['ids'][i],nx['ids'][j]],
                margins=[nx['contact'][i][j]-nx['cutoff'],ny['contact'][i][j]-ny['cutoff']],
                conductance=[float(ga[i,j]),float(gb[i,j])]))
    for name,expected in list(r.WATCHED.items()):assert r.sha(name)==expected,name
    result=dict(review_completed=True,scoped_refinement_accepted=True,original_cohort_accepted=False,
        qualification_scope='History1103 feedback-off at grid120, dt0.000234375 versus dt0.0001171875, t0–18; inherited screens only.',
        reviewed_utc=datetime.now(timezone.utc).isoformat(),seconds=time.perf_counter()-started,
        protocol_sha256=ph,new_scientific_histories=0,new_numerical_paths=1,reused_numerical_paths=1,
        source_count=len(p['source_sha256']),input_count=len(p['input_sha256']),paths=paths,
        recorded_actual_start_prefix=prefix,pair=row,checkpoints=checkpoints,live_checkpoints=live,endpoints=endpoints,
        original_failure_retained=True,original_dense_failure=old_dense,
        prior_failure_time_comparison=former,new_peak_cutoff_edges=peak_edges,
        spectral_time_series=spectral_series,
        general_developmental_convergence_accepted=False,biological_identity_qualified=False,
        general_continuum_qualified=False,post_cleavage_started=False,
        scope='Independent NumPy/SciPy saved-evidence review: both histories, all dense graphs/pressure/events, all seven density checkpoints per path and all six live checkpoint states checked. Recorded actual-start and eight-step restart checks verified through saved evidence; no fresh GPU stepping, full fine CPU ladder or full-horizon restart replay. One existing history; remaining controls and long windows untested.',
        original_files_sha256=dict(r.WATCHED))
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('review_completed','scoped_refinement_accepted','original_cohort_accepted','seconds')},indent=2),flush=True)


if __name__=='__main__':
    # The established read-only migration alias retains canonical provenance names.
    alias=module('verified_study_path_alias',Path('scripts/with_study_output_path_alias.py'))
    with alias.base.validated_alias(Path('archive/output-relocation/2026-10-10-output/studies-manifest.json'),
            '8fb2545d8bae35a1a673ec3a43d3cc191148f03cac58892f29bd368182a8aff6',
            source_root=alias.SOURCE_ROOT,target_root=alias.TARGET_ROOT):
        main()
