"""Durable, bounded CPU-only stages for completed-study provenance."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import tempfile
import textwrap
import time

BASE=Path('/tmp/verify_redistribution_response_provenance.py')
original=BASE.read_text()
exec(compile(original[:original.index('report=dict')],str(BASE),'exec'),globals())
PROGRESS=Path('/tmp/redistribution-response-provenance-progress.json')
REPORT=Path('/tmp/redistribution-response-provenance-verification.json')


def save(progress):
 progress['updated_utc']=datetime.now(timezone.utc).isoformat()
 progress['cache']=cache;progress['initial_stat']=initial_stat;progress['evidence_sha256']=evidence;progress['archive_streams']=archived
 with tempfile.NamedTemporaryFile(mode='w',dir=PROGRESS.parent,delete=False) as stream:
  path=Path(stream.name);json.dump(progress,stream,indent=2);stream.write('\n');stream.flush();os.fsync(stream.fileno())
 os.replace(path,PROGRESS)


def load():
 progress=read(PROGRESS);cache.update(progress['cache']);initial_stat.update(progress['initial_stat']);evidence.update(progress['evidence_sha256']);archived.extend(progress['archive_streams']);return progress


def main():
 parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['initial','paths-7','paths-8','paths-9','reports','final']);parser.add_argument('--chunk',type=int,default=400);args=parser.parse_args()
 assert not torch.cuda.is_available();torch.set_num_threads(1)
 if args.stage=='initial':
  assert not PROGRESS.exists(),'Initial progress already exists'
  body=original[original.index(' assert not torch.cuda.is_available()'):original.index(' with alias.validated_alias(RECORD,RECORD_SHA):')]
  exec(compile(textwrap.dedent(body),str(BASE),'exec'),globals())
  p=globals()['p'];ph=globals()['ph'];independent=globals()['independent'];errors=globals()['errors']
  progress=dict(passed=False,stage='initial_complete',started_utc=datetime.now(timezone.utc).isoformat(),protocol_sha256=ph,source_hashes_verified=len(p['source_sha256']),input_hashes_verified=len(p['input_sha256']),runtime_binary_hashes=2,torch_runtime=p['device']['torch'],torch_cuda_runtime=p['device']['torch_cuda'],manuscript_pins_in_protocol=0,prepared_paths=independent,prepared_paths_verified=66,maximum_relative_pulse_amount_error=max(errors),new_paths=48,reused_paths=18,independent_histories=3,new_histories=0,completed_paths=[],context_gates=[],paired_reports=[],pilot_quality=[],history_stages_completed=[],all_mismatches=[],original_files_written=0,permanent_restoration=False,scope='Frozen science/input/runtime hashes, independently reconstructed preparations, exact reuse, completed-path and stored-prefix/pilot/full provenance; root separately recomputes physical/response metrics.',initial_script=str(BASE),verification_script=str(Path(__file__)),approved_outputs_alias=str(DEST),migration_manifest_sha256=RECORD_SHA)
  save(progress);print('Initial verification durably saved',flush=True);return
 progress=load();p=read(FOLDER/'protocol.json');ph=sha(FOLDER/'protocol.json');assert ph==progress['protocol_sha256'];summary=read(FOLDER/'summary.json');sha(FOLDER/'summary.json');sha(FOLDER/'status.json')
 if args.stage.startswith('paths-'):
  seed=int(args.stage.split('-')[-1]);assert seed not in progress['history_stages_completed']
  with alias.validated_alias(RECORD,RECORD_SHA):
   from embryo import phase_carry_redistribution_response as study
   from embryo import phase_carry_exchange_response as moving
   from embryo.validate_gpu_backend import CRITERIA
   for index,job in enumerate([j for j in p['jobs'] if j['seed']==seed],1):
    reused=job['key'] in p['reused']
    if reused:
     item=p['reused'][job['key']];oldroot=Path(p['reuse_parent']);oldp=read(oldroot/'protocol.json');oldjob=item['job'];folder=oldroot/oldjob['key']
     sha(oldroot/'protocol.json',item['protocol_sha256'])
     for file,expected in item['evidence_sha256'].items():sha(file,expected)
     result,history=study.current_evidence(FOLDER,job,p)
     assert (folder/'prefix/result.json').exists()
     gate=moving.context_gate(oldroot,oldjob,oldp)
    else:
     folder=FOLDER/job['key'];result=read(folder/'result.json');history=read(folder/'history.json')
     assert result['job']==job and result['protocol_sha256']==ph and result['quality_pass']
     sha(job['source'],result['source_sha256']);sha(folder/'history.json',result['history_sha256']);sha(folder/'latest_state.npz',result['checkpoint_sha256'])
     with np.load(folder/'latest_state.npz',allow_pickle=False) as z:
      meta=json.loads(str(z['metadata']));saved=json.loads(str(z['long_experiment']));ids=z['ids'].tolist();chemistry=np.array([z['activator'],z['inhibitor']])
      assert str(z['precision_arm'])=='phase_carry' and z['rounding'].dtype==np.int32
     assert saved['job']==job and saved['protocol_hash']==ph and saved['history']==history and saved['audit']==result['audit']
     assert meta['config']==p['accepted_configs'][f'{seed}_{job["level"]}'] and abs(meta['time']-570.)<1e-9
     assert np.array_equal(history[-1]['chemistry'],chemistry) and history[-1]['ids']==ids
     assert len(history)==401 and np.allclose([r['elapsed'] for r in history],np.arange(401)*.15,rtol=0,atol=1e-9)
     assert all(r['ids']==ids and abs(r['time']-510.-r['elapsed'])<1e-9 for r in history)
     with np.load(job['source'],allow_pickle=False) as z:assert np.array_equal(history[0]['chemistry'],[z['activator'],z['inhibitor']])
     assert (folder/'prefix/result.json').exists()
     gate=moving.context_gate(FOLDER,job,p)
    assert next(r for r in summary['results'] if r['job']['key']==job['key'])==dict(job=job,reused=reused,evidence=result)
    assert len(history)==401 and abs(history[-1]['time']-570.)<1e-9 and abs(history[-1]['elapsed']-60.)<1e-9
    assert gate['passed'] and gate['independent_chemical_steps']==8 and gate['restart_exact_steps']==4
    assert gate['cpu_chemical_log_max']<=1e-11 and gate['amount_accounting_max']<=p['amount_error_max']
    assert all(v<=CRITERIA[k] for k,v in gate['errors'].items()) and gate['final_phi_abs_max']<=CRITERIA['final_phi_abs_max']
    for file,expected in gate['evidence_sha256'].items():sha(file,expected)
    for file in ('result.json','history.json','latest_state.npz'):sha(folder/file)
    progress['completed_paths'].append(dict(job=job['key'],reused=reused,history_rows=len(history),final_time=history[-1]['time'],final_elapsed=history[-1]['elapsed'],result_sha256=sha(folder/'result.json'),quality_pass=result['quality_pass']))
    progress['context_gates'].append(dict(job=job['key'],reused=reused,prefix_sha256=sha(folder/'prefix/result.json'),protocol_sha256=gate['protocol_sha256'],errors=gate['errors'],final_phi_abs_max=gate['final_phi_abs_max'],cpu_chemical_log_max=gate['cpu_chemical_log_max'],amount_accounting_max=gate['amount_accounting_max']))
    if index%11==0:print('Completed path/prefix provenance',seed,index,'/22',flush=True)
  progress['history_stages_completed'].append(seed);progress['stage']=args.stage+'_complete';save(progress);print('History provenance durably saved',seed,flush=True);return
 if args.stage=='reports':
  assert sorted(progress['history_stages_completed'])==[7,8,9]
  with alias.validated_alias(RECORD,RECORD_SHA):
   from embryo import phase_carry_exchange_response as moving
   superseded=[]
   for seed in (7,8,9):
    by_history=next(row for row in summary['by_history'] if row['seed']==seed);assert by_history['qualified']
    for label,horizon in (('pilot',6.),('full',60.)):
     path=FOLDER/'pairs'/f'seed-{seed}'/f'{label}-refinement.json';saved=read(path);sha(path)
     assert saved==by_history['reports'][label] and saved['passed'] and saved['protocol_sha256']==ph and saved['horizon']==horizon
     assert len(saved['mechanical'])==11 and all(row['passed'] for row in saved['mechanical']) and saved['responses']['passed']
     for file,expected in saved['evidence_sha256'].items():sha(file,expected)
     progress['paired_reports'].append(dict(seed=seed,label=label,sha256=sha(path),mechanical_pairs=11,responses_passed=True))
    path=FOLDER/'pairs'/f'seed-{seed}'/'pilot-quality.json';quality=read(path);sha(path)
    assert quality['passed'] and quality['protocol_sha256']==ph and len(quality['rows'])==22
    for row in quality['rows']:
     assert row['passed'] and moving.quality_pass(row['audit'],p)
     job=next(j for j in p['jobs'] if j['key']==row['job']);current=next(r['evidence']['audit'] for r in summary['results'] if r['job']['key']==job['key'])
     assert current['min_radius']<=row['audit']['min_radius']+1e-12
     assert all(current[k]+1e-15>=row['audit'][k] for k in ('max_volume_error','max_clipping','boundary_max','dilution_error_max','wall_seconds'))
     if job['key'] in p['reused']:
      for file,expected in row['evidence_sha256'].items():sha(file,expected)
     else:
      assert len(row['evidence_sha256'])==1
      superseded.append(dict(job=job['key'],pilot_checkpoint_sha256=next(iter(row['evidence_sha256'].values())),interpretation='Pilot checkpoint intentionally replaced by full continuation; immutable pilot history/field/quality records retained, and full cumulative audit bounds pilot.'))
    progress['pilot_quality'].append(dict(seed=seed,rows=22,passed=True,sha256=sha(path)))
   pre=read(FOLDER/'preflight.json');assert pre['passed'] and pre['protocol_sha256']==ph and pre['full_size_contexts']==4
   for file,expected in pre['evidence_sha256'].items():sha(file,expected)
   for gate in pre['contexts']:assert gate==read(FOLDER/gate['job']['key']/'prefix/result.json')
  impl=read(FOLDER/'implementation-verification.json');assert impl['passed'] and impl['total_regression_cases_passed']==110
  for file,expected in impl['new_files_sha256'].items():sha(ROOT/file,expected)
  sha(FOLDER/'regression-tests.log',impl['log_sha256'])
  snapshot=read(ROOT/'archive/study-document-snapshots/2026-10-07-before-redistribution-response-launch/manifest.json');historical=[]
  for file,expected in snapshot['unchanged_scientific_documents'].items():
   actual=direct_sha(ROOT/file)
   if actual!=expected:historical.append(dict(path=file,historical_sha256=expected,current_sha256=actual,scope='Historical launch document snapshot; user-authorized manuscript edit is outside current science/protocol hash pins.'))
  assert all(row['path'].startswith('manuscript/') for row in historical),historical
  progress.update(pilot_checkpoints_intentionally_superseded=superseded,historical_document_snapshot_guards=len(snapshot['unchanged_scientific_documents']),guard_mismatch_summary=dict(frozen_science_source_input_runtime_mismatches=[],historical_document_differences=historical,historical_differences_are_scientific_failures=False),archive_resolution_summary=dict(streams=archived,permanent_restoration=False),stage='reports_complete',final_rechecked=[])
  save(progress);print('All report provenance durably saved',flush=True);return
 assert progress['stage'] in ('reports_complete','final_in_progress')
 done=set(progress['final_rechecked']);remaining=[name for name in cache if name not in done]
 for index,name in enumerate(remaining[:args.chunk],1):
  assert metadata(name)==initial_stat[name],name
  assert direct_sha(name)==cache[name],name
  progress['final_rechecked'].append(name)
  if index%100==0:print('Post-audit immutable recheck',len(progress['final_rechecked']),'/',len(cache),flush=True)
 progress['stage']='final_in_progress'
 for file,expected in archive_pins.items():assert direct_sha(file)==expected
 assert direct_sha(RECORD)==RECORD_SHA
 if len(progress['final_rechecked'])==len(cache):
  assert len(progress['completed_paths'])==len(progress['context_gates'])==len(progress['prepared_paths'])==66
  assert len(progress['paired_reports'])==6 and len(progress['pilot_quality'])==3
  progress.update(passed=True,stage='completed',finished_utc=datetime.now(timezone.utc).isoformat(),immutable_before_after_registry={name:dict(sha256=cache[name],before=initial_stat[name],after=metadata(name)) for name in cache},immutable_files_rechecked=len(cache),verification_script_sha256=direct_sha(Path(__file__)),initial_script_sha256=direct_sha(BASE))
 save(progress)
 if progress['passed']:
  final={k:v for k,v in progress.items() if k not in ('cache','initial_stat','final_rechecked')};REPORT.write_text(json.dumps(final,indent=2)+'\n');print('Final provenance passed',REPORT,flush=True)
 else:print('Post-audit progress durably saved',len(progress['final_rechecked']),'/',len(cache),flush=True)

if __name__=='__main__':main()
