"""CPU-only read-only completed-study provenance and exact-start audit."""
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback

sys.dont_write_bytecode=True
ROOT=Path('/home/sji/Documents/Programming/Python/embryogenesis-simulation')
FOLDER=ROOT/'outputs/phase-carry-redistribution-response'
RECORD=ROOT/'archive/output-relocation/2026-10-07/manifest.json'
RECORD_SHA='edcf6b433ed7d3db70adb89c8ca701ec49b0d73a5b981380c61a966bf89a3c02'
DEST=Path('/mnt/data/Documents/Python/embryogenesis-simulation/outputs')
REPORT=Path('/tmp/redistribution-response-provenance-verification.json')
sys.path.insert(0,str(ROOT))
import numpy as np
import torch

read=lambda p:json.loads(Path(p).read_text())

def metadata(p):
 s=Path(p).stat();return dict(device=s.st_dev,inode=s.st_ino,bytes=s.st_size,mtime_ns=s.st_mtime_ns,mode=s.st_mode,nlink=s.st_nlink)

def direct_sha(p):
 before=metadata(p)
 with Path(p).open('rb') as stream:digest=hashlib.file_digest(stream,'sha256').hexdigest()
 assert before==metadata(p),p
 return digest

assert direct_sha(RECORD)==RECORD_SHA
spec=importlib.util.spec_from_file_location('approved_output_alias',ROOT/'scripts/with_output_path_alias.py')
alias=importlib.util.module_from_spec(spec);spec.loader.exec_module(alias)
assert (ROOT/'outputs').is_symlink() and os.readlink(ROOT/'outputs')==str(DEST)
archives={}
archive_pins={}
for date in ('2026-10-06','2026-10-07-redistribution-launch'):
 file=ROOT/'archive/output-compression'/date/'manifest.json';archive_pins[str(file)]=direct_sha(file)
 for row in read(file)['files']:archives.setdefault(row['path'],[]).append((file,row))
cache={};evidence={};initial_stat={};archived=[]

def sha(path,expected=None):
 path=Path(path);name=str(path)
 if not path.is_absolute():path=ROOT/path;name=str(path)
 if not path.exists():
  relative=path.relative_to(ROOT).as_posix()
  rows=[r for _,r in archives.get(relative,[]) if expected is None or r['sha256']==expected]
  assert len(rows)==1,(relative,'archive mapping ambiguous')
  row=rows[0];assert row.get('codec','zstd')=='zstd',relative
  blob=ROOT/row['blob'];assert direct_sha(blob)==row['compressed_sha256']
  digest=hashlib.sha256();count=0
  with subprocess.Popen(['zstd','-q','-d','-c',str(blob)],stdout=subprocess.PIPE) as process:
   while block:=process.stdout.read(4<<20):digest.update(block);count+=len(block)
   assert process.wait()==0
  assert count==row['bytes'] and digest.hexdigest()==row['sha256']
  actual=digest.hexdigest();archived.append(dict(path=name,blob=row['blob'],sha256=actual))
 else:
  current=metadata(path)
  if name in cache:
   assert current==initial_stat[name],path
   actual=cache[name]
  else:
   actual=direct_sha(path);cache[name]=actual;initial_stat[name]=current
 if expected is not None:assert actual==expected,(name,expected,actual)
 evidence[name]=actual
 return actual

def state(p):
 with np.load(p,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}

report=dict(passed=False,started_utc=datetime.now(timezone.utc).isoformat(),scope='Independent exact preparations, frozen science/runtime/input hashes, completed checkpoint/history and stored gate provenance; no scientific outcome recomputation.',gpu_disabled=True,original_files_written=0,permanent_restoration=False)
started=time.monotonic()
try:
 assert not torch.cuda.is_available(),'GPU must be disabled for this audit'
 torch.set_num_threads(1)
 p=read(FOLDER/'protocol.json');ph=sha(FOLDER/'protocol.json')
 status=read(FOLDER/'status.json');summary=read(FOLDER/'summary.json')
 sha(FOLDER/'status.json');sha(FOLDER/'summary.json')
 assert status==dict(state='completed',passed=True,completed=66,total=66,new_jobs=48,reused_jobs=18,failed_histories=[],independent_histories=3,new_histories=0)
 assert summary['passed'] and summary['completed']==summary['total']==66 and summary['missing']==[] and summary['protocol_sha256']==ph
 assert len(p['jobs'])==66 and len(p['reused'])==18 and len(p['contexts'])==9 and p['histories']==[7,8,9]
 assert not any('/manuscript/' in f for key in ('source_sha256','input_sha256') for f in p[key])
 for key in ('source_sha256','input_sha256'):
  for index,(file,expected) in enumerate(p[key].items(),1):
   sha(file,expected)
   if index%300==0:print('Verified',key,index,'/',len(p[key]),flush=True)
  print('Verified',key,len(p[key]),flush=True)
 for key in ('experimental_binary','accepted_binary'):sha(p['device'][key],p['device'][key+'_sha256'])
 assert torch.__version__==p['device']['torch'] and torch.version.cuda==p['device']['torch_cuda']
 launch=read(FOLDER/'launch.json');sha(FOLDER/'launch.json')
 for file,key in (('preparation-verification.json','preparation_verification_sha256'),('preflight.json','preflight_sha256'),('implementation-verification.json','implementation_verification_sha256')):
  sha(FOLDER/file,launch[key])
 assert launch['protocol_sha256']==ph and launch['total_paths']==66 and launch['new_paths']==48 and launch['reused_paths']==18
 prep=read(FOLDER/'preparation-verification.json');saved_preparations={r['job']:r for r in prep['paths']}
 assert prep['passed'] and prep['protocol_sha256']==ph and prep['prepared_paths_verified']==66
 independent=[];errors=[]
 for seed in (7,8,9):
  contexts=[c for c in p['contexts'] if c['seed']==seed]
  originals={c['name']:state(c['source']) for c in contexts}
  masses={c['name']:np.asarray(p['preparations'][f'{seed}_{c["name"]}']['measured_gpu_masses']) for c in contexts}
  for ctx in contexts:
   original=originals[ctx['name']];history=read(Path(ctx['source']).parent/'history.json')
   assert history[-1]['time']==510. and np.array_equal(history[-1]['ids'],original['ids'])
   assert np.array_equal(history[-1]['chemistry'],[original['activator'],original['inhibitor']])
   assert np.array_equal(history[-1]['volumes'],masses[ctx['name']])
   assert json.loads(str(original['metadata']))['time']==510.
  common={}
  for target in p['selections'][str(seed)]['selected_ids']:
   amounts=[]
   for background in ('sham',f'redistributed-cell-{target}'):
    original=originals[background];index=original['ids'].tolist().index(target)
    amounts.append(float(original['activator'][index]*masses[background][index]))
   common[target]=.1*min(amounts)
   assert p['common_doses'][f'{seed}_{target}']==dict(background_target_amounts=amounts,removed_amount=common[target])
  for job in [j for j in p['jobs'] if j['seed']==seed]:
   original=originals[job['background']];actual=state(job['source'])
   chemistry=np.array([original['activator'],original['inhibitor']]);amount=0.;before=None
   if job['target'] is not None:
    index=original['ids'].tolist().index(job['target']);before=float(chemistry[0,index]*masses[job['background']][index])
    requested=.1*before if job['dose']=='fractional' else common[job['target']]
    factor=.9 if job['dose']=='fractional' else 1-requested/before
    assert 0<factor<1 and job['factor']==factor and job['requested_removed_amount']==requested
    chemistry[0,index]*=factor
    amount=float((original['activator'][index]-chemistry[0,index])*masses[job['background']][index]);errors.append(abs(amount/requested-1))
    assert amount==job['removed_activator_amount'] and abs(np.log(factor))==job['initial_log_amplitude']
   else:assert job['factor']==1. and job['initial_log_amplitude']==0. and job['removed_activator_amount']==0.
   meta=json.loads(str(original['metadata']));old_dt=meta['config']['dt'];meta['step_number']=round(meta['time']/job['dt'])
   for key in ('steps','save_every'):
    value=meta['config'][key]*old_dt/job['dt'];assert abs(value-round(value))<1e-8;meta['config'][key]=round(value)
   meta['config']['dt']=job['dt'];assert meta['config']==p['accepted_configs'][f'{seed}_{job["level"]}']
   assert actual.keys()=={k for k in original if k!='long_experiment'}
   for key in actual:
    expected=np.array(json.dumps(meta)) if key=='metadata' else chemistry[0] if key=='activator' else chemistry[1] if key=='inhibitor' else original[key]
    assert actual[key].dtype==expected.dtype and np.array_equal(actual[key],expected),(job['key'],key)
   entry=dict(job=job['key'],reused=job['key'] in p['reused'],source_sha256=sha(job['source']),removed_activator_amount=amount,initial_target_amount=before,initial_log_amplitude=job['initial_log_amplitude'],all_fields_exact=True)
   assert entry==saved_preparations[job['key']],job['key']
   independent.append(entry)
  print('Independently verified all preparations for history',seed,flush=True)
 assert max(errors)<=2e-13 and max(errors)==prep['maximum_relative_pulse_amount_error']
 with alias.validated_alias(RECORD,RECORD_SHA):
  from embryo import phase_carry_redistribution_response as study
  from embryo import phase_carry_exchange_response as moving
  from embryo.validate_gpu_backend import CRITERIA
  logical=[];gates=[]
  for index,job in enumerate(p['jobs'],1):
   result,history=study.current_evidence(FOLDER,job,p)
   assert next(r for r in summary['results'] if r['job']['key']==job['key'])==dict(job=job,reused=job['key'] in p['reused'],evidence=result)
   if job['key'] in p['reused']:
    reused=p['reused'][job['key']];oldroot=Path(p['reuse_parent']);oldp=read(oldroot/'protocol.json');oldjob=reused['job']
    sha(oldroot/'protocol.json',reused['protocol_sha256'])
    for file,expected in reused['evidence_sha256'].items():sha(file,expected)
    gate=moving.context_gate(oldroot,oldjob,oldp);folder=oldroot/oldjob['key']
   else:
    gate=moving.context_gate(FOLDER,job,p);folder=FOLDER/job['key']
   sha(folder/'prefix/result.json')
   assert gate['passed'] and gate['independent_chemical_steps']==8 and gate['restart_exact_steps']==4
   assert gate['cpu_chemical_log_max']<=1e-11 and gate['amount_accounting_max']<=p['amount_error_max']
   assert all(v<=CRITERIA[k] for k,v in gate['errors'].items()) and gate['final_phi_abs_max']<=CRITERIA['final_phi_abs_max']
   for file,expected in gate['evidence_sha256'].items():sha(file,expected)
   for file in ('result.json','history.json','latest_state.npz'):sha(folder/file)
   logical.append(dict(job=job['key'],reused=job['key'] in p['reused'],history_rows=len(history),final_time=history[-1]['time'],final_elapsed=history[-1]['elapsed'],result_sha256=sha(folder/'result.json'),quality_pass=result['quality_pass']))
   gates.append(dict(job=job['key'],reused=job['key'] in p['reused'],prefix_sha256=sha(folder/'prefix/result.json'),protocol_sha256=gate['protocol_sha256'],errors=gate['errors'],final_phi_abs_max=gate['final_phi_abs_max'],cpu_chemical_log_max=gate['cpu_chemical_log_max'],amount_accounting_max=gate['amount_accounting_max']))
   if index%11==0:print('Verified completed paths/context gates',index,'/66',flush=True)
  pairs=[];pilot=[];superseded=[]
  for seed in (7,8,9):
   by_history=next(row for row in summary['by_history'] if row['seed']==seed);assert by_history['qualified']
   for label,horizon in (('pilot',6.),('full',60.)):
    path=FOLDER/'pairs'/f'seed-{seed}'/f'{label}-refinement.json';saved=read(path);sha(path)
    assert saved==by_history['reports'][label] and saved['passed'] and saved['protocol_sha256']==ph and saved['horizon']==horizon
    assert len(saved['mechanical'])==11 and all(row['passed'] for row in saved['mechanical']) and saved['responses']['passed']
    for file,expected in saved['evidence_sha256'].items():sha(file,expected)
    pairs.append(dict(seed=seed,label=label,sha256=sha(path),mechanical_pairs=len(saved['mechanical']),responses_passed=True))
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
     superseded.append(dict(job=job['key'],pilot_checkpoint_sha256=next(iter(row['evidence_sha256'].values())),interpretation='Original pilot checkpoint is deliberately overwritten by qualified full continuation; immutable pilot history/field/quality records verify the gate. Full cumulative physical audit bounds pilot.'))
   pilot.append(dict(seed=seed,rows=22,passed=True,sha256=sha(path)))
  pre=read(FOLDER/'preflight.json');assert pre['passed'] and pre['protocol_sha256']==ph and pre['full_size_contexts']==4
  for file,expected in pre['evidence_sha256'].items():sha(file,expected)
  for gate in pre['contexts']:
   job=gate['job'];assert gate==read(FOLDER/job['key']/'prefix/result.json')
 impl=read(FOLDER/'implementation-verification.json');assert impl['passed'] and impl['total_regression_cases_passed']==110
 for file,expected in impl['new_files_sha256'].items():sha(ROOT/file,expected)
 sha(FOLDER/'regression-tests.log',impl['log_sha256'])
 snapshot=read(ROOT/'archive/study-document-snapshots/2026-10-07-before-redistribution-response-launch/manifest.json')
 historical=[]
 for file,expected in snapshot['unchanged_scientific_documents'].items():
  actual=direct_sha(ROOT/file)
  if actual!=expected:historical.append(dict(path=file,historical_sha256=expected,current_sha256=actual,scope='Historical launch document snapshot; authorized manuscript edit is not a scientific protocol hash mismatch.'))
 assert all(row['path'].startswith('manuscript/') for row in historical),historical
 for file,expected in archive_pins.items():assert direct_sha(file)==expected
 assert direct_sha(RECORD)==RECORD_SHA
 print('Rechecking immutable file metadata and hashes',len(cache),flush=True)
 for index,(file,expected) in enumerate(cache.items(),1):
  assert metadata(file)==initial_stat[file],file
  assert direct_sha(file)==expected,file
  if index%400==0:print('Immutable recheck',index,'/',len(cache),flush=True)
 report.update(passed=True,protocol_sha256=ph,scientific_source_hashes=len(p['source_sha256']),input_hashes=len(p['input_sha256']),runtime_binary_hashes=2,torch_runtime=p['device']['torch'],torch_cuda_runtime=p['device']['torch_cuda'],manuscript_pins_in_protocol=0,prepared_paths_verified=66,new_paths=48,reused_paths=18,independent_histories=3,new_histories=0,maximum_relative_pulse_amount_error=max(errors),preparations=independent,completed_paths=logical,context_gates=gates,paired_reports=pairs,pilot_quality=pilot,pilot_checkpoints_intentionally_superseded=superseded,historical_document_snapshot_guards=len(snapshot['unchanged_scientific_documents']),historical_document_differences=historical,historical_differences_are_scientific_failures=False,evidence_sha256=evidence,evidence_files_rechecked=len(cache),archive_streams=archived,approved_outputs_alias=str(DEST),migration_manifest_sha256=RECORD_SHA,all_mismatches=[])
except BaseException as exc:
 report.update(error=str(exc),traceback=traceback.format_exc(),all_mismatches=[str(exc)])
 raise
finally:
 report['finished_utc']=datetime.now(timezone.utc).isoformat();report['elapsed_seconds']=time.monotonic()-started
 REPORT.write_text(json.dumps(report,indent=2)+'\n')
 print('Read-only provenance passed',report['passed'],'report',REPORT,flush=True)
