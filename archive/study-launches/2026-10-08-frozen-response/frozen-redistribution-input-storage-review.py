"""Read-only nine-context frozen-assay input and storage review."""
from datetime import datetime,timezone
import hashlib,json,os,shutil,sys,time
from pathlib import Path
sys.dont_write_bytecode=True
import numpy as np
ROOT=Path('/home/sji/Documents/Programming/Python/embryogenesis-simulation')
FOLDER=ROOT/'outputs/phase-carry-redistribution-response'
APPROVED=Path('/mnt/data/Documents/Python/embryogenesis-simulation/outputs')
REPORT=Path('/tmp/frozen-redistribution-input-storage-review.json')
read=lambda path:json.loads(Path(path).read_text())
evidence={};stats={}

def metadata(path):
 s=Path(path).stat();return dict(device=s.st_dev,inode=s.st_ino,bytes=s.st_size,mode=s.st_mode,mtime_ns=s.st_mtime_ns,nlink=s.st_nlink)

def sha(path,expected=None):
 path=Path(path);before=metadata(path)
 with path.open('rb') as stream:value=hashlib.file_digest(stream,'sha256').hexdigest()
 assert before==metadata(path),path
 if expected is not None:assert value==expected,(str(path),expected,value)
 key=str(path);assert key not in evidence or evidence[key]==value,path;evidence[key]=value;stats[key]=before
 return value

def arraysha(array):
 value=np.ascontiguousarray(array);return dict(dtype=str(value.dtype),shape=list(value.shape),sha256=hashlib.sha256(value.tobytes()).hexdigest())

report=dict(passed=False,started_utc=datetime.now(timezone.utc).isoformat(),scope='Read-only exact t510 inputs and runtime/storage readiness for frozen chemical response; not a new scientific outcome.',files_written_under_outputs=0,permanent_decompression=False,phase_arrays_loaded_or_copied=False)
started=time.monotonic()
try:
 assert (ROOT/'outputs').is_symlink() and os.readlink(ROOT/'outputs')==str(APPROVED)
 migration=ROOT/'archive/output-relocation/2026-10-07/manifest.json';sha(migration,'edcf6b433ed7d3db70adb89c8ca701ec49b0d73a5b981380c61a966bf89a3c02')
 p=read(FOLDER/'protocol.json');ph=sha(FOLDER/'protocol.json');status=read(FOLDER/'status.json');summary=read(FOLDER/'summary.json')
 assert status['state']=='completed' and status['passed'] and status['completed']==status['total']==66
 assert summary['passed'] and summary['completed']==summary['total']==66 and summary['missing']==[] and summary['protocol_sha256']==ph
 compact_path=ROOT/'docs/phase_carry_redistribution_response_assessment.json';compact=read(compact_path);sha(compact_path)
 assert compact['full_study_accepted'] and compact['protocol_sha256']==ph
 sha(FOLDER/'status.json',compact['original_status_sha256']);sha(FOLDER/'summary.json',compact['original_summary_sha256'])
 assert compact['independent_histories']==3 and compact['new_histories']==0 and compact['paths']==66
 reviewroot=ROOT/compact['review_archive'];review_path=reviewroot/'review-provenance.json';review=read(review_path);sha(review_path)
 assert review['passed'] and review['protocol_sha256']==ph and review['original_summary_sha256']==compact['original_summary_sha256'] and review['original_status_sha256']==compact['original_status_sha256']
 excluded=[]
 for name,expected in review['files_sha256'].items():
  # These were captured as reporting state. A subsequent authorized launch can
  # update README/plan; they are not sources for the frozen chemical equations.
  if name in ('README.md','docs/plan.md'):excluded.append(name);continue
  sha(ROOT/name,expected)
 proofpath=reviewroot/'provenance_verification.json';proof=read(proofpath)
 assert proof['passed'] and proof['protocol_sha256']==ph and len(proof['prepared_paths'])==66 and len(proof['context_gates'])==66
 for name,expected in p['source_sha256'].items():sha(name,expected)
 for key in ('experimental_binary','accepted_binary'):sha(p['device'][key],p['device'][key+'_sha256'])
 controls=[j for j in p['jobs'] if j['level']=='fine' and j['target'] is None];assert len(controls)==9
 contexts=[]
 for j in controls:
  reused=j['key'] in p['reused']
  if reused:
   item=p['reused'][j['key']];rawjob=item['job'];rawroot=Path(p['reuse_parent']);sha(rawroot/'protocol.json',item['protocol_sha256']);rawph=item['protocol_sha256']
  else:rawjob=j;rawroot=FOLDER;rawph=ph
  folder=rawroot/rawjob['key'];source=Path(j['source']);historypath=folder/'history.json';resultpath=folder/'result.json'
  for file in (source,historypath,resultpath):sha(file,proof['evidence_sha256'][str(file)])
  for file in (source,):sha(file,p['input_sha256'][str(file)])
  result=read(resultpath);history=read(historypath);row=history[0]
  assert result['job']==rawjob and result['protocol_sha256']==rawph and result['quality_pass']
  sha(source,result['source_sha256']);sha(historypath,result['history_sha256'])
  assert len(history)==401 and row['elapsed']==0. and row['time']==510. and history[-1]['time']==570.
  with np.load(source,allow_pickle=False) as z:
   meta=json.loads(str(z['metadata']));ids=z['ids'].copy();chemistry=np.array([z['activator'],z['inhibitor']]);precision=str(z['precision_arm'])
  assert meta['time']==510. and meta['config']==p['accepted_configs'][f'{j["seed"]}_fine'] and precision=='phase_carry'
  assert np.array_equal(row['ids'],ids) and np.array_equal(row['chemistry'],chemistry)
  bg=p['preparations'][f'{j["seed"]}_{j["background"]}'];volumes=np.asarray(row['volumes'],dtype=np.float64);delta=np.asarray(row['delta'],dtype=np.float64)
  assert np.array_equal(bg['ids'],ids) and np.array_equal(bg['chemistry'],chemistry) and np.array_equal(bg['measured_gpu_masses'],volumes)
  assert chemistry.shape==(2,16) and delta.shape==(16,16) and volumes.shape==(16,) and np.all(chemistry>0) and np.all(volumes>0)
  assert all(np.isfinite(x).all() for x in (chemistry,volumes,delta))
  parent_context=next(c for c in p['contexts'] if c['seed']==j['seed'] and c['name']==j['background']);parent_original=Path(parent_context['source']);parent_history=parent_original.parent/'history.json';sha(parent_history,proof['evidence_sha256'][str(parent_history)])
  prior=read(parent_history)[-1];assert prior['time']==510.
  assert all(np.array_equal(prior[k],row[k]) for k in ('ids','chemistry','volumes','delta','polarity','centers'))
  coarse=next(c for c in p['jobs'] if c['seed']==j['seed'] and c['background']==j['background'] and c['level']=='coarse' and c['target'] is None)
  if coarse['key'] in p['reused']:coarsefolder=Path(p['reuse_parent'])/p['reused'][coarse['key']]['job']['key']
  else:coarsefolder=FOLDER/coarse['key']
  coarsehistory=coarsefolder/'history.json';sha(coarsehistory,proof['evidence_sha256'][str(coarsehistory)]);crow=read(coarsehistory)[0]
  assert all(np.array_equal(crow[k],row[k]) for k in ('ids','chemistry','volumes','delta','polarity','centers'))
  k=-volumes[:,None]*delta;s=-np.sqrt(volumes)[:,None]*delta/np.sqrt(volumes)[None,:]
  rowzero=float(np.abs(delta@np.ones(16)).max());masszero=float(np.abs(volumes@delta).max());symerror=float(np.abs(k-k.T).max());lambdas=np.linalg.eigvalsh((s+s.T)/2)
  assert rowzero<=1e-11 and masszero<=1e-12 and symerror<=1e-12 and lambdas[0]>=-1e-11
  np.testing.assert_allclose(lambdas,row['laplacian_lambdas'],rtol=1e-10,atol=1e-11)
  array_bytes=ids.nbytes+chemistry.nbytes+volumes.nbytes+delta.nbytes
  contexts.append(dict(seed=j['seed'],background=j['background'],recipient=j['recipient'],targets=p['selections'][str(j['seed'])]['selected_ids'] if j['recipient'] is None else [j['recipient']],job=j['key'],reused_sham=reused,source=str(source),source_sha256=evidence[str(source)],history=str(historypath),history_sha256=evidence[str(historypath)],result=str(resultpath),result_sha256=evidence[str(resultpath)],parent_endpoint_history=str(parent_history),parent_endpoint_history_sha256=evidence[str(parent_history)],coarse_initial_graph_verified=True,time=510.,elapsed=0.,ids=ids.tolist(),chemistry=chemistry.tolist(),volumes=volumes.tolist(),delta=delta.tolist(),config={key:meta['config'][key] for key in ('signal_beta','signal_da','signal_dh','dt')},precision=row['precision'],array_hashes={key:arraysha(value) for key,value in (('ids',ids),('chemistry',chemistry),('volumes',volumes),('delta',delta))},small_array_bytes=array_bytes,source_compressed_bytes=source.stat().st_size,row_sum_error=rowzero,mass_conservation_error=masszero,conductance_symmetry_error=symerror,min_laplacian_eigenvalue=float(lambdas[0])))
  print('Verified exact frozen context',j['seed'],j['background'],flush=True)
 # Tiny state/operator records replace phase-field copies. 66 trajectories at
 # 401 observations and 32 chemical values require only 6.46 MiB uncompressed.
 frozen_context_bytes=sum(c['small_array_bytes'] for c in contexts);trajectory_bytes=66*401*2*16*8
 prior_review_pins={name:digest for name,digest in evidence.items() if name.startswith(str(reviewroot)) or name==str(compact_path)}
 selected_input_pins={name:digest for name,digest in evidence.items() if name not in p['source_sha256'] and name not in prior_review_pins}
 for name,digest in list(evidence.items()):assert sha(name,digest)==digest
 for name,value in stats.items():assert metadata(name)==value
 report.update(passed=True,protocol_sha256=ph,parent_status=status,parent_summary_accepted=True,compact_assessment=str(compact_path),compact_has_full_assessment_key='full_assessment' in compact,review_archive=str(reviewroot),review_schema='Accepted compact assessment + immutable archive review-provenance, provenance_verification, numerical_verification and response_reconstruction; no full_assessment field.',source_hashes_verified=len(p['source_sha256']),runtime_binary_hashes_verified=2,input_evidence_paths_verified=len(selected_input_pins),scientific_review_pins=prior_review_pins,source_pins=p['source_sha256'],selected_parent_input_pins=selected_input_pins,excluded_reporting_pins=excluded,contexts=contexts,contexts_verified=9,initial_coarse_fine_equivalence_contexts=9,initial_parent_endpoint_equivalence_contexts=9,estimated_small_state_operator_bytes=frozen_context_bytes,estimated_66_raw_chemical_trajectories_bytes=trajectory_bytes,minimum_state_operator_total_mib=frozen_context_bytes/2**20,raw_trajectory_total_mib=trajectory_bytes/2**20,workspace_disk=shutil.disk_usage(ROOT)._asdict(),relocated_outputs_disk=shutil.disk_usage(ROOT/'outputs')._asdict(),suggested_workspace_stage=str(ROOT/'archive/study-experiments/2026-10-08-frozen-redistribution-response'),stage_is_within_workspace=True,stage_requires_external_output_write=False,original_outputs_alias=str(APPROVED),runtime=dict(torch='2.7.1+cu118',torch_cuda='11.8',qualified_gpu=1,qualified_device_name='NVIDIA GeForce GTX 1080 Ti',qualified_gpu_uuid='GPU-2991ab3d-dcbd-41f3-eac1-133be278f56e',free_gpu_mib=11159,compute_jobs_before_check=0,cuda_available=True,float64_16x16_matrix_exp_availability_check_passed=True,scope='Read-only nvidia-smi and one tiny PyTorch GPU operation; no assay launched.'),evidence_sha256=evidence,immutable_after_recheck=True,all_mismatches=[])
except BaseException as exc:
 import traceback;report.update(error=str(exc),traceback=traceback.format_exc(),all_mismatches=[str(exc)]);raise
finally:
 report['finished_utc']=datetime.now(timezone.utc).isoformat();report['elapsed_seconds']=time.monotonic()-started
 REPORT.write_text(json.dumps(report,indent=2)+'\n');print('Read-only frozen input/storage review passed',report['passed'],REPORT,flush=True)
