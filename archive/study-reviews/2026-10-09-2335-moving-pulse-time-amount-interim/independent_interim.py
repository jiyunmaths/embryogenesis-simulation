"""Read-only interim reconstruction; no embryo imports or job mutations."""
from pathlib import Path
from datetime import datetime,timezone,timedelta
from zoneinfo import ZoneInfo
import json,hashlib,re,time
import numpy as np
ROOT=Path('/home/sji/Documents/Programming/Python/embryogenesis-simulation')
R=ROOT/'output/studies/moving-pulse-time-amount-response'
L=ROOT/'archive/study-launches/2026-10-09-moving-pulse-time-amount-response'
WATCHED={}
def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
def read(f,stable=True):
 f=Path(f);raw=f.read_bytes();h=hashlib.sha256(raw).hexdigest()
 if stable:WATCHED[str(f.resolve())]=h
 return json.loads(raw)
def history(job,pilot=True):
 h=read(R/job['key']/('pilot-history.json' if pilot else 'history.json'))
 u=np.asarray([x['chemistry'] for x in h]);times=np.asarray([x['elapsed'] for x in h])
 expected=6. if pilot else 60.
 assert len(h)==round(expected/.15)+1 and np.allclose(times,np.arange(len(h))*.15,atol=1e-10,rtol=0)
 assert all(abs(x['time']-570.-t)<1e-9 for x,t in zip(h,times))
 assert np.all(u>0) and np.isfinite(u).all()
 with np.load(job['source'],allow_pickle=False) as z:
  assert np.array_equal(u[0],z['initial']) and h[0]['ids']==z['ids'].tolist()
 return h,u,times

def main():
 now=datetime.now(timezone.utc);p=read(R/'protocol.json');ph=sha(R/'protocol.json');status=read(R/'status.json',False);launch=read(L/'process-launch.json');q=read(R/'qualification.json');qa=read(L/'qualification-audit.json');pilot=read(R/'pilot-gate.json')
 assert q['passed'] and q['jobs']==48 and qa['passed'] and qa['actual_start_gates']==48 and q['protocol_sha256']==qa['protocol_sha256']==ph
 assert pilot['passed'] and pilot['contexts']==12 and pilot['contrast_pairs']==6 and pilot['protocol_sha256']==ph
 for name in ('source_sha256','input_sha256'):
  for f,h in p[name].items():assert sha(f)==h,f
 for f,h in pilot['evidence_sha256'].items():assert sha(f)==h,f
 rows=[];max_response_error=0.;max_signed_error=0.;max_mechanical={};checks=0
 for c in p['contexts']:
  rep=read(R/'pairs'/c['key']/'response-pilot.json');assert rep['passed'] and rep['protocol_sha256']==ph
  for row in rep['mechanical']:
   for k,v in row['errors'].items():max_mechanical[k]=max(max_mechanical.get(k,0),v)
  paths={}
  for level in ('coarse','fine'):
   for dose in (None,'fractional'):
    job=next(j for j in p['jobs'] if (j['context'],j['level'],j['dose'])==(c['key'],level,dose));h,u,t=history(job);paths[level,dose]=u;checks+=1
  measured=float(abs(np.log(paths['coarse','fractional']/paths['coarse',None])-np.log(paths['fine','fractional']/paths['fine',None])).max()/abs(np.log(.9)))
  assert abs(measured-rep['response'][0]['normalized_response_error'])<1e-14 and measured<=.01;max_response_error=max(max_response_error,measured)
 for pair in p['pairs']:
  report=read(R/'contrasts'/f'history-{pair["seed"]}_cell-{pair["cell"]}_pilot.json');assert report['passed'];waves={};effects={}
  for level in ('coarse','fine'):
   w=[]
   for C in ('S','R'):
    key=pair['contexts'][C];jp=next(j for j in p['jobs'] if (j['context'],j['level'],j['dose'])==(key,level,'fractional'));jc=next(j for j in p['jobs'] if (j['context'],j['level'],j['dose'])==(key,level,None))
    h,u,t=history(jp);_,b,_=history(jc);i=h[0]['ids'].index(pair['cell']);w.append(np.log(u[:,:,i]/b[:,:,i])/abs(np.log(.9)))
   waves[level]=w[1]-w[0];effects[level]=float(np.sqrt(np.trapezoid(np.mean(waves[level]**2,axis=1),t)/(t[-1]-t[0])))
  discrepancy=float(abs(waves['coarse']-waves['fine']).max());row=report['rows'][0]
  assert effects==row['effects'] and discrepancy==row['direct_signed_waveform_error'];assert discrepancy<=.01 and (effects['coarse']>.01)==(effects['fine']>.01)
  max_signed_error=max(max_signed_error,discrepancy);rows.append(dict(history=pair['seed'],cell=pair['cell'],pilot_rms_response_contrast=effects['fine'],above_pilot_screen=effects['fine']>.01,timestep_signed_error=discrepancy))
 completed=[];quality=[];progress=[]
 for j in p['jobs']:
  f=R/j['key']/'result.json';live=read(R/j['key']/'status.json',False)
  elapsed=float(live['elapsed']);progress.append(dict(job=j['key'],level=j['level'],elapsed=elapsed))
  if f.exists():
   rr=read(f);assert rr['quality_pass'] and rr['job']==j and rr['protocol_sha256']==ph
   for k,name in [('history_sha256','history.json'),('checkpoint_sha256','latest_state.npz')]:assert sha(R/j['key']/name)==rr[k]
   h,u,t=history(j,False);completed.append(dict(job=j['key'],level=j['level'],wall_seconds=rr['audit']['wall_seconds'],final_log_activator_sd=float(np.std(np.log(u[-1,0]))),whole_window_min_log_activator_sd=float(np.std(np.log(u[:,0]),axis=1).min())))
   a=rr['audit']
  else:a=live['audit']
  assert a['max_volume_error']<.05 and a['min_radius']>=4 and a['max_clipping']==0 and a['boundary_max']<.01 and a['dilution_error_max']<=2e-14
  quality.append(a)
 # Independently reconstruct each completed full background and each fully
 # complete eight-path recipient comparison. Never read unfinished histories.
 full_pairs=[];full_response_rows=[];recovery=[];full_normalized_max=0.;full_raw_max=0.;full_response_auc_max=0.;full_mechanical_max={}
 def recovered(t,dev,threshold):
  suffix=np.logical_and.accumulate((dev<=threshold)[::-1])[::-1]
  eligible=np.flatnonzero(suffix & (t<=t[-1]-24.))
  return float(t[eligible[0]]) if len(eligible) else None
 for f in sorted((R/'pairs').glob('*/response-full.json')):
  rep=read(f);assert rep['passed'] and rep['protocol_sha256']==ph;key=rep['context'];paths={}
  for level in ('coarse','fine'):
   for dose in (None,'fractional'):
    job=next(j for j in p['jobs'] if (j['context'],j['level'],j['dose'])==(key,level,dose));hh,u,tt=history(job,False);paths[level,dose]=(hh,u,tt)
  a,b,c,d=[paths[level,dose][1] for level,dose in [('coarse','fractional'),('coarse',None),('fine','fractional'),('fine',None)]]
  error=float(abs(np.log(a/b)-np.log(c/d)).max()/abs(np.log(.9)))
  assert abs(error-rep['response'][0]['normalized_response_error'])<1e-14 and error<=.01;full_normalized_max=max(full_normalized_max,error)
  for measured in rep['mechanical']:
   dose=measured['dose'];h0,x,t=paths['coarse',dose];h1,y,_=paths['fine',dose]
   raw=float(abs(np.log(x/y)).max());assert abs(raw-measured['errors']['chemical_log_max'])<1e-14;full_raw_max=max(full_raw_max,raw)
   for k,v in measured['errors'].items():assert v<=p['mechanical_limits'][k];full_mechanical_max[k]=max(full_mechanical_max.get(k,0.),v)
  context=next(c for c in p['contexts'] if c['key']==key);cell=context['recipient'];target_index=paths['fine',None][0][0]['ids'].index(cell);auc=[];recover_levels={}
  for level in ('coarse','fine'):
   hh,u,tt=paths[level,'fractional'];bb,base,_=paths[level,None];v=np.asarray(hh[0]['volumes']);log=np.log(u/base);amp=abs(np.log(.9));network=np.sqrt(np.sum((np.log(u)-np.log(base))**2*v,axis=(-1,-2))/(2*v.sum()));target=np.abs(log[:,:,target_index]).max(1)
   auc.append(float(np.trapezoid(abs(log[:,0,target_index]),tt)/amp));rec=dict(target=recovered(tt,target,.1*amp),network=recovered(tt,network,.1*network[0]));recover_levels[level]=rec
   if level=='fine':recovery.append(dict(history=context['seed'],cell=cell,chemical_origin=context['chemical_origin'],**rec,final_network_log_rms=float(network[-1])))
  aucerr=abs(auc[1]/auc[0]-1);assert abs(aucerr-rep['response'][0]['relative_auc_error'])<1e-14 and aucerr<=.02;full_response_auc_max=max(full_response_auc_max,aucerr)
  for kind in ('target','network'):
   x,y=[recover_levels[l][kind] for l in ('coarse','fine')];e=0. if x is None and y is None else None if x is None or y is None else abs(x-y)
   assert e==rep['response'][0]['recovery_time_errors'][kind+'_recovery_time'] and e is not None and e<=.3+1e-9
  full_pairs.append(dict(context=rep['context'],normalized_response_error=error,relative_auc_error=aucerr,recovery_times=recover_levels,mechanical_errors=[x['errors'] for x in rep['mechanical']]))
 frozen_review_file=ROOT/'archive/study-reviews/2026-10-09-pulse-time-amount-completed/independent_review.json';fr=read(frozen_review_file)
 assert fr['passed'] and fr['accepted_at_tested_scope']
 for pair in p['pairs']:
  relevant=[j for j in p['jobs'] if j['context'] in pair['contexts'].values()]
  if len(relevant)!=8 or not all((R/j['key']/'result.json').exists() for j in relevant):continue
  # Only classify if both contexts' complete timestep reports passed.
  if not all(any(x['context']==key for x in full_pairs) for key in pair['contexts'].values()):continue
  waves={};effects={}
  for level in ('coarse','fine'):
   responses=[]
   for C in ('S','R'):
    jp=next(j for j in relevant if (j['context'],j['level'],j['dose'])==(pair['contexts'][C],level,'fractional'));jc=next(j for j in relevant if (j['context'],j['level'],j['dose'])==(pair['contexts'][C],level,None))
    hh,u,tt=history(jp,False);_,b,_=history(jc,False);i=hh[0]['ids'].index(pair['cell']);responses.append(np.log(u[:,:,i]/b[:,:,i])/abs(np.log(.9)))
   waves[level]=responses[1]-responses[0];effects[level]=float(np.sqrt(np.trapezoid(np.mean(waves[level]**2,axis=1),tt)/(tt[-1]-tt[0])))
  raw=float(abs(waves['coarse']-waves['fine']).max());err=abs(effects['fine']-effects['coarse']);same=(effects['fine']>.01)==(effects['coarse']>.01)
  assert same and raw<=.01 and err<=.01
  previous=next(x for x in fr['scientific_comparisons'] if (x['seed'],x['cell'],x['matching_reference'],x['budget_reference'],x['physics'])==(pair['seed'],pair['cell'],'S','S','S'))
  frozen=previous['full_effects']['fine'];full_response_rows.append(dict(history=pair['seed'],cell=pair['cell'],moving_effects=effects,above_screen=effects['fine']>.01,frozen_effect=frozen,frozen_above_screen=previous['above_screen'],classification_agrees_with_frozen=(effects['fine']>.01)==previous['above_screen'],descriptive_relative_rms_change=effects['fine']/frozen-1,signed_timestep_error=raw,rms_timestep_error=err,numerical_gates_pass=True))
 log=(L/'run.log').read_text();timings={}
 for line in log.splitlines():
  match=re.fullmatch(r'(seed-\S+) ([0-9.]+) ([0-9.]+) ([0-9.]+)',line)
  if match:
   key,e,h,w=match.groups();timings.setdefault(key,{})[float(h)]=dict(elapsed=float(e),wall_seconds=float(w))
 calibration={}
 for level in ('coarse','fine'):
  ratios=[]
  for j in p['jobs']:
   if j['level']!=level or j['key'] not in timings:continue
   z=timings[j['key']];base=z.get(6.);full=z.get(60.)
   if base and full and full['elapsed']>12.:
    pilot_rate=base['wall_seconds']/6.;recent_rate=(full['wall_seconds']-base['wall_seconds'])/(full['elapsed']-6.)
    ratios.append(recent_rate/pilot_rate)
  calibration[level]=float(np.median(ratios)) if ratios else 1.
 seconds_remaining=0.;estimates=[]
 for j in p['jobs']:
  if (R/j['key']/'result.json').exists():continue
  live=read(R/j['key']/'status.json',False);elapsed=float(live['elapsed']);base=timings[j['key']][6.];rate=base['wall_seconds']/6.*calibration[j['level']]
  remain=max(0.,60.-elapsed)*rate;seconds_remaining+=remain;estimates.append(dict(job=j['key'],elapsed=elapsed,seconds_per_model_unit=rate,remaining_seconds=remain))
 # Initialization, checkpoint compression, paired assessments/array review add
 # overhead not fully represented by recorded per-step wall timings.
 overhead_seconds=5.*len(estimates)+300.;central=seconds_remaining+overhead_seconds;eta=now+timedelta(seconds=central)
 quality_extrema=dict(maximum_volume_error=max(x['max_volume_error'] for x in quality),minimum_radius_grid_intervals=min(x['min_radius'] for x in quality),maximum_boundary_occupancy=max(x['boundary_max'] for x in quality),maximum_dilution_amount_error=max(x['dilution_error_max'] for x in quality),maximum_clipping=max(x['max_clipping'] for x in quality))
 for f,h in WATCHED.items():assert sha(f)==h,f
 report=dict(observed_utc=now.isoformat(),protocol_sha256=ph,status_snapshot=status,independent_histories=3,new_histories=0,qualification_passed=48,qualification_seconds=q['seconds'],independent_prefix_audit_passed=True,source_pins_verified=len(p['source_sha256']),input_pins_verified=len(p['input_sha256']),pilot_paths_verified=checks,pilot_context_pairs_passed=12,pilot_signed_contrast_pairs_passed=6,pilot_response_rows=rows,pilot_max_normalized_response_error=max_response_error,pilot_max_signed_contrast_error=max_signed_error,pilot_mechanical_maxima=max_mechanical,completed_paths=len(completed),total_paths=48,completed_path_details=completed,completed_full_context_pairs=full_pairs,full_response_rows=full_response_rows,completed_fine_pulse_recovery=recovery,full_max_normalized_response_error=full_normalized_max,full_max_raw_chemical_error=full_raw_max,full_max_relative_auc_error=full_response_auc_max,full_mechanical_maxima=full_mechanical_max,quality_extrema=quality_extrema,simulated_work_fraction=sum(x['elapsed'] for x in progress)/(48*60),trajectory_progress=progress,timing_calibration=calibration,remaining_integrator_seconds=seconds_remaining,additional_overhead_allowance_seconds=overhead_seconds,estimated_remaining_hours=central/3600,estimated_completion_utc=eta.isoformat(),estimated_completion_local=eta.astimezone(ZoneInfo('America/New_York')).isoformat(),eta_assumptions='Single existing GPU worker, current speed retained, pilot-rate estimates calibrated against measured full paths; all jobs retain 6-unit pilot progress. Explicit overhead allowance is heuristic, not a confidence interval.',elapsed_wall_hours=(now-datetime.fromisoformat(launch['utc'])).total_seconds()/3600,source_script_sha256=sha(__file__),evidence_sha256=WATCHED,full_scientific_conclusion_pending=True,scope='Independent literal interim response reconstruction from all stable pilot chemistry histories; completed path checks and live-status quality snapshots. No full-horizon S/R classification until both complete controls and pulses at both timesteps exist; no new histories or identity acceptance.')
 out=ROOT/'archive/study-reviews/2026-10-09-2335-moving-pulse-time-amount-interim';out.mkdir(exist_ok=True);(out/'interim.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');(out/'independent_interim.py').write_text(Path(__file__).read_text())
 print(json.dumps({k:report[k] for k in ['observed_utc','status_snapshot','qualification_passed','pilot_paths_verified','pilot_response_rows','pilot_max_normalized_response_error','pilot_max_signed_contrast_error','completed_paths','full_response_rows','completed_fine_pulse_recovery','full_max_normalized_response_error','full_max_raw_chemical_error','full_max_relative_auc_error','quality_extrema','simulated_work_fraction','estimated_remaining_hours','estimated_completion_local','elapsed_wall_hours']},indent=2))
main()
