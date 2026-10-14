"""Read-only NumPy interim assessment; no scientific solver imports or stepping."""
from datetime import datetime, timezone, timedelta
from zoneinfo import ZoneInfo
from pathlib import Path
import hashlib,json,time,shutil
import numpy as np

ROOT=Path.cwd(); STUDY=ROOT/'output/studies/fresh-zygote-maintenance'
EXPECTED='78f7a98bdb2a3d0ecdaa825f94898d9111b664ba197654014d2ecc02ed656bfd'
inputs={}
def read(path):
 path=Path(path)
 for i in range(20):
  raw=path.read_bytes()
  try:result=json.loads(raw)
  except json.JSONDecodeError:
   if i==19:raise
   time.sleep(.025);continue
  inputs[str(path.absolute())]=hashlib.sha256(raw).hexdigest()
  return result
p=read(STUDY/'protocol.json'); assert inputs[str((STUDY/'protocol.json').absolute())]==EXPECTED
status=read(STUDY/'status.json'); launch=read(ROOT/'archive/study-launches/2026-10-13-fresh-zygote-maintenance/process-launch.json')
now=datetime.now(timezone.utc); local=now.astimezone(ZoneInfo('America/New_York'))
completed=[]; active=[];pairs=[];prefixes=[];max_diagnostic_errors={}
def measure(rows,cfg,complete):
 times=np.array([x['elapsed'] for x in rows]);u=np.array([x['chemistry'] for x in rows]);v=np.array([x['volumes'] for x in rows]);delta=np.array([x['delta'] for x in rows]);a=np.array([x['axis_ratio'] for x in rows])
 assert u.shape==(len(rows),2,16) and np.isfinite(u).all() and u.min()>0 and (v>0).all()
 sd=np.std(np.log(u[:,0,:]),axis=1); assert np.max(np.abs(sd-np.array([x['log_activator_sd'] for x in rows])))<1e-12
 assert np.allclose(times,np.arange(len(rows))*.15,atol=1e-9,rtol=0)
 assert all(abs(x['time']-(78+x['elapsed']))<1e-9 and abs(x['total_elapsed']-(60+x['elapsed']))<1e-9 for x in rows)
 assert all(x['ids']==rows[0]['ids'] and x['parents']==rows[0]['parents'] and x['nondividing_components']==[1]*16 and not x['dividing_ids'] for x in rows)
 assert max(x['boundary_occupancy'] for x in rows)<p['physical_limits']['boundary_max']
 G=v[:,:,None]*delta;G[:,np.arange(16),np.arange(16)]=0
 symmetry=float(abs(G-G.transpose(0,2,1)).max());constant=float(abs(delta.sum(axis=2)).max());conservation=float(abs(np.einsum('ti,tij->tj',v,delta)).max())
 assert symmetry<1e-12 and constant<1e-10 and conservation<1e-12 and G.min()>-1e-12
 spectrum=np.linalg.eigvalsh(-np.sqrt(v)[:,:,None]*delta/np.sqrt(v)[:,None,:]);assert spectrum.min()>-1e-10
 recorded=np.array([x['spectrum'] for x in rows]);assert np.max(abs(spectrum-recorded))<1e-9
 jac=np.array([[1.,-1.],[2*cfg['signal_beta'],-cfg['signal_beta']]])
 blocks=np.broadcast_to(jac,(len(rows),16,2,2)).copy();blocks[:,:,0,0]-=cfg['signal_da']*spectrum;blocks[:,:,1,1]-=cfg['signal_dh']*spectrum
 rates=np.linalg.eigvals(blocks).real.max(axis=2);support=spectrum>1e-10;counts=((rates>1e-10)&support).sum(axis=1);growth=np.where(support,rates,-np.inf).max(axis=1)
 assert np.array_equal(counts,[x['frozen_uniform_unstable_modes'] for x in rows])
 assert np.max(abs(growth-np.array([x['frozen_uniform_max_growth'] for x in rows])))<1e-9
 conductance=G.sum(axis=(1,2))/2
 outcomes=[]
 for theta in p['thresholds']:
  above=sd>theta;late=times>=156-1e-9;transitions=np.flatnonzero(above[1:]!=above[:-1])+1
  outcomes.append(dict(threshold=theta,continuous_retention_so_far=bool(above.all()),late_persistent=bool(above[late].all()) if complete else None,sampled_transitions=[dict(elapsed=float(times[i]),direction='up' if above[i] else 'down') for i in transitions]))
 nonpositive=np.flatnonzero(counts==0)
 startzero=float(times[nonpositive[0]]) if len(nonpositive) else None
 zeroepisodes=[]
 for i in np.flatnonzero(np.r_[True,counts[1:]!=counts[:-1]]):
  if counts[i]==0:zeroepisodes.append(float(times[i]))
 return dict(rows=len(rows),physical_end=float(rows[-1]['time']),maintenance_elapsed=float(times[-1]),total_elapsed_since_t18=float(rows[-1]['total_elapsed']),initial_log_sd=float(sd[0]),minimum_log_sd=float(sd.min()),final_log_sd=float(sd[-1]),late_min_log_sd=float(sd[late].min()) if complete else None,outcomes=outcomes,initial_undirected_conductance=float(conductance[0]),final_undirected_conductance=float(conductance[-1]),conductance_fractional_change=float(conductance[-1]/conductance[0]-1),initial_growing_modes=int(counts[0]),final_growing_modes=int(counts[-1]),minimum_growing_modes=int(counts.min()),first_sampled_zero_growing_modes_elapsed=startzero,zero_mode_episode_starts_elapsed=zeroepisodes,final_frozen_uniform_max_growth=float(growth[-1]),initial_axis_ratio=float(a[0]),final_axis_ratio=float(a[-1]),relative_axis_change=float(a[-1]/a[0]-1),max_conductance_asymmetry=symmetry,max_constant_error=constant,max_amount_conservation_error=conservation)
for job in p['jobs']:
 folder=STUDY/job['key'];rp=folder/'gpu-result.json'; pp=folder/'actual-start-prefix.json'
 if pp.exists():
  q=read(pp);assert q['passed'] and q['exact_gpu_restart'] and q['protocol_sha256']==EXPECTED;assert all(q['errors'][k]<=p['prefix_limits'][k] for k in q['errors']);prefixes.append(dict(case=job['key'],seconds=q['seconds'],errors=q['errors']))
 if not (folder/'gpu-history.json').exists():continue
 rows=read(folder/'gpu-history.json')
 if rp.exists():
  report=read(rp);assert report['passed'] and report['protocol_sha256']==EXPECTED and report['final_time']==258
  assert report['exact_restart_steps']=={'36':8,'108':8} and report['lineage_cycle_rng_unchanged']
  assert len(rows)==1201
  m=measure(rows,job['config'],True)
  for x,y in zip(m['outcomes'],report['maintenance']):
   assert x['continuous_retention_so_far']==y['continuous_retention'] and x['late_persistent']==y['late_persistent'] and x['sampled_transitions']==[dict(elapsed=z['elapsed'],direction=z['direction']) for z in y['sampled_transitions']]
  completed.append(dict(case=job['key'],history=job['seed'],branch=job['branch'],level=job['level'],seconds=report['seconds'],audit=report['audit'],measurement=m))
 else:
  progress=read(folder/'gpu-progress.json');assert progress['state']=='running'
  active.append(dict(case=job['key'],history=job['seed'],branch=job['branch'],level=job['level'],progress=progress,measurement=measure(rows,job['config'],False)))
for pair in p['time_pairs']:
 path=STUDY/f"time-h{pair['seed']}-{pair['branch']}.json"
 if path.exists():
  report=read(path);assert report['passed'] and not report['failed_gates'] and report['protocol_sha256']==EXPECTED
  for name,val in report['errors'].items():assert val<=p['time_limits'][name]
  for name,val in report['homologous_errors'].items():assert val<=p['homologous_limits'][name]
  assert report['field_max']<=p['time_limits']['field_max'] and report['dense_peak']['relative_spectrum']<=p['dense_time_spectrum_limit'] and report['growth_abs_max']<=p['growth_abs_limit']
  assert all(x['passed'] for x in report['maintenance_checks'])
  pairs.append(dict(history=pair['seed'],branch=pair['branch'],field_max=report['field_max'],dense_peak=report['dense_peak'],homologous_errors=report['homologous_errors'],growth_abs_max=report['growth_abs_max'],topology_mismatch_times=report['topology_mismatch_times']))
assert len(completed)==status['completed'] and len(pairs)==status['pairs_completed'] and len(active)==1
assert status['state']=='running_maintenance' and active[0]['case']==status['active']
means={level:float(np.mean([x['seconds'] for x in completed if x['level']==level])) for level in ('coarse','fine')}
prefixmean=float(np.mean([x['seconds'] for x in prefixes])); completedseconds=sum(x['seconds'] for x in completed)
# Status is written after the active prefix. Residual overhead includes restore, cleanup, comparison and bookkeeping; it is not isolated comparison runtime.
comparison_overhead=status['seconds']-completedseconds-sum(x['seconds'] for x in prefixes)
assert comparison_overhead>=0
perpair=comparison_overhead/len(pairs)
bykey={x['case']:x for x in completed}; activekey=active[0]['case'];ap=active[0]['progress']
remaining=0
for job in p['jobs']:
 if job['key'] in bykey:continue
 if job['key']==activekey:
  remaining+=max(0.,means[job['level']]-ap['seconds'])
 else:remaining+=means[job['level']]+prefixmean
remaining+=(12-len(pairs))*perpair
eta=now+timedelta(seconds=remaining)
launchat=datetime.fromisoformat(launch['utc']);elapsed=(now-launchat).total_seconds()
result=dict(observed_utc=now.isoformat(),observed_local=local.isoformat(),protocol_sha256=EXPECTED,state=status['state'],completed_paths=len(completed),total_paths=24,completed_pairs=len(pairs),total_pairs=12,completed_history_contexts=sorted(set((x['history'],x['branch']) for x in completed)),independent_histories_planned=3,independent_histories_with_completed_pairs=sorted(set(x['history'] for x in pairs)),new_scientific_histories=0,prefixes_checked=len(prefixes),completed=completed,active=active,paired_reports=pairs,elapsed_wall_hours=elapsed/3600,mean_path_seconds=means,mean_prefix_seconds=prefixmean,completed_unattributed_coordinator_overhead_seconds=comparison_overhead,mean_overhead_seconds_per_completed_pair=perpair,estimated_remaining_hours=remaining/3600,estimated_completion_local=eta.astimezone(ZoneInfo('America/New_York')).isoformat(),estimated_total_hours=(elapsed+remaining)/3600,conditional_remaining_range_hours=[20,23],external_free_gib=shutil.disk_usage(STUDY).free/2**30,scientific_accepted=False,independent_final_review_pending=True,scope='Read-only interim reconstruction of saved regular chemistry, spectra, geometry and maintenance outcomes; completed production comparison reports/prefix/restart records checked without rerunning them. Dense-spectrum/native-field gates are committed production reports, not independently reconstructed here. Active data are a timestamped prefix of the trajectory. No final qualification, full source/input rehash, checkpoint-native-field review, biological identity or whole-system stability claim.',consumed_file_snapshots_sha256=inputs)
Path('/tmp/fresh_maintenance_interim.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['observed_local','completed_paths','completed_pairs','prefixes_checked','estimated_remaining_hours','estimated_completion_local','estimated_total_hours','external_free_gib']},indent=2))
for x in completed:
 if x['level']=='fine':print('FINE',x['branch'],json.dumps(x['measurement']))
print('ACTIVE',json.dumps(active))
