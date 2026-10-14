from pathlib import Path
import json, math, numpy as np
from scipy.integrate import trapezoid
ROOT=Path('/home/sji/Documents/Programming/Python/embryogenesis-simulation/outputs/phase-carry-redistribution-response')
P=json.loads((ROOT/'protocol.json').read_text()); JOBS=P['jobs']; TIME=np.arange(401)*.15
OUT=Path('/tmp/redistribution_numerical_audit.json')
failures=[]; disagreements=[]; observations={}; path_reports=[]; mechanical=[]; pulse_reports=[]; effects=[]

def require(flag, label):
 if not bool(flag): failures.append(label)

def numeric_match(actual, stored, label, atol=1e-12, rtol=1e-11):
 if actual is None or stored is None:
  ok=actual is None and stored is None
 elif isinstance(actual,(bool,np.bool_)) or isinstance(stored,(bool,np.bool_)):
  ok=bool(actual)==bool(stored)
 else:
  ok=np.isclose(actual,stored,atol=atol,rtol=rtol)
 if not ok: disagreements.append(dict(label=label,reconstructed=actual,recorded=stored))

def folder(j):
 if j['key'] in P['reused']:
  return Path(P['reuse_parent'])/P['reused'][j['key']]['job']['key']
 return ROOT/j['key']

def phi_file(j,label):
 return Path(P['reused'][j['key']][label+'_phi']) if j['key'] in P['reused'] else ROOT/'pairs'/f'seed-{j["seed"]}'/f'{label}_{j["key"]}.npz'

def key(j): return (j['seed'],j['level'],j['background'],j['target'],j['dose'])

def recovery(t,d,threshold):
 valid=np.logical_and.accumulate((np.asarray(d)<=threshold)[::-1])[::-1]
 index=np.flatnonzero(valid & (t<=t[-1]-24.))
 return float(t[index[0]]) if len(index) else None

def contrast(o,end):
 count=round(end/.15)+1; times=o['elapsed'][:count]; sd=np.std(np.log(o['chemistry'][:count,0,:]),axis=1)
 idx=np.flatnonzero((sd[:-1]>.1)&(sd[1:]<=.1)); loss=None
 if len(idx):
  i=idx[0]; loss=float(times[i]+(.1-sd[i])*(times[i+1]-times[i])/(sd[i+1]-sd[i]))
 late=times>=end-min(24.,end)-1e-9
 return dict(final_log_sd=float(sd[-1]),late_min_log_sd=float(sd[late].min()),persistent_contrast=bool(sd[late].min()>.1),whole_window_min_log_sd=float(sd.min()),whole_window_retained=bool(sd.min()>.1),first_loss=loss)

def waveform_distance(a,b,t): return float(np.sqrt(trapezoid(np.mean((a-b)**2,axis=1),t)/(t[-1]-t[0])))

def response(o,control,j,end):
 count=round(end/.15)+1; t=o['elapsed'][:count]; x=o['chemistry'][:count]; c=control['chemistry'][:count]
 i=o['ids'][0].tolist().index(j['target']); f=float(x[0,0,i]/c[0,0,i]); amp=abs(np.log(j['factor']))
 require(np.array_equal(x[0,0,i],c[0,0,i]*j['factor']),j['key']+' exact prescribed initial target multiplication')
 numeric_match(f,j['factor'],j['key']+' measured factor',atol=1e-16,rtol=1e-15)
 raw=np.log(x/c); volume=control['volumes'][0]; target=np.max(np.abs(raw[:,:,i]),axis=1)
 # Match the inherited network metric's separate-log arithmetic explicitly.
 network=np.sqrt(np.sum((np.log(x)-np.log(c))**2*volume,axis=(1,2))/(2*volume.sum()))
 others=np.arange(len(volume))!=i
 neighbor=np.sqrt(np.sum(raw[:,:,others]**2*volume[others],axis=(1,2))/(2*volume[others].sum()))
 metrics=dict(target_activator_peak_gain=float(np.max(np.abs(raw[:,0,i]))/amp),target_inhibitor_peak_gain=float(np.max(np.abs(raw[:,1,i]))/amp),target_activator_log_auc_per_log_pulse=float(trapezoid(np.abs(raw[:,0,i]),t)/amp),other_cells_peak_log_rms_per_log_pulse=float(neighbor.max()/amp),target_recovery_time=recovery(t,target,.1*amp),network_recovery_time=recovery(t,network,.1*network[0]),final_log_rms_from_control=float(network[-1]),late_max_log_rms_from_control=float(network[t>=t[-1]-24].max()))
 return dict(metrics=metrics,raw=raw,local=raw[:,:,i]/amp,amp=amp,factor=f,times=t)

require(len(JOBS)==66 and len({j['key'] for j in JOBS})==66,'66 unique logical jobs')
require(len(P['reused'])==18,'18 reused paths')
for number,j in enumerate(JOBS):
 k=j['key']; f=folder(j); hist=json.loads((f/'history.json').read_text()); result=json.loads((f/'result.json').read_text()); status=json.loads((f/'status.json').read_text())
 require(len(hist)==401,k+' 401 samples'); require(status['state']=='completed' and status['elapsed']==60,k+' completion')
 o={name:np.array([row[name] for row in hist]) for name in ('elapsed','time','ids','chemistry','volumes','centers','polarity','delta','axis_ratio','boundary_occupancy','log_activator_sd')}
 observations[key(j)]=o
 require(np.max(np.abs(o['elapsed']-TIME))<1e-9,k+' elapsed times')
 require(np.max(np.abs(o['time']-510-TIME))<1e-9,k+' physical times')
 for name,a in o.items():require(np.isfinite(a).all(),k+' finite '+name)
 require(np.all(o['chemistry']>0) and np.all(o['volumes']>0) and np.all(o['axis_ratio']>0),k+' positive state')
 require(o['chemistry'].shape==(401,2,16) and o['polarity'].shape==(401,16,3) and o['delta'].shape==(401,16,16),k+' state shapes')
 require(np.all(o['ids']==o['ids'][0]),k+' fixed identities')
 require(np.max(np.abs(np.std(np.log(o['chemistry'][:,0]),axis=1)-o['log_activator_sd']))<1e-12,k+' raw contrast')
 with np.load(j['source'],allow_pickle=False) as z:
  sourcechem=np.array([z['activator'],z['inhibitor']]); target=z['target'].copy(); srcmeta=json.loads(str(z['metadata'])); srcids=z['ids'].copy()
  require(np.array_equal(sourcechem,o['chemistry'][0]),k+' checkpoint initial chemistry')
  require(np.array_equal(srcids,o['ids'][0]),k+' checkpoint identities')
  require(srcmeta['time']==510 and srcmeta['config']['dt']==j['dt'],k+' source clock')
 with np.load(f/'latest_state.npz',allow_pickle=False) as z:
  cpmeta=json.loads(str(z['metadata'])); saved=json.loads(str(z['long_experiment'])); audit=saved['audit']
  require(cpmeta['time']==570 and cpmeta['step_number']==round(570/j['dt']),k+' endpoint clock')
  require(cpmeta['config']==P['accepted_configs'][f'{j["seed"]}_{j["level"]}'],k+' accepted endpoint configuration')
  require(not cpmeta['divisions'] and cpmeta['next_id']==srcmeta['next_id'] and cpmeta['lineage']==srcmeta['lineage'],k+' no division')
  require(saved['history']==hist and audit==result['audit'],k+' embedded history/audit identity')
  require(np.array_equal(np.array([z['activator'],z['inhibitor']]),o['chemistry'][-1]) and np.array_equal(z['polarity'],o['polarity'][-1]),k+' endpoint chemical/polarity identity')
  require(np.array_equal(target,z['target']) and np.array_equal(srcids,z['ids']),k+' endpoint targets/ids')
  phi=z['phi']; require(phi.dtype==np.float32 and phi.shape==(16,72,72,72) and np.isfinite(phi).all() and phi.min()>=0 and phi.max()<=1,k+' phase-field validity')
  carry=z['phase_carry']; require(carry.dtype==np.float64 and carry.shape==phi.shape and np.isfinite(carry).all(),k+' carry validity'); del carry
  require(z['rounding'].dtype==np.int32 and np.all(z['rounding']>=0),k+' rounding counter validity')
  with np.load(phi_file(j,'full'),allow_pickle=False) as q: require(np.array_equal(phi,q['phi']),k+' full phi snapshot identity')
  dx=2*cpmeta['config']['extent']/cpmeta['config']['grid']
  # Independent endpoint volume from stored float32 fields, allowing CPU/GPU rounding in this diagnostic only.
  measured=np.sum(phi*phi*(3-2*phi),axis=(1,2,3),dtype=np.float64)*dx**3
  endpoint_volume_recompute=float(np.max(np.abs(measured/o['volumes'][-1]-1)))
  del phi
 with np.load(phi_file(j,'pilot'),allow_pickle=False) as q:require(q['phi'].shape==(16,72,72,72) and np.isfinite(q['phi']).all(),k+' pilot phi validity')
 for label,count in [('full',401),('pilot',41)]:
  paired=json.loads((ROOT/'pairs'/f'seed-{j["seed"]}'/f'{label}_{k}.json').read_text());require(paired==hist[:count],k+' '+label+' immutable paired history')
 del paired,saved
 V=o['volumes']; D=o['delta']; off=D.copy(); ii=np.arange(16); off[:,ii,ii]=0; scale=np.maximum(np.max(np.abs(D),axis=(1,2)),1.)
 rowres=np.max(np.abs(D.sum(axis=2)),axis=1)/scale
 massres=np.max(np.abs(np.einsum('ti,tij->tj',V,D)),axis=1)/(scale*V.sum(axis=1))
 reciprocal=V[:,:,None]*D; recip=np.max(np.abs(reciprocal-np.swapaxes(reciprocal,1,2)),axis=(1,2))
 require(np.all(off>=0) and np.all(D[:,ii,ii]<=0) and rowres.max()<=1e-12 and massres.max()<=1e-12,k+' conservative positive transport')
 require(np.all(np.isclose(reciprocal,np.swapaxes(reciprocal,1,2),rtol=1e-12,atol=1e-14*scale[:,None,None])),k+' reciprocal flux')
 sampled_volume=float(np.max(np.abs(V/target-1))); sampled_radius=float(np.min((3*V/(4*np.pi))**(1/3)/dx)); sampled_boundary=float(o['boundary_occupancy'].max())
 require(audit['max_volume_error']<.05 and audit['min_radius']>=4 and audit['max_clipping']==0 and audit['boundary_max']<.01 and audit['dilution_error_max']<=2e-14,k+' full every-step physical audit')
 require(sampled_volume<=audit['max_volume_error']+1e-14 and sampled_radius>=audit['min_radius']-1e-13,k+' cumulative audit bounds sampled values')
 numeric_match(sampled_boundary,audit['boundary_max'],k+' independently sampled boundary',atol=1e-18)
 contrasts=contrast(o,60.)
 for name,value in contrasts.items():numeric_match(value,result[name],k+' result '+name)
 path_reports.append(dict(key=k,reused=k in P['reused'],audit=audit,sampled_max_volume_error=sampled_volume,sampled_min_radius=sampled_radius,sampled_boundary=sampled_boundary,endpoint_CPU_volume_relative_discrepancy=endpoint_volume_recompute,transport_row_residual=float(rowres.max()),transport_mass_residual=float(massres.max()),transport_reciprocal_absolute_residual=float(recip.max()),contrast=contrasts))
 if (number+1)%11==0:print('audited paths',number+1,flush=True)

JOB_BY={key(j):j for j in JOBS}
for seed in (7,8,9):
 for horizon,label in [(6.,'pilot'),(60.,'full')]:
  recorded=json.loads((ROOT/'pairs'/f'seed-{seed}'/f'{label}-refinement.json').read_text()); count=round(horizon/.15)+1
  for j in [q for q in JOBS if q['seed']==seed and q['level']=='coarse']:
   fine=JOB_BY[(seed,'fine',j['background'],j['target'],j['dose'])]; co=observations[key(j)]; fi=observations[key(fine)]
   require(np.array_equal(co['elapsed'],fi['elapsed']) and np.array_equal(co['time'],fi['time']),j['key']+' exact coarse/fine time alignment')
   for name in ('chemistry','volumes','centers','polarity','delta','axis_ratio'):
    require(np.array_equal(co[name][0],fi[name][0]),j['key']+' same '+name+' start')
   def a(name):return co[name][:count]
   def b(name):return fi[name][:count]
   errs=dict(chemical_log_max=float(np.max(np.abs(np.log(a('chemistry')/b('chemistry'))))),polarity_abs_max=float(np.max(np.abs(a('polarity')-b('polarity')))),relative_axis_max=float(np.max(np.abs(a('axis_ratio')/b('axis_ratio')-1))),relative_volume_max=float(np.max(np.abs(a('volumes')/b('volumes')-1))),relative_transport_max=float(np.max(np.linalg.norm(a('delta')-b('delta'),axis=(1,2))/np.maximum(np.linalg.norm(b('delta'),axis=(1,2)),1e-300))))
   with np.load(phi_file(j,label),allow_pickle=False) as c,np.load(phi_file(fine,label),allow_pickle=False) as f:errs['final_phi_abs_max']=float(np.max(np.abs(c['phi']-f['phi'])))
   con=[contrast(o,horizon) for o in (co,fi)];loss=[c['first_loss'] for c in con];loss_error=0. if loss==[None,None] else None if None in loss else abs(loss[0]-loss[1]);same=all(con[0][x]==con[1][x] for x in ('persistent_contrast','whole_window_retained'))
   passed=all(v<=P['mechanical_limits'][k] for k,v in errs.items()) and same and loss_error is not None and loss_error<=.3
   require(passed,j['key']+' '+label+' mechanical gates')
   old=next(q for q in recorded['mechanical'] if (q['context'],q['target'],q['dose'])==(j['background'],j['target'],j['dose']))
   for n,v in errs.items():numeric_match(v,old['errors'][n],j['key']+' '+label+' '+n)
   for oldc,newc in zip(old['contrast_summaries'],con):
    for n,v in newc.items():numeric_match(v,oldc[n],j['key']+' '+label+' '+n)
   mechanical.append(dict(key=j['key'],horizon=horizon,passed=bool(passed),errors=errs,contrast=con,loss_error=loss_error))
  caches={}
  for target in P['selections'][str(seed)]['selected_ids']:
   for bg in ('sham',f'redistributed-cell-{target}'):
    for dose in ('fractional','amount'):
     pair=[]
     for lev in ('coarse','fine'):
      j=JOB_BY[(seed,lev,bg,target,dose)];o=observations[key(j)];ctrl=observations[(seed,lev,bg,None,None)]; r=response(o,ctrl,j,horizon); caches[(lev,bg,target,dose)]=r; pair.append(r)
      oldrow=next(x for x in recorded['responses'][lev] if x['cell']==target)['contexts'][bg]['doses'][dose]['metrics']
      for name,value in r['metrics'].items():numeric_match(value,oldrow[name],j['key']+' '+label+' metric '+name)
      if horizon==60:
       i=o['ids'][0].tolist().index(target);v=ctrl['volumes'][0,i];before=ctrl['chemistry'][0,0,i]*v;actual=(ctrl['chemistry'][0,0,i]-o['chemistry'][0,0,i])*v
       changed=np.argwhere(o['chemistry'][0]!=ctrl['chemistry'][0]);require(changed.tolist()==[[0,i]],j['key']+' target-only pulse')
       req=.1*before if dose=='fractional' else .1*min(observations[(seed,lev,'sham',None,None)]['chemistry'][0,0,i]*observations[(seed,lev,'sham',None,None)]['volumes'][0,i],observations[(seed,lev,f'redistributed-cell-{target}',None,None)]['chemistry'][0,0,i]*observations[(seed,lev,f'redistributed-cell-{target}',None,None)]['volumes'][0,i])
       rel=abs(actual/req-1);require(rel<=2e-13,j['key']+' actual dose')
       for name,value in [('initial_target_activator_amount',before),('removed_activator_amount',actual),('requested_removed_amount',req),('initial_log_amplitude',r['amp'])]:numeric_match(value,j[name],j['key']+' '+name,atol=1e-18,rtol=2e-13)
       r['dose_report']=dict(initial_target_amount=float(before),requested=float(req),removed=float(actual),relative_amount_error=float(rel),factor=j['factor'],normalization=r['amp'])
     err=float(np.max(np.abs(pair[0]['raw']-pair[1]['raw']))/pair[0]['amp']); auc=abs(pair[1]['metrics']['target_activator_log_auc_per_log_pulse']/pair[0]['metrics']['target_activator_log_auc_per_log_pulse']-1)
     rec={}
     for name in ('target_recovery_time','network_recovery_time'):
      x,y=[r['metrics'][name] for r in pair]; rec[name]=0. if x is None and y is None else None if None in (x,y) else abs(x-y)
     passed=err<=.01 and auc<=.02 and all(x is not None and x<=.3+1e-9 for x in rec.values());require(passed,f'{seed} {bg} {target} {dose} {label} response gates')
     old=next(q for q in recorded['responses']['rows'] if (q['context'],q['cell'],q['dose'])==(bg,target,dose));numeric_match(err,old['normalized_response_error'],f'{seed} {bg} {target} {dose} {label} response error');numeric_match(auc,old['relative_auc_error'],f'{seed} {bg} {target} {dose} {label} areaerror')
     for name,value in rec.items():numeric_match(value,old['recovery_time_errors'][name],f'{seed} {bg} {target} {dose} {label} '+name)
     pulse_reports.append(dict(seed=seed,background=bg,target=target,dose=dose,horizon=horizon,passed=bool(passed),normalized_response_error=err,relative_auc_error=auc,recovery_time_errors=rec,coarse_metrics=pair[0]['metrics'],fine_metrics=pair[1]['metrics'],coarse_dose=pair[0].get('dose_report'),fine_dose=pair[1].get('dose_report')))
   for dose in ('fractional','amount'):
    shift=[]
    for lev in ('coarse','fine'):
     first=caches[(lev,'sham',target,dose)];second=caches[(lev,f'redistributed-cell-{target}',target,dose)];t=first['times'];e=waveform_distance(first['local'],second['local'],t);raw=waveform_distance(first['local']*first['amp'],second['local']*second['amp'],t);shift.append(e)
     old=next(x for x in recorded['responses'][lev] if x['cell']==target)['differences'][dose]
     numeric_match(e,old['normalized_response_shift_rms'],f'{seed} {target} {dose} {label} {lev} effect');numeric_match(raw,old['raw_log_response_shift_rms'],f'{seed} {target} {dose} {label} {lev} raw effect');numeric_match(e>.01,old['response_shift_detected'],f'{seed} {target} {dose} {label} {lev} decision')
    passed=(shift[0]>.01)==(shift[1]>.01) and abs(shift[0]-shift[1])<=.01;require(passed,f'{seed} {target} {dose} {label} effect gate')
    effects.append(dict(seed=seed,target=target,dose=dose,horizon=horizon,coarse=shift[0],fine=shift[1],error=abs(shift[0]-shift[1]),above_screen=bool(shift[1]>.01),passed=bool(passed)))
   for lev in ('coarse','fine'):
    for bg in ('sham',f'redistributed-cell-{target}'):
     f,a=[caches[(lev,bg,target,d)] for d in ('fractional','amount')]; dist=waveform_distance(f['local'],a['local'],f['times']);old=next(x for x in recorded['responses'][lev] if x['cell']==target)['contexts'][bg]['dose_dependence_normalized_rms'];numeric_match(dist,old,f'{seed} {target} {lev} {bg} {label} dose dependence')
  print('paired gates reconstructed',seed,label,flush=True)

fullm=[r for r in mechanical if r['horizon']==60]; fullp=[r for r in pulse_reports if r['horizon']==60]; fulle=[r for r in effects if r['horizon']==60]
summary=dict(logical_paths=len(path_reports),new_paths=sum(not x['reused'] for x in path_reports),reused_paths=sum(x['reused'] for x in path_reports),observations_per_path=401,window=[510.,570.],mechanical_pairs=len(fullm),response_pairs=len(fullp),effects=len(fulle),all_reconstructed_gates_pass=not failures,recorded_metric_agreement=not disagreements,mechanical_max={k:max(r['errors'][k] for r in fullm) for k in fullm[0]['errors']},max_normalized_response_error=max(r['normalized_response_error'] for r in fullp),max_relative_auc_error=max(r['relative_auc_error'] for r in fullp),max_recovery_time_error=max(v for r in fullp for v in r['recovery_time_errors'].values() if v is not None),max_effect_error=max(r['error'] for r in fulle),max_volume_error=max(r['audit']['max_volume_error'] for r in path_reports),min_radius=min(r['audit']['min_radius'] for r in path_reports),max_clipping=max(r['audit']['max_clipping'] for r in path_reports),max_boundary=max(r['sampled_boundary'] for r in path_reports),max_dilution_error=max(r['audit']['dilution_error_max'] for r in path_reports),max_dose_error=max(r['fine_dose']['relative_amount_error'] for r in fullp),minimum_sampled_contrast=min(r['contrast']['whole_window_min_log_sd'] for r in path_reports),minimum_late_contrast=min(r['contrast']['late_min_log_sd'] for r in path_reports),endpoint_CPU_volume_relative_discrepancy=max(r['endpoint_CPU_volume_relative_discrepancy'] for r in path_reports),max_transport_mass_residual=max(r['transport_mass_residual'] for r in path_reports))
limitations=['Every-step physical maxima are cumulative audits embedded in checkpoints and cross-checked against path result records; only stored 0.15-unit histories can be independently resampled without rerunning.','Checkpoint endpoint field volumes recomputed on CPU are diagnostic; GPU/CPU rounding differences do not replace actual recorded GPU volumes used in dose matching.','Only one timestep halving, three inherited mature histories, and 60-unit followup; no spatial convergence, stationarity, autonomous identity, challenge-time amount matching, or mechanical mediation claim.']
report=dict(summary=summary,failures=failures,recorded_metric_disagreements=disagreements,paths=path_reports,mechanical=mechanical,responses=pulse_reports,effects=effects,limitations=limitations)
OUT.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
print(json.dumps(summary,indent=2)); print('failures',failures);print('disagreements',disagreements[:10]);print('effects',json.dumps(fulle));print('REPORT',OUT)
