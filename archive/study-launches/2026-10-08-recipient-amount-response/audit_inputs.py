"""Literal saved-array preparation audit; no project scientific helpers."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
WATCHED={}
def digest(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
def pin(f,h=None):
 f=str(Path(f).resolve());v=digest(f)
 assert h is None or v==h,('hash',f)
 assert f not in WATCHED or WATCHED[f]==v,('changed',f)
 WATCHED[f]=v;return f
def read(f):return json.loads(Path(pin(f)).read_text())
def arrays(f):
 with np.load(pin(f),allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
def signature(z,u=None):return tuple((a.dtype.str,a.shape,a.tobytes()) for a in [z[k] for k in ('delta','volumes','ids')]+[z['initial'] if u is None else u])
def close(a,b):np.testing.assert_allclose(a,b,rtol=2e-13,atol=0)
def audit(root):
 root=Path(root).resolve();p=read(root/'protocol.json');old=read(Path(p['parent'])/'protocol.json')
 for field in ('source_sha256','input_sha256'):
  for f,h in p[field].items():pin(f,h)
 phase=p['phase'];n=96 if phase=='immediate' else 144
 assert (len(p['contexts']),len(p['trials']),len(p['factorials']))==(48,n,12)
 assert p['start']==(510. if phase=='immediate' else 570.) and p['duration']==60. and p['interval']==.15
 prepared={c['key']:arrays(c['source']) for c in p['contexts']};old_contexts={c['key']:c for c in old['contexts']}
 amounts=[];allocation_separations=[];sources=[];groups={};physical=[]
 if phase=='immediate':
  donors={c['key']:arrays(c['source']) for c in old['contexts']}
  for f in p['factorials']:
   parent=next(x for x in old['factorials'] if all(x[k]==f[k] for k in ('seed','cell','matching_reference')))
   for P in ('S','R'):
    first,second=[donors[parent['contexts'][C+'_'+P]] for C in ('S','R')]
    i=first['ids'].tolist().index(f['cell']);v=first['volumes'];others=np.arange(len(v))!=i
    for k in ('delta','volumes','ids'):np.testing.assert_array_equal(first[k],second[k])
    np.testing.assert_array_equal(first['initial'][:,i],second['initial'][:,i])
    budgets=[z['initial'][:,others]@v[others] for z in (first,second)];common=np.sqrt(budgets[0])*np.sqrt(budgets[1])
    for C,z,before in zip(('S','R'),(first,second),budgets):
     key=f['contexts'][C+'_'+P];c=next(x for x in p['contexts'] if x['key']==key);new=prepared[key];u=z['initial'];expected=u.copy();expected[:,others]*=(common/before)[:,None]
     assert c['parent_context']==parent['contexts'][C+'_'+P]
     for k in ('delta','volumes','ids'):np.testing.assert_array_equal(new[k],z[k])
     np.testing.assert_array_equal(new['initial'],expected);np.testing.assert_array_equal(new['initial'][:,i],u[:,i])
     after=expected[:,others]@v[others];fb=u[:,others]*v[others]/before[:,None];fa=expected[:,others]*v[others]/after[:,None]
     record=c['amount_matching']
     vals=dict(target_pair=u[:,i],surrounding_indices=np.flatnonzero(others),pre_surrounding_species_amounts=before,common_surrounding_species_amounts=common,scale_factors=common/before,post_surrounding_species_amounts=after,pre_per_cell_species_amounts=u*v,post_per_cell_species_amounts=expected*v,per_cell_species_amount_change=(expected-u)*v,pre_global_species_amounts=u@v,post_global_species_amounts=expected@v,pre_surrounding_allocation_fractions=fb,post_surrounding_allocation_fractions=fa)
     for name,value in vals.items():np.testing.assert_array_equal(record[name],value)
     err=float(np.max(abs(after/common-1.)));alloc=float(np.max(abs(fa/fb-1.)))
     assert max(err,alloc)<=2e-13 and err==record['maximum_relative_amount_error'] and alloc==record['maximum_relative_allocation_error']
     amounts.append(dict(context=key,amount_relative_error=err,allocation_relative_error=alloc,scale_factors=(common/before).tolist()))
    a,b=[prepared[f['contexts'][C+'_'+P]]['initial'] for C in ('S','R')]
    close(a@v,b@v);close(a[:,others]@v[others],b[:,others]@v[others]);np.testing.assert_array_equal(a[:,i],b[:,i])
    allocation_separations.append(dict(seed=f['seed'],cell=f['cell'],matching_reference=f['matching_reference'],physics=P,maximum_surrounding_log_separation=float(np.max(abs(np.log(a[:,others]/b[:,others]))))))
 else:
  times=np.arange(401)*.15;paths={}
  for label in ('coarse','fine','diagnostic'):
   z=arrays(Path(p['parent'])/'production'/label/'paths.npz');np.testing.assert_array_equal(z['times'],times)
   paths[label]={str(k):a for k,a in zip(z['trial_keys'],z['trajectory'].transpose(1,0,2,3))}
  for c in p['contexts']:
   prior=old_contexts[c['parent_context']];z=arrays(prior['source']);ctrl=next(t for t in old['trials'] if t['context']==prior['key'] and t['target'] is None)
   endpoint=paths['fine'][ctrl['key']][-1];new=prepared[c['key']]
   for k in ('delta','volumes','ids'):np.testing.assert_array_equal(new[k],z[k])
   np.testing.assert_array_equal(new['initial'],endpoint)
   refs=arrays(Path(p['parent'])/'references'/ctrl['key']/'paths.npz');np.testing.assert_array_equal(refs['times'],times)
   errors={label:float(np.max(abs(np.log(endpoint/other)))) for label,other in [(l,paths[l][ctrl['key']][-1]) for l in ('coarse','diagnostic')]+[(l,refs[l][-1]) for l in ('dop853','radau')]}
   assert max(errors.values())<=1e-8
   drift=float(np.max(abs(np.log(paths['fine'][ctrl['key']][times>=36]/endpoint))))
   a,b=endpoint;d=new['delta'];rhs=np.stack([a*a/b-a+p['da']*(d@a),a*a-p['beta']*b+p['db']*(d@b)])
   rate=float(np.max(abs(rhs/endpoint)));settled=drift<=.001 and rate<=.0001
   record=c['settling'];assert record['qualified'] and record['settled']==settled
   close(record['maximum_reference_endpoint_log_error'],max(errors.values()))
   close(record['late_max_log_drift'],drift);close(record['maximum_relative_rhs'],rate)
   close(c['challenge_global_species_amounts'],endpoint@new['volumes']);close(c['challenge_per_cell_species_amounts'],endpoint*new['volumes'])
   sources.append(dict(context=c['key'],errors=errors,late_max_log_drift=drift,maximum_relative_rhs=rate,settled=bool(settled)))
 starts={};controls={};unique=set();label_unique=set();pulse_errors=[]
 for c in p['contexts']:
  z=prepared[c['key']];d,v,u,ids=[z[k] for k in ('delta','volumes','initial','ids')]
  assert d.dtype==v.dtype==u.dtype==np.float64 and d.shape==(16,16) and u.shape==(2,16) and len(set(ids.tolist()))==16
  assert np.isfinite(u).all() and np.all(u>0) and np.isfinite(v).all() and np.all(v>0)
  off=d.copy();np.fill_diagonal(off,0);assert np.all(off>=0) and np.max(abs(d.sum(axis=1)))<1e-11
  assert np.max(abs(v@d))<1e-12 and np.max(abs(v[:,None]*d-(v[:,None]*d).T))<1e-12
  physical.append(dict(context=c['key'],capacity_conservation_error=float(np.max(abs(v@d)))))
 for t in p['trials']:
  z=prepared[t['context']];u=arrays(t['source'])['initial'];expected=z['initial'].copy();factor=1.
  if t['target'] is None:controls[t['context']]=t['key']
  else:
   i=z['ids'].tolist().index(t['target']);f=next(f for f in p['factorials'] if all(f[k]==t[k] for k in ('seed','matching_reference')) and f['cell']==t['target'])
   C=next(C for C in ('S','R') if t['context'] in (f['contexts'][C+'_S'],f['contexts'][C+'_R']));P=next(P for P in ('S','R') if t['context']==f['contexts'][C+'_'+P])
   values=[prepared[f['contexts'][cc+'_'+P]]['initial'][0,i] for cc in ('S','R')]
   requested=.1*float(z['initial'][0,i]*z['volumes'][i]) if t['dose']=='fractional' else .1*min(values)*float(z['volumes'][i])
   factor=.9 if t['dose']=='fractional' else 1-requested/float(z['initial'][0,i]*z['volumes'][i]);expected[0,i]*=factor
   np.testing.assert_array_equal(expected,u);assert factor==t['factor']
   close(t['requested_removed_amount'],requested);removed=float((z['initial'][0,i]-u[0,i])*z['volumes'][i]);close(t['removed_activator_amount'],removed)
   err=abs(removed/requested-1);assert err<=2e-13;pulse_errors.append(err)
  np.testing.assert_array_equal(u,expected);assert t['factor']==factor and t['initial_log_amplitude']==abs(np.log(factor))
  starts[t['key']]=u;unique.add(signature(z,u));label_unique.add((t['context'],signature(z,u)))
  sig=tuple((z[k].dtype.str,z[k].shape,z[k].tobytes()) for k in ('delta','volumes','ids'));groups.setdefault(sig,[]).append(t['key'])
 assert len(controls)==48 and len(pulse_errors)==n-48 and len(unique)==p['unique_physical_initial_value_problems'] and len(label_unique)==p['unique_initial_value_problems']
 assert len(groups)==9 and max(map(len,groups.values()))==p['batching']['width']
 assert 9*p['batching']['width']-n==p['batching']['padding_paths']
 for f in p['factorials']:
  for P in ('S','R'):
   a,b=[prepared[f['contexts'][C+'_'+P]] for C in ('S','R')];i=a['ids'].tolist().index(f['cell'])
   for k in ('delta','volumes','ids'):np.testing.assert_array_equal(a[k],b[k])
   if phase=='immediate':np.testing.assert_array_equal(a['initial'][:,i],b['initial'][:,i]);close(a['initial']@a['volumes'],b['initial']@b['volumes'])
 for f,h in list(WATCHED.items()):assert digest(f)==h,f
 return dict(passed=True,phase=phase,protocol_sha256=digest(root/'protocol.json'),contexts=48,logical_cases=n,unique_initial_value_problems=len(label_unique),unique_physical_initial_value_problems=len(unique),controls=48,pulses=n-48,batching=p['batching'],amount_matching=amounts,allocation_separations=allocation_separations,waiting_sources=sources,physical_checks=physical,maximum_pulse_relative_amount_error=max(pulse_errors),inputs_unchanged=True,inputs_sha256=WATCHED,script_sha256=digest(__file__))
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--root',required=True);parser.add_argument('--output',required=True);args=parser.parse_args();out=Path(args.output)
 assert not out.exists();result=audit(args.root);out.write_text(json.dumps(result,indent=2)+'\n');print({k:result[k] for k in ('passed','phase','contexts','logical_cases','unique_initial_value_problems','unique_physical_initial_value_problems')})
