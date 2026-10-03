"""Frozen chemical-state exchange followed by matched local response assays."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp, trapezoid
from .attribute_exchange import exchange, classify
from .attribute_persistence import rhs, distance
from .cell_response import TIMES, pulse, response_metrics
from .feedback_endpoint_bistability import chemical_jacobian
from .feedback_long import digest
from .resolution import write_json

ARMS=('untouched','sham','exact','conservative')
RELAX_TIMES=np.arange(0.,1201.,2.)


def transplant(initial,masses,pair,arm):
    x=np.array(initial,float,copy=True);i,j=pair
    if arm=='untouched':return x
    if arm=='sham':
        x[:,[i,j]]=x[:,[i,j]].copy()
        return x
    if arm=='exact':x[:,[i,j]]=x[:,[j,i]];return x
    if arm=='conservative':return exchange(x,masses,i,j)
    raise ValueError(arm)


def waveform_distance(first,second,times):
    """Time RMS of two signed, pulse-normalized local species responses."""
    first,second=np.asarray(first),np.asarray(second)
    if first.shape!=second.shape or first.shape!=(len(times),2):raise ValueError('Misaligned response')
    return float(np.sqrt(trapezoid(np.mean((first-second)**2,axis=1),times)/(times[-1]-times[0])))


def checked_solve(fun,jac,initial,times):
    paths=[]
    for method in ('DOP853','Radau'):
        kwargs={'jac':jac} if method=='Radau' else {}
        sol=solve_ivp(fun,(times[0],times[-1]),initial.ravel(),method=method,t_eval=times,rtol=1e-12,atol=1e-14,**kwargs)
        if not sol.success or not np.isfinite(sol.y).all() or np.any(sol.y<=0):raise RuntimeError('Solver positivity/finiteness failure')
        paths.append(sol.y.T.reshape(len(times),*initial.shape))
    error=float(abs(np.log(paths[0]/paths[1])).max())
    if error>=1e-5:raise RuntimeError(f'Independent solver disagreement {error}')
    return paths[0],error


def prepare(root):
    root=Path(root).resolve()
    if root.exists():raise FileExistsError(root)
    inputs={};graphs=[]
    for seed in (7,8,9):
        source=Path('outputs/feedback-endpoint-bistability'+('' if seed==7 else f'-seed-{seed}'))
        prior=json.loads((source/'protocol.json').read_text())
        for f,h in prior['sha256'].items():
            if digest(f)!=h:raise ValueError('Changed validated input: '+f)
        results=json.loads((source/'results.json').read_text())
        checkpoint_root=Path('outputs/feedback-survival' if seed==7 else f'outputs/feedback-survival-validation/seed-{seed}/survival')
        from .attribute_development import AttributeSimulation
        for branch in ('switch_on','keep_off'):
            if not results[branch]['chemical_bistability_supported']:raise ValueError('Unvalidated endpoint')
            file=(source/f'{branch}.npz').resolve()
            with np.load(file) as d:
                initial=d['trajectories'][0,-1];ids=d['ids'];delta=d['delta'];m=d['volumes']
                order=sorted(range(len(ids)),key=lambda i:(initial[0,i],int(ids[i])))
                pair=[order[0],order[-1]]
            checkpoint=checkpoint_root/branch/'latest_state.npz';sim=AttributeSimulation.restore(checkpoint);g=sim.signaling_graph()
            if not np.array_equal(sim.ids,ids) or not np.allclose(g.delta,delta,rtol=1e-12,atol=1e-12) or not np.allclose(g.masses,m):raise ValueError('Source graph mismatch')
            c=sim.config
            graphs.append(dict(key=f'seed-{seed}_{branch}',seed=seed,branch=branch,source=str(file),pair=pair,ids=[int(ids[i]) for i in pair],beta=c.signal_beta,da=c.signal_da,db=c.signal_dh))
            inputs[str(checkpoint.resolve())]=digest(checkpoint);inputs[str(file)]=digest(file)
        for filename in ('protocol.json','results.json'):inputs[str((source/filename).resolve())]=digest(source/filename)
    root.mkdir(parents=True)
    p=dict(graphs=graphs,arms=ARMS,relaxation_times=RELAX_TIMES.tolist(),response_times=TIMES.tolist(),factors=[.9,1.1],
        criteria=dict(max_solver_log_error=1e-5,stationarity_rhs=1e-6,stable_growth_max=0.,pair_informative_min=.1,pair_return_ratio=.1,response_reference_separation_min=.01),
        selection='Min/max initial patterned activator, ties by cell ID; selected before exchange results.',
        design='Six frozen endpoint graphs, three histories. Four arms per graph. Relax 1200; pulse each selected cell +/-10% for 240 against its own unperturbed continuation. External pulse doses are fractional and differ in absolute amount.',
        interpretation='Signed two-species target log-response waveforms normalized by absolute log pulse. Compare to untouched donor and destination waveforms at same pulse factor. Continuous distances, descriptive nearest reference only if baseline separation exceeds .01. Nearer donor is not equivalence, autonomy, or a discrete identity. Exchanges can reorganize the whole interacting network.',
        input_sha256=inputs,source_sha256={str(f.resolve()):digest(f) for f in (Path(__file__),*[Path(__file__).with_name(x) for x in ('attribute_exchange.py','attribute_persistence.py','cell_response.py','feedback_endpoint_bistability.py','feedback_long.py','resolution.py')])})
    write_json(root/'protocol.json',p);write_json(root/'status.json',dict(state='prepared',completed=0,total=24))


def run(root):
    root=Path(root);p=json.loads((root/'protocol.json').read_text());ph=digest(root/'protocol.json')
    for f,h in {**p['input_sha256'],**p['source_sha256']}.items():
        if digest(f)!=h:raise ValueError('Frozen input changed: '+f)
    done=0;all_results={}
    try:
        for graph in p['graphs']:
            key=graph['key'];pair=graph['pair'];folder=root/key;folder.mkdir(exist_ok=True)
            with np.load(graph['source']) as d:initial=d['trajectories'][0,-1].copy();delta=d['delta'];m=d['volumes'];ids=d['ids']
            fun=rhs(delta,graph['beta'],graph['da'],graph['db'])
            jac=lambda t,y:chemical_jacobian(y.reshape(2,-1),delta,graph['beta'],graph['da'],graph['db'])
            paths={};records={}
            for arm in p['arms']:
                result_file=folder/f'{arm}.json';path_file=folder/f'{arm}.npz'
                write_json(root/'status.json',dict(state='running',completed=done,total=24,current=f'{key}/{arm}'))
                if result_file.exists():
                    record=json.loads(result_file.read_text())
                    if record['protocol_sha256']!=ph or digest(path_file)!=record['trajectory_sha256']:raise ValueError('Saved result mismatch')
                    with np.load(path_file) as d:paths[arm]={k:d[k] for k in d.files}
                    records[arm]=record;done+=1;continue
                start=transplant(initial,m,pair,arm)
                relax,error=checked_solve(fun,jac,start,RELAX_TIMES)
                end=relax[-1];residual=float(abs(fun(1200,end.ravel())).max());growth=float(np.linalg.eigvals(jac(1200,end.ravel())).real.max())
                control,err=checked_solve(fun,jac,end,TIMES);errors=[error,err];trials=[];responses=[];pulse_paths=[]
                for cell in pair:
                    for factor in p['factors']:
                        path,err=checked_solve(fun,jac,pulse(end,cell,factor),TIMES);errors.append(err)
                        wave=np.log(path[:,:,cell]/control[:,:,cell])/abs(np.log(factor))
                        responses.append(wave);pulse_paths.append(path)
                        trials.append(dict(cell=int(ids[cell]),index=cell,factor=factor,activator_amount_added=float((factor-1)*end[0,cell]*m[cell]),**response_metrics(path,control,TIMES,m,cell,factor)))
                arrays=dict(relaxation=relax,control=control,responses=np.array(responses),pulse_paths=np.array(pulse_paths),ids=ids,masses=m,delta=delta,relaxation_times=RELAX_TIMES,response_times=TIMES)
                np.savez_compressed(path_file,**arrays)
                record=dict(protocol_sha256=ph,trajectory_sha256=digest(path_file),max_solver_log_error=max(errors),final_rhs_max=residual,final_growth=growth,settled=bool(residual<1e-6 and growth<0),initial_amount_change=(start@m-initial@m).tolist(),trials=trials)
                write_json(result_file,record);paths[arm]=arrays;records[arm]=record;done+=1
                print(f'{key}/{arm}: settled={record["settled"]}, max solver error={max(errors):.3g}',flush=True)
            comparisons=[];baseline=paths['untouched'];late=RELAX_TIMES>=1000
            for arm in ('exact','conservative'):
                transplant_ref=np.array([transplant(x,m,pair,arm) for x in baseline['relaxation']]);path=paths[arm]['relaxation']
                sep=float(distance(path[0][:,pair],initial[:,pair],m[pair]))
                dest=float(distance(path[:,:,pair],baseline['relaxation'][:,:,pair],m[pair])[late].max()/max(sep,1e-30))
                donor=float(distance(path[:,:,pair],transplant_ref[:,:,pair],m[pair])[late].max()/max(sep,1e-30))
                response_rows=[]
                for index,trial in enumerate(records[arm]['trials']):
                    donor_index=(index+2)%4;wave=paths[arm]['responses'][index];destination_wave=baseline['responses'][index];donor_wave=baseline['responses'][donor_index]
                    separation=waveform_distance(destination_wave,donor_wave,TIMES);dd=waveform_distance(wave,destination_wave,TIMES);ds=waveform_distance(wave,donor_wave,TIMES)
                    eligible=records[arm]['settled'] and records['untouched']['settled'] and separation>.01
                    response_rows.append(dict(cell=trial['cell'],factor=trial['factor'],baseline_separation=separation,destination_distance=dd,donor_distance=ds,donor_minus_destination_distance=ds-dd,nearest_reference=('donor' if ds<dd else 'destination') if eligible else 'unresolved'))
                comparisons.append(dict(arm=arm,state_outcome=classify(dest,donor,records[arm]['settled'],sep>.1),destination_ratio=dest,transferred_ratio=donor,global_destination_distance=float(distance(path[-1],baseline['relaxation'][-1],m)),responses=response_rows))
            sham_error=float(abs(np.log(paths['sham']['relaxation']/baseline['relaxation'])).max())
            sham_response_error=float(abs(paths['sham']['responses']-baseline['responses']).max())
            if sham_error>1e-10 or sham_response_error>1e-10:raise RuntimeError('Sham differs from untouched')
            summary=dict(graph=graph,sham_max_log_error=sham_error,sham_response_error=sham_response_error,comparisons=comparisons)
            write_json(folder/'comparison.json',summary);all_results[key]=summary
            write_json(root/'results.json',all_results)
        write_json(root/'status.json',dict(state='completed',completed=done,total=24))
    except Exception as exc:
        write_json(root/'status.json',dict(state='failed',completed=done,total=24,error=str(exc)));raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('command',choices=('prepare','run'));parser.add_argument('--output',type=Path,default=Path('outputs/cell-response-exchange'))
    a=parser.parse_args()
    if a.command=='prepare':prepare(a.output)
    else:run(a.output)
