"""Gated moving-geometry formation and retention after conservative exchange."""
import argparse
import json
from pathlib import Path
import shutil
import time
import numpy as np
from . import cell_response_moving as moving
from .attribute_exchange import exchange
from .attribute_persistence import rhs, distance
from .cell_response_exchange import checked_solve
from .feedback_endpoint_bistability import chemical_jacobian
from .feedback_long import digest
from .resolution import write_json

FAMILIES=('fresh_exchange','relaxed_exchange')


def gate_status(dependency):
    dependency=Path(dependency)
    status=json.loads((dependency/'status.json').read_text())
    if status['state']=='failed':return 'blocked'
    report=dependency/'refinement.json'
    if status['state']=='completed' and report.exists():
        return 'ready' if json.loads(report.read_text())['passed'] else 'blocked'
    return 'waiting'


def prepare(root,dependency):
    root,dependency=Path(root).resolve(),Path(dependency).resolve()
    if root.exists():raise FileExistsError(root)
    previous=json.loads((dependency/'protocol.json').read_text())
    for f,h in {**previous['source_sha256'],**previous['input_sha256']}.items():
        if digest(f)!=h:raise ValueError('Changed refinement input: '+f)
    evidence=Path('outputs/cell-response-exchange').resolve()
    ep=json.loads((evidence/'protocol.json').read_text())
    for f,h in {**ep['source_sha256'],**ep['input_sha256']}.items():
        if digest(f)!=h:raise ValueError('Changed exchange input: '+f)
    graph=next(g for g in ep['graphs'] if g['key']=='seed-7_switch_on')
    folder=evidence/graph['key'];result=json.loads((folder/'conservative.json').read_text())
    if not result['settled'] or digest(folder/'conservative.npz')!=result['trajectory_sha256']:raise ValueError('Requires validated relaxed exchange state')
    with np.load(dependency/'initial_states.npz') as d:base=d['pattern'].copy();ids=d['ids'];m=d['masses']
    with np.load(folder/'conservative.npz') as d:
        if not np.array_equal(ids,d['ids']) or not np.allclose(m,d['masses'],rtol=1e-12,atol=1e-12):raise ValueError('Source mismatch')
        relaxed=d['relaxation'][-1].copy();delta=d['delta'].copy()
    pair=graph['pair'];fresh=exchange(base,m,*pair)
    with np.load(graph['source']) as d:
        if not np.allclose(base,d['trajectories'][0,-1],rtol=1e-12,atol=1e-12):raise ValueError('Initial chemistry mismatch')
    root.mkdir(parents=True);shutil.copy2(previous['checkpoint'],root/'source.npz')
    np.savez_compressed(root/'initial_states.npz',fresh_exchange=fresh,relaxed_exchange=relaxed,baseline=base,ids=ids,masses=m)
    times=np.arange(401)*.15;fun=rhs(delta,graph['beta'],graph['da'],graph['db'])
    jac=lambda t,y:chemical_jacobian(y.reshape(2,-1),delta,graph['beta'],graph['da'],graph['db'])
    frozen={};errors={}
    for key,x in dict(baseline=base,fresh_exchange=fresh,relaxed_exchange=relaxed).items():
        frozen[key],errors[key]=checked_solve(fun,jac,x,times)
    np.savez_compressed(root/'frozen_reference.npz',times=times,**frozen)
    files=[root/'source.npz',root/'initial_states.npz',root/'frozen_reference.npz',dependency/'protocol.json',dependency/'pattern_control/history.json',dependency/'pattern_control/result.json',evidence/'protocol.json',folder/'conservative.json',folder/'conservative.npz',folder/'comparison.json']
    p=dict(checkpoint=str(root/'source.npz'),start=previous['start'],dt=previous['dt'],duration=60.,interval=.15,
        dependency=str(dependency),dependency_protocol_sha256=digest(dependency/'protocol.json'),pair=pair,selected_ids=[int(ids[i]) for i in pair],
        jobs=[dict(family=key,target=None,factor=1.) for key in FAMILIES],
        frozen_solver_errors=errors,
        criteria=dict(volume_max=.05,radius_min=4.,clipping_max=0.,boundary_max=.01,late_start=36.,pair_return_ratio=.1),
        scope='Seed 7 full-coupling mature t=150 geometry. Conservative fresh exchange tests formation opportunity; pre-relaxed exchanged chemistry tests retention after transplantation onto the same starting geometry and polarity. Both compared to validated unexchanged moving control. No cell displacement, new cleavage, pulse-response transfer, or autonomous identity claim.',
        analysis='Late maximum volume-weighted pair log distances normalized by initial fresh-exchange separation, comparing moving unexchanged control and its conservatively exchanged reference; additionally compare the pre-relaxed pattern. Report continuous distances, pair ordering, overall contrast and geometric differences. No static-equilibrium claim for moving states.',
        time_refinement='Uses dt=.00375; dependency validates small-pulse responses, not the larger exchange intervention. Exchange-specific refinement remains necessary.',
        input_sha256={str(f):digest(f) for f in files},
        source_sha256={str(f.resolve()):digest(f) for f in [Path(__file__),*[Path(__file__).with_name(x) for x in ('cell_response_moving.py','cell_response_exchange.py','attribute_exchange.py','attribute_persistence.py','feedback_endpoint_bistability.py','feedback_long.py','attribute_development.py','model.py','transport.py','signaling.py','polarity.py','resolution.py')]]})
    write_json(root/'protocol.json',p);write_json(root/'status.json',dict(state='queued',completed=0,total=2,reason='Awaiting completed passing moving-response timestep refinement'))


def verify(p):
    for f,h in {**p['source_sha256'],**p['input_sha256']}.items():
        if digest(f)!=h:raise ValueError('Frozen input changed: '+f)


def compare(root,p):
    baseline=json.loads((Path(p['dependency'])/'pattern_control/history.json').read_text())
    with np.load(root/'initial_states.npz') as d:m=d['masses'];ids=d['ids'];initial={k:d[k] for k in ('baseline',*FAMILIES)}
    pair=p['pair'];times=np.array([x['elapsed'] for x in baseline]);base=np.array([x['chemistry'] for x in baseline]);late=times>=p['criteria']['late_start']
    sep=float(distance(initial['fresh_exchange'][:,pair],initial['baseline'][:,pair],m[pair]));donor=np.array([exchange(x,m,*pair) for x in base])
    rows=[]
    with np.load(root/'frozen_reference.npz') as frozen:
        for job in p['jobs']:
            folder=root/moving.name(job);r=json.loads((folder/'result.json').read_text())
            if not r['quality_pass'] or r['protocol_sha256']!=digest(root/'protocol.json'):raise ValueError('Invalid completed result')
            h=json.loads((folder/'history.json').read_text())
            if len(h)!=len(baseline) or not np.allclose([x['elapsed'] for x in h],times,rtol=0,atol=1e-9) or any(not np.array_equal(x['ids'],ids) for x in h+baseline):raise ValueError('Unaligned histories')
            path=np.array([x['chemistry'] for x in h]);end=h[-1];control=baseline[-1]
            destination=float(distance(path[:,:,pair],base[:,:,pair],m[pair])[late].max()/sep)
            transferred=float(distance(path[:,:,pair],donor[:,:,pair],m[pair])[late].max()/sep)
            rows.append(dict(family=job['family'],late_destination_ratio=destination,late_transferred_ratio=transferred,
                destination_like=destination<.1,transferred_like=transferred<.1,
                late_pair_order_reversed=bool(np.all(path[late,0,pair[0]]>path[late,0,pair[1]])),
                final_global_distance_from_unexchanged=float(distance(path[-1],base[-1],m)),
                final_global_distance_from_frozen=float(distance(path[-1],frozen[job['family']][-1],m)),
                final_global_distance_from_pre_relaxed=float(distance(path[-1],initial['relaxed_exchange'],m)),
                final_log_activator_sd=end['log_activator_sd'],
                final_relative_operator_difference=float(np.linalg.norm(np.array(end['delta'])-control['delta'])/np.linalg.norm(control['delta'])),
                final_relative_axis_ratio_difference=float(abs(end['axis_ratio']/control['axis_ratio']-1))))
    result=dict(trials=rows,scope=p['scope'],limitation=p['time_refinement'])
    write_json(root/'comparison.json',result)


def run(root,wait=False):
    root=Path(root);p=json.loads((root/'protocol.json').read_text());verify(p)
    while True:
        state=gate_status(p['dependency'])
        if state=='blocked':
            write_json(root/'status.json',dict(state='blocked',completed=0,total=2,reason='Prerequisite numerical validation failed'));return
        if state=='ready':break
        if not wait:raise RuntimeError('Prerequisite still running; use --wait')
        time.sleep(30)
    verify(p)
    report=Path(p['dependency'])/'refinement.json'
    write_json(root/'dependency_check.json',dict(refinement_sha256=digest(report),refinement=json.loads(report.read_text())))
    done=0
    try:
        for job in p['jobs']:
            write_json(root/'status.json',dict(state='running',completed=done,total=2,current=job['family'],workers=1))
            moving.worker((str(root),job));done+=1
        compare(root,p);write_json(root/'status.json',dict(state='completed',completed=done,total=2))
    except Exception as exc:
        write_json(root/'status.json',dict(state='failed',completed=done,total=2,error=str(exc)));raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('command',choices=('prepare','run'));parser.add_argument('--output',type=Path,default=Path('outputs/cell-exchange-moving'));parser.add_argument('--dependency',type=Path,default=Path('outputs/cell-response-moving-refined'));parser.add_argument('--wait',action='store_true')
    args=parser.parse_args()
    if args.command=='prepare':prepare(args.output,args.dependency)
    else:run(args.output,args.wait)
