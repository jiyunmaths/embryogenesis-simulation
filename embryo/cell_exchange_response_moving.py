"""Matched moving pulse responses after moving chemical-state exchange."""
import argparse
import json
from pathlib import Path
import shutil
import time
import numpy as np
from . import cell_response_moving as moving
from .attribute_development import AttributeSimulation
from .cell_response import response_metrics
from .cell_response_exchange import waveform_distance
from .feedback_long import digest
from .resolution import write_json

BACKGROUNDS=('unexchanged','fresh_exchange','relaxed_exchange')


def local_response(path,control,cell,factor):
    path,control=np.asarray(path),np.asarray(control)
    if path.shape!=control.shape or path.ndim!=3 or path.shape[1]!=2 or not np.isfinite(path).all() or not np.isfinite(control).all() or np.any(path<=0) or np.any(control<=0):raise ValueError('Unaligned or invalid chemistry')
    if factor<=0 or factor==1:raise ValueError('A nonzero positive pulse is required')
    return np.log(path[:,:,cell]/control[:,:,cell])/abs(np.log(factor))


def donor_id(cell,targets):
    if len(targets)!=2 or len(set(targets))!=2 or cell not in targets:raise ValueError('Expected a distinct exchange pair')
    return targets[1-targets.index(cell)]


def verify(p):
    for f,h in {**p.get('source_sha256',{}),**p.get('input_sha256',{})}.items():
        if digest(f)!=h:raise ValueError('Frozen input changed: '+f)


def prepare(root):
    root=Path(root).resolve()
    if root.exists():raise FileExistsError(root)
    exchange=Path('outputs/cell-exchange-moving').resolve();refined=Path('outputs/cell-response-moving-refined').resolve()
    ep=json.loads((exchange/'protocol.json').read_text());verify(ep)
    rp=json.loads((refined/'protocol.json').read_text());verify(rp)
    if not json.loads((refined/'refinement.json').read_text())['passed']:raise ValueError('Requires accepted preceding response refinement')
    sources={'unexchanged':dict(root=str(refined),folder='pattern_control')}
    sources.update({key:dict(root=str(exchange),folder=key+'_control') for key in BACKGROUNDS[1:]})
    root.mkdir(parents=True)
    files=[exchange/'protocol.json',refined/'protocol.json',refined/'refinement.json']
    code=[Path(__file__),*[Path(__file__).with_name(f) for f in ('cell_response_moving.py','cell_response.py','cell_response_exchange.py','attribute_development.py','attribute_persistence.py','feedback_endpoint_bistability.py','feedback_long.py','model.py','polarity.py','transport.py','signaling.py','resolution.py')]]
    p=dict(backgrounds=BACKGROUNDS,sources=sources,targets=ep['selected_ids'],start=ep['start']+ep['duration'],dt=ep['dt'],duration=60.,interval=.15,factors=[.9,1.1],
        criteria=dict(volume_max=.05,radius_min=4.,clipping_max=0.,boundary_max=.01,response_reference_separation_min=.01,recovery_ratio=.1,recovery_hold=24.),
        scope='Seed 7, same-age t=210 endpoints from unexchanged, fresh-exchange, and pre-relaxed-exchange moving trajectories. Preserve each native geometry, chemistry, polarity and random streams. Five moving continuations per endpoint: matched unperturbed control and +/-10% activator pulses in each original exchange cell. No re-equilibration and no new chemical transplant.',
        interpretation='Compare signed two-species target log-response waveforms relative to each background\'s own moving control, normalized by absolute log pulse, over 60 units. Untouched same-age donor and destination responses are references; nearest donor is descriptive, not equivalence or autonomy. Starting geometries and collective chemical states differ intentionally. Three backgrounds are interventions in one history, not replicates.',
        limitations='No new pulse-specific timestep or spatial refinement; preceding small-pulse check at t=150 does not establish convergence at these exchange endpoints. No full moving-system stability or autonomous cell-type claim. All compartments remain well mixed and capped at 16 cells.',
        scheduling='One worker; each source must have completed with passing quality. Start available same-age endpoints now and wait for unfinished source branches; do not use partial checkpoints.',
        source_sha256={str(f.resolve()):digest(f) for f in code},input_sha256={str(f):digest(f) for f in files})
    write_json(root/'protocol.json',p);write_json(root/'status.json',dict(state='prepared',completed=0,total=15))


def source_ready(source):
    root=Path(source['root']);folder=root/source['folder']
    if (folder/'result.json').exists():
        result=json.loads((folder/'result.json').read_text())
        if not result['quality_pass'] or result['protocol_sha256']!=digest(root/'protocol.json'):raise ValueError('Invalid source result')
        return True
    for f in (folder/'status.json',root/'status.json'):
        if f.exists() and json.loads(f.read_text())['state'] in ('failed','blocked'):raise RuntimeError('Source study failed')
    return False


def materialize(root,p,key):
    child=root/key
    if (child/'protocol.json').exists():verify(json.loads((child/'protocol.json').read_text()));return
    source=p['sources'][key];folder=Path(source['root'])/source['folder']
    if not source_ready(source):raise ValueError('Source incomplete')
    checkpoint=folder/'latest_state.npz';sim=AttributeSimulation.restore(checkpoint)
    if abs(sim.time-p['start'])>1e-9 or sim.config.dt!=p['dt'] or sim.divisions or len(sim.ids)!=16:raise ValueError('Source age/dt/cleavage mismatch')
    if any(cell not in sim.ids for cell in p['targets']):raise ValueError('Target cells missing')
    history=json.loads((folder/'history.json').read_text())
    if not np.array_equal(history[-1]['ids'],sim.ids) or not np.allclose(history[-1]['chemistry'],[sim.activator,sim.inhibitor],atol=0,rtol=1e-12):raise ValueError('Checkpoint does not match source endpoint')
    if child.exists():raise ValueError('Incomplete prior preparation; preserve and inspect '+str(child))
    child.mkdir();shutil.copy2(checkpoint,child/'source.npz')
    np.savez_compressed(child/'initial_states.npz',**{key:np.array([sim.activator,sim.inhibitor])},ids=sim.ids,masses=sim.volumes())
    jobs=[dict(family=key,target=None,factor=1.)]+[dict(family=key,target=cell,factor=factor) for cell in p['targets'] for factor in p['factors']]
    files=[root/'protocol.json',checkpoint,folder/'result.json',folder/'history.json',child/'source.npz',child/'initial_states.npz']
    cp=dict(checkpoint=str(child/'source.npz'),start=p['start'],dt=p['dt'],duration=p['duration'],interval=p['interval'],jobs=jobs,source_sha256=p['source_sha256'],input_sha256={str(f):digest(f) for f in files})
    write_json(child/'protocol.json',cp)


def assess(root,p):
    datasets={};summaries={}
    times=np.arange(round(p['duration']/p['interval'])+1)*p['interval']
    for key in p['backgrounds']:
        child=root/key;cp=json.loads((child/'protocol.json').read_text());verify(cp);ph=digest(child/'protocol.json');paths={}
        with np.load(child/'initial_states.npz') as d:ids=d['ids'];m=d['masses']
        for job in cp['jobs']:
            folder=child/moving.name(job);result=json.loads((folder/'result.json').read_text())
            if not result['quality_pass'] or result['protocol_sha256']!=ph:raise ValueError('Invalid pulse result')
            h=json.loads((folder/'history.json').read_text())
            if len(h)!=len(times) or not np.allclose([x['elapsed'] for x in h],times,rtol=0,atol=1e-9) or any(not np.array_equal(x['ids'],ids) for x in h):raise ValueError('Misaligned observations')
            paths[(job['target'],job['factor'])]=np.array([x['chemistry'] for x in h])
        control=paths[(None,1.)];waves={};rows=[]
        for cell in p['targets']:
            index=list(ids).index(cell)
            for factor in p['factors']:
                path=paths[(cell,factor)];waves[(cell,factor)]=local_response(path,control,index,factor)
                rows.append(dict(cell=cell,factor=factor,injected_activator_amount=float((factor-1)*control[0,0,index]*m[index]),**response_metrics(path,control,times,m,index,factor)))
        summaries[key]=rows;datasets[key]=waves
    rows=[]
    for key in BACKGROUNDS[1:]:
        for cell in p['targets']:
            for factor in p['factors']:
                donor=donor_id(cell,p['targets']);wave=datasets[key][(cell,factor)]
                destination_wave=datasets['unexchanged'][(cell,factor)];donor_wave=datasets['unexchanged'][(donor,factor)]
                separation=waveform_distance(destination_wave,donor_wave,times);dd=waveform_distance(wave,destination_wave,times);ds=waveform_distance(wave,donor_wave,times)
                rows.append(dict(background=key,cell=cell,donor=donor,factor=factor,baseline_separation=separation,destination_distance=dd,donor_distance=ds,donor_minus_destination_distance=ds-dd,nearest_reference=('donor' if ds<dd else 'destination') if separation>.01 else 'unresolved'))
    report=dict(response_metrics=summaries,comparisons=rows,scope=p['scope'],interpretation=p['interpretation'],limitations=p['limitations'])
    write_json(root/'comparison.json',report);return report


def run(root):
    root=Path(root);p=json.loads((root/'protocol.json').read_text());verify(p);done=0
    try:
        for key in p['backgrounds']:
            while not source_ready(p['sources'][key]):
                write_json(root/'status.json',dict(state='waiting',completed=done,total=15,source=key));time.sleep(30)
            verify(p);materialize(root,p,key)
            cp=json.loads((root/key/'protocol.json').read_text());verify(cp)
            for job in cp['jobs']:
                write_json(root/'status.json',dict(state='running',completed=done,total=15,background=key,current=moving.name(job),workers=1))
                moving.worker((str(root/key),job));done+=1
        assess(root,p);write_json(root/'status.json',dict(state='completed',completed=done,total=15))
    except Exception as exc:
        write_json(root/'status.json',dict(state='failed',completed=done,total=15,error=str(exc)));raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('command',choices=('prepare','run','assess'));parser.add_argument('--output',type=Path,default=Path('outputs/cell-exchange-response-moving'))
    a=parser.parse_args()
    if a.command=='prepare':prepare(a.output)
    elif a.command=='run':run(a.output)
    else:assess(a.output,json.loads((a.output/'protocol.json').read_text()))
