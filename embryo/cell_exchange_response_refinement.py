"""Targeted timestep halving of the completed moving exchange-response assay."""
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import fcntl
import json
import multiprocessing as mp
from pathlib import Path
import shutil
import numpy as np
from . import cell_exchange_response_moving as study
from .cell_response_moving import name
from .cell_response import response_metrics
from .cell_response_exchange import waveform_distance
from .cell_response_moving_refinement import CRITERIA,response_error
from .feedback_survival_validation import retime
from .native_response_resume import worker
from .feedback_long import digest
from .resolution import write_json

BACKGROUNDS=('unexchanged','fresh_exchange')


def verify(root):
    p=json.loads((Path(root)/'protocol.json').read_text());study.verify(p);return p


def prepare(root,baseline,workers=3,threads=2):
    root,baseline=Path(root).resolve(),Path(baseline).resolve()
    if root.exists():raise FileExistsError(root)
    if workers<1 or threads<1 or workers*threads>6:raise ValueError('Use at most six native threads')
    old=verify(baseline)
    if json.loads((baseline/'status.json').read_text())['state']!='completed':raise ValueError('Requires completed baseline')
    transition=baseline/'native-backend-transition.json';manifest=json.loads(transition.read_text())
    for f,h in {**manifest['source_sha256'],**manifest['validation_sha256'],**manifest['completed_results']}.items():
        if digest(f)!=h:raise ValueError('Backend evidence changed: '+f)
    jobs={key:[dict(family=key,target=None,factor=1.)]+[dict(family=key,target=cell,factor=.9) for cell in old['targets']] for key in BACKGROUNDS}
    files=[baseline/'protocol.json',baseline/'status.json',baseline/'comparison.json',transition]
    # Verify all selected original results before preparing any fine input.
    for key in BACKGROUNDS:
        child=baseline/key;cp=verify(child);ph=digest(child/'protocol.json')
        files.extend(child/f for f in ('protocol.json','source.npz','initial_states.npz'))
        for job in jobs[key]:
            folder=child/name(job);r=json.loads((folder/'result.json').read_text())
            if not r['quality_pass'] or r['protocol_sha256']!=ph:raise ValueError('Invalid baseline result')
            files.extend(folder/f for f in ('result.json','history.json'))
    root.mkdir(parents=True)
    for key in BACKGROUNDS:
        child=root/key;child.mkdir();source=baseline/key
        sim=retime(source/'source.npz',child/'refined-source.npz',old['dt']/2)
        if sim.divisions or abs(sim.time-old['start'])>1e-9:raise ValueError('Invalid retimed physical state')
        shutil.copy2(source/'initial_states.npz',child/'initial_states.npz')
        files.extend(child/f for f in ('refined-source.npz','initial_states.npz'))
        cp=dict(checkpoint=str(child/'refined-source.npz'),start=old['start'],dt=old['dt']/2,
                duration=old['duration'],interval=old['interval'],jobs=jobs[key],source_sha256={},input_sha256={})
        write_json(child/'protocol.json',cp);files.append(child/'protocol.json')
    code=[Path(__file__),*[Path(__file__).with_name(f) for f in ('native_response_resume.py','fast_response_resume.py','native_mechanics.py','native_mechanics.cpp','fast_mechanics.py','fast_mechanics.c','cell_exchange_response_moving.py','cell_response_moving_refinement.py','cell_response.py','cell_response_exchange.py','feedback_survival_validation.py','attribute_development.py','model.py','polarity.py','signaling.py','transport.py','feedback_long.py','resolution.py')]]
    p=dict(baseline=str(baseline),backgrounds=BACKGROUNDS,targets=old['targets'],factor=.9,start=old['start'],dt=old['dt']/2,duration=old['duration'],interval=old['interval'],workers=workers,threads_per_worker=threads,
        refinement_criteria=CRITERIA,classification_criteria=dict(reference_separation_min=.01,nearest_reference_unchanged=True),
        scope='Same-state mature t=210 continuation: one timestep halving for untouched and fresh-exchange backgrounds, matched controls and -10% activator pulses in both original exchange targets. Six fine runs, full 60-unit horizon. Preserves physical fields, chemistry, polarity, IDs and random streams; retimes integer clocks only.',
        limits='Single developmental history; excludes +10% pulses, pre-relaxed exchange and spatial refinement. Starting exchange-generated t=210 states are held fixed, so this does not refine their prior t=150-to-210 formation trajectory or establish autonomous identity.',
        source_sha256={str(f.resolve()):digest(f) for f in code},input_sha256={str(f):digest(f) for f in files})
    write_json(root/'protocol.json',p);write_json(root/'status.json',dict(state='prepared',completed=0,total=6))
    return p


def load(root,key,jobs,times):
    child=Path(root)/key;p=verify(child);ph=digest(child/'protocol.json');paths={};histories={}
    with np.load(child/'initial_states.npz') as d:ids=d['ids'];masses=d['masses']
    for job in jobs:
        folder=child/name(job);r=json.loads((folder/'result.json').read_text())
        if not r['quality_pass'] or r['protocol_sha256']!=ph:raise ValueError('Invalid result: '+str(folder))
        h=json.loads((folder/'history.json').read_text())
        if len(h)!=len(times) or not np.allclose([x['elapsed'] for x in h],times,rtol=0,atol=1e-9) or any(not np.array_equal(x['ids'],ids) for x in h):raise ValueError('History alignment mismatch')
        x=np.array([row['chemistry'] for row in h])
        if not np.isfinite(x).all() or np.any(x<=0):raise ValueError('Invalid chemical path')
        paths[job['target']]=x;histories[job['target']]=h
    return paths,histories,ids,masses


def recover_error(a,b):
    return 0. if a is None and b is None else (None if a is None or b is None else abs(a-b))


def assess(root):
    root=Path(root);p=verify(root);baseline=Path(p['baseline']);criteria=p['refinement_criteria']
    times=np.arange(round(p['duration']/p['interval'])+1)*p['interval'];datasets={};waves={};rows=[];raw=[]
    for parent in (baseline,root):
        waves[parent]={};datasets[parent]={}
        for key in p['backgrounds']:
            jobs=json.loads((root/key/'protocol.json').read_text())['jobs']
            data,h,ids,masses=load(parent,key,jobs,times);datasets[parent][key]=(data,h,ids,masses);waves[parent][key]={}
            for cell in p['targets']:
                i=list(ids).index(cell);waves[parent][key][cell]=study.local_response(data[cell],data[None],i,p['factor'])
    for key in p['backgrounds']:
        coarse,ch,ids,masses=datasets[baseline][key];fine,fh,fine_ids,fine_masses=datasets[root][key]
        if not np.array_equal(ids,fine_ids) or not np.array_equal(masses,fine_masses):raise ValueError('Refinement altered IDs or starting masses')
        for target in (None,*p['targets']):raw.append(float(abs(np.log(fine[target]/coarse[target])).max()))
        for cell in p['targets']:
            i=list(ids).index(cell);a=response_metrics(coarse[cell],coarse[None],times,masses,i,p['factor']);b=response_metrics(fine[cell],fine[None],times,masses,i,p['factor'])
            error=response_error(coarse[cell],coarse[None],fine[cell],fine[None],p['factor'])
            metric='target_activator_log_auc_per_log_pulse';auc=abs(b[metric]/a[metric]-1)
            recovery={k:recover_error(a[k],b[k]) for k in ('target_recovery_time','network_recovery_time')}
            passed=error<=criteria['max_normalized_response_error'] and auc<=criteria['max_relative_auc_error'] and all(v is not None and v<=criteria['max_recovery_time_error']+1e-9 for v in recovery.values())
            rows.append(dict(background=key,cell=cell,normalized_response_error=error,relative_auc_error=auc,recovery_time_errors=recovery,coarse_metrics=a,fine_metrics=b,passed=bool(passed)))
    classification=[]
    for cell in p['targets']:
        donor=study.donor_id(cell,p['targets']);values=[]
        for parent in (baseline,root):
            dest=waves[parent]['unexchanged'][cell];donor_wave=waves[parent]['unexchanged'][donor];wave=waves[parent]['fresh_exchange'][cell]
            separation=waveform_distance(dest,donor_wave,times);dd=waveform_distance(wave,dest,times);ds=waveform_distance(wave,donor_wave,times)
            label=('donor' if ds<dd else 'destination') if separation>p['classification_criteria']['reference_separation_min'] else 'unresolved'
            values.append(dict(separation=separation,destination_distance=dd,donor_distance=ds,nearest_reference=label))
        passed=all(v['nearest_reference']!='unresolved' for v in values) and values[0]['nearest_reference']==values[1]['nearest_reference']
        classification.append(dict(cell=cell,donor=donor,coarse=values[0],fine=values[1],passed=passed))
    passed=all(row['passed'] for row in rows+classification) and max(raw)<=criteria['max_raw_chemical_log_error']
    report=dict(passed=bool(passed),max_raw_chemical_log_error=max(raw),criteria=criteria,trials=rows,classifications=classification,scope=p['scope'],limits=p['limits'])
    write_json(root/'refinement.json',report);return report


def run(root):
    root=Path(root).resolve()
    with (root/'run.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);p=verify(root);rows=[]
        for key in p['backgrounds']:
            cp=verify(root/key)
            for job in cp['jobs']:rows.append((str(root/key),job,p['threads_per_worker']))
        done=sum((Path(child)/name(job)/'result.json').exists() for child,job,_ in rows)
        try:
            write_json(root/'status.json',dict(state='running',completed=done,total=len(rows),workers=p['workers'],threads_per_worker=p['threads_per_worker']))
            pending=[row for row in rows if not (Path(row[0])/name(row[1])/'result.json').exists()]
            with ProcessPoolExecutor(max_workers=p['workers'],mp_context=mp.get_context('spawn')) as pool:
                futures={pool.submit(worker,row):row for row in pending}
                for future in as_completed(futures):
                    future.result();done+=1
                    write_json(root/'status.json',dict(state='running',completed=done,total=len(rows),workers=p['workers'],threads_per_worker=p['threads_per_worker'],last_completed=name(futures[future][1])))
            result=assess(root);write_json(root/'status.json',dict(state='completed',completed=done,total=len(rows),passed=result['passed']))
        except Exception as exc:
            write_json(root/'status.json',dict(state='failed',completed=done,total=len(rows),error=str(exc)));raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('command',choices=('prepare','run','assess'));parser.add_argument('--output',type=Path,default=Path('outputs/cell-exchange-response-refined'));parser.add_argument('--baseline',type=Path,default=Path('outputs/cell-exchange-response-moving'));a=parser.parse_args()
    if a.command=='prepare':prepare(a.output,a.baseline)
    elif a.command=='run':run(a.output)
    else:assess(a.output)
