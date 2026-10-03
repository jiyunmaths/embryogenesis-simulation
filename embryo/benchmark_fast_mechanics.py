"""Paired opt-in mechanics validation and throughput benchmark on a checkpoint."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import platform
import time
import numpy as np
from .attribute_development import AttributeSimulation
from .fast_mechanics import FastAttributeSimulation,kernel
from .feedback_long import digest
from .resolution import write_json


def trial(args):
    checkpoint,steps=args;sim=FastAttributeSimulation.restore(checkpoint);kernel();sim.step()
    start=time.perf_counter()
    for _ in range(steps):sim.step()
    return dict(seconds=time.perf_counter()-start,steps=steps)


def run(checkpoint,output,steps=48,concurrency=False):
    output=Path(output)
    if output.exists():raise FileExistsError(output)
    output.mkdir(parents=True);checkpoint=Path(checkpoint).resolve();kernel()
    a=AttributeSimulation.restore(checkpoint);b=FastAttributeSimulation.restore(checkpoint)
    a.step();b.step() # warm both paths before timing; both states advanced equally
    elapsed={'original':[],'optimized':[]};maxima=dict(phi_abs=0.,chemical_log_abs=0.,polarity_abs=0.,relative_volume=0.,clipping=0.)
    for k in range(steps):
        for label,sim in ([('original',a),('optimized',b)] if k%2==0 else [('optimized',b),('original',a)]):
            start=time.perf_counter();sim.step();elapsed[label].append(time.perf_counter()-start)
        maxima['phi_abs']=max(maxima['phi_abs'],float(abs(a.phi-b.phi).max()))
        maxima['chemical_log_abs']=max(maxima['chemical_log_abs'],float(abs(np.log(np.array([a.activator,a.inhibitor])/np.array([b.activator,b.inhibitor]))).max()))
        maxima['polarity_abs']=max(maxima['polarity_abs'],float(abs(a.polarity-b.polarity).max()))
        maxima['relative_volume']=max(maxima['relative_volume'],float(abs(a.volumes()/b.volumes()-1).max()))
        maxima['clipping']=max(maxima['clipping'],a.clipped_fraction,b.clipped_fraction)
    tolerances=dict(phi_abs=2e-6,chemical_log_abs=1e-5,polarity_abs=1e-6,relative_volume=1e-6,clipping=0.)
    passed=all(maxima[k]<=v for k,v in tolerances.items()) and a.ids.tolist()==b.ids.tolist() and a.rng.bit_generator.state==b.rng.bit_generator.state and a.signal_rng.bit_generator.state==b.signal_rng.bit_generator.state
    report=dict(checkpoint=str(checkpoint),checkpoint_sha256=digest(checkpoint),steps=steps,model_time_advanced=(steps+1)*a.config.dt,dt=a.config.dt,grid=a.config.grid,cells=len(a.ids),
        timing_seconds_per_step={key:float(np.mean(value)) for key,value in elapsed.items()},speedup=float(sum(elapsed['original'])/sum(elapsed['optimized'])),
        maximum_discrepancies=maxima,tolerances=tolerances,passed=passed,
        provenance={str(f.resolve()):digest(f) for f in [Path(__file__),Path(__file__).with_name('fast_mechanics.py'),Path(__file__).with_name('fast_mechanics.c'),Path(__file__).with_name('attribute_development.py'),Path(__file__).with_name('model.py')]},
        platform=platform.platform(),scope='Short local trajectory agreement, not full developmental or long-horizon scientific equivalence. Interleaved timings under concurrent production load; observations and checkpoint IO excluded. Existing simulations were not interrupted.')
    write_json(output/'comparison.json',report)
    print(json.dumps(report,indent=2),flush=True)
    if not passed:raise RuntimeError('Optimized dynamics failed local agreement')
    if concurrency:
        rows=[]
        for workers in (1,2,3,4):
            with ProcessPoolExecutor(max_workers=workers) as pool:
                start=time.perf_counter();results=list(pool.map(trial,[(str(checkpoint),8)]*workers));seconds=time.perf_counter()-start
            row=dict(benchmark_workers=workers,total_step_throughput=8*workers/seconds,wall_seconds=seconds,individual_runs=results)
            rows.append(row);write_json(output/'concurrency.json',dict(trials=rows,scope='Benchmark worker counts only; two production simulations remain active in addition. Includes startup, checkpoint reads and one warmup per worker; short noisy throughput screen, not isolated scaling evidence.'))
            print(json.dumps(row),flush=True)
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--checkpoint',type=Path,default=Path('outputs/cell-exchange-response-moving/unexchanged/unexchanged_control/latest_state.npz'));parser.add_argument('--output',type=Path,required=True);parser.add_argument('--steps',type=int,default=48);parser.add_argument('--concurrency',action='store_true');a=parser.parse_args()
    if a.steps<1:raise ValueError('Positive benchmark step count required')
    run(a.checkpoint,a.output,a.steps,a.concurrency)
