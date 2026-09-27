"""Compare identical mature-state workloads against the frozen reference kernel.

Run from the repository: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python
benchmarks/mechanics.py --checkpoint outputs/shape-persistence/full/final_state.npz
--output outputs/kernel-repeat
"""
import argparse
from copy import deepcopy
import hashlib
import json
import os
from pathlib import Path
import platform
import statistics
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'tests')]
import numpy as np
import scipy
from embryo.model import Simulation
from reference_mechanics import ReferenceSimulation


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checkpoint',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--steps',type=int,default=30)
    parser.add_argument('--repeats',type=int,default=3)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError('choose a fresh benchmark directory')
    if args.steps<1 or args.repeats<1:raise ValueError('positive steps and repeats required')
    source=Simulation.restore(args.checkpoint)
    fast=deepcopy(source);reference=deepcopy(source);reference.__class__=ReferenceSimulation
    for _ in range(args.steps):fast.step();reference.step()
    for name in ('phi','fate','activator','inhibitor','polarity','target'):
        np.testing.assert_array_equal(getattr(fast,name),getattr(reference,name))
    timings={'reference':[],'optimized':[]}
    for repeat in range(args.repeats):
        order=('reference','optimized') if repeat%2==0 else ('optimized','reference')
        for name in order:
            sim=deepcopy(source)
            if name=='reference':sim.__class__=ReferenceSimulation
            sim.step();start=time.perf_counter()
            for _ in range(args.steps):sim.step()
            timings[name].append(time.perf_counter()-start)
    result={'checkpoint':str(args.checkpoint),'checkpoint_sha256':hashlib.sha256(args.checkpoint.read_bytes()).hexdigest(),
            'steps_per_repeat':args.steps,'repeats':args.repeats,'timings_seconds':timings,
            'median_speedup':statistics.median(timings['reference'])/statistics.median(timings['optimized']),
            'exact_state_match':True,'reference':'tests/reference_mechanics.py',
            'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
                           'OPENBLAS_NUM_THREADS':os.environ.get('OPENBLAS_NUM_THREADS'),
                           'OMP_NUM_THREADS':os.environ.get('OMP_NUM_THREADS')},
            'scope':'Single-process checkpoint workload; not a GPU or scaling claim.'}
    args.output.mkdir(parents=True)
    (args.output/'analysis.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
