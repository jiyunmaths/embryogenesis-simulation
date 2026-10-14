"""Group identical physical operators without coupling chemical trajectories."""
import argparse
import fcntl
import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import get_context
import os
from pathlib import Path
import time
import numpy as np
import torch
from scipy import sparse
from . import phase_carry_frozen_response as frozen
from . import phase_carry_recipient_state_response as base
from .feedback_long import digest
from .neighbor_context import read
from .parameter_robustness import verify
write_json=frozen.write_json
save_npz=frozen.save_npz


def batch_layout(p):
    groups={}
    for trial in p['trials']:
        _,delta,volumes,initial,ids=frozen.load_context(p,trial)
        key=tuple((a.dtype.str,a.shape,a.tobytes()) for a in (delta,volumes,ids))
        if key not in groups:groups[key]=dict(delta=delta,values=[])
        groups[key]['values'].append((trial,initial))
    width=max(len(g['values']) for g in groups.values())
    matrices=[];states=[];mapping=[];count=None
    for gi,g in enumerate(groups.values()):
        if count is None:count=g['delta'].shape[0]
        if g['delta'].shape!=(count,count):raise ValueError('Equal compartment counts required')
        values=[u for t,u in g['values']]
        if g['values'][0][0]['target'] is not None:raise ValueError('Padding must duplicate a control')
        states.append(np.stack(values+[values[0].copy() for _ in range(width-len(values))],axis=-1))
        matrices.append(g['delta'])
        mapping.extend(dict(trial=t['key'],context_index=gi,column=i) for i,(t,u) in enumerate(g['values']))
    # Restore logical trial order for prefix checks and immutable observation keys.
    by={m['trial']:m for m in mapping};mapping=[by[t['key']] for t in p['trials']]
    return np.concatenate(states,axis=1),sparse.block_diag(matrices).toarray(),mapping,count


def prepare(root,parent=Path('output/studies/phase-carry-shared-context-response'),
            review_file=Path('docs/phase_carry_shared_context_response_assessment.json'),
            protocol_doc=Path('docs/phase_carry_recipient_state_response_grouped.md')):
    p=base.prepare(root,parent,review_file,protocol_doc)
    initial,matrix,mapping,count=batch_layout(p)
    p['batching']=dict(mode='identical physical operators with independent chemical columns',
        physical_blocks=matrix.shape[0]//count,width=initial.shape[-1],logical_paths=len(mapping),
        padding_paths=initial.shape[1]*initial.shape[2]//count-len(mapping))
    if (p['batching']['physical_blocks'],p['batching']['width'],p['batching']['padding_paths'])!=(9,15,42):
        raise ValueError('Unexpected measured operator grouping')
    p['source_sha256'].update(frozen.pin([Path(__file__),Path('tests/test_phase_carry_recipient_state_response_grouped.py')]))
    p['input_sha256'].update(frozen.pin([Path('docs/phase_carry_recipient_state_response.md')]))
    write_json(Path(root)/'protocol.json',p)
    return p


def production(root, p, level):
    root = Path(root); ph = digest(root/'protocol.json'); folder = root/'production'/level['level']
    file, report, checkpoint = folder/'paths.npz', folder/'result.json', folder/'checkpoint.npz'
    if report.exists():
        r = read(report)
        if (r['protocol_sha256'] != ph or r['level'] != level or r['paths_sha256'] != digest(file)
                or r['checkpoint_sha256'] != digest(checkpoint)):
            raise ValueError('Changed completed production evidence')
        return r
    initial, delta, mapping, count = batch_layout(p); device = frozen.runtime_check(p)
    dt = level['dt']; total = frozen._steps(p['duration'], dt); every = frozen._steps(p['interval'], dt)
    save_every = frozen._steps(p['checkpoint_interval'], dt)
    for c in p['contexts']:
        with np.load(c['source']) as z:
            if frozen.substeps(z['delta'], dt, p['beta'], p['da'], p['db']) != frozen.substeps(delta, dt, p['beta'], p['da'], p['db']):
                raise ValueError('Batched step changes positivity substeps')
    completed, wall, history = 0, 0., [frozen.unpack(initial, mapping, count)]
    state = initial
    if checkpoint.exists():
        commit = read(folder/'checkpoint-commit.json')
        if commit['protocol_sha256'] != ph or commit['checkpoint_sha256'] != digest(checkpoint):
            raise ValueError('Changed or incompletely committed production checkpoint; retain for review')
        with np.load(checkpoint) as z:
            meta = json.loads(str(z['metadata'])); completed = meta['step']; wall = meta['wall_seconds']
            state, history = z['state'].copy(), list(z['history'].copy())
            if (meta['protocol_sha256'] != ph or meta['level'] != level
                    or completed % save_every != 0 or not 0 <= completed <= total
                    or state.shape != initial.shape or np.asarray(history).shape != (completed//every+1, len(mapping), 2, count)
                    or len(history) != completed//every+1 or not np.array_equal(history[0], frozen.unpack(initial, mapping, count))
                    or not np.array_equal(history[-1], frozen.unpack(state, mapping, count))):
                raise ValueError('Checkpoint clock/preparation/history mismatch')
            frozen.valid_state(state); frozen.valid_state(np.array(history))
    a, b, matrix = [torch.as_tensor(x, dtype=torch.float64, device=device) for x in (state[0], state[1], delta)]
    batch_substeps = frozen.substeps(delta, dt, p['beta'], p['da'], p['db'])
    started = time.perf_counter()
    for step in range(completed, total+1):
        if step % save_every == 0:
            current = np.stack((a.detach().cpu().numpy(), b.detach().cpu().numpy()))
            meta = dict(protocol_sha256=ph, level=level, step=step, wall_seconds=wall+time.perf_counter()-started)
            save_npz(checkpoint, state=current, history=np.array(history), metadata=np.array(json.dumps(meta)))
            write_json(folder/'checkpoint-commit.json', dict(protocol_sha256=ph,
                checkpoint_sha256=digest(checkpoint), level=level, step=step))
            write_json(folder/'status.json', dict(state='running', elapsed=step*dt, horizon=p['duration']))
            print('Grouped recipient production:', level['level'], step*dt, flush=True)
        if step == total:
            break
        a, b = frozen.checked_torch_step(a, b, matrix, dt, p, batch_substeps)
        if (step+1) % every == 0:
            history.append(frozen.unpack(np.stack((a.detach().cpu().numpy(), b.detach().cpu().numpy())), mapping, count))
    paths = np.array(history); times = np.arange(len(paths))*p['interval']
    if paths.shape != (len(times), len(mapping), 2, count) or len(times) != frozen._steps(p['duration'], p['interval'])+1:
        raise ValueError('Incomplete batched production')
    save_npz(file, times=times, trajectory=paths, trial_keys=np.array([m['trial'] for m in mapping]))
    r = dict(protocol_sha256=ph, level=level, paths_sha256=digest(file), checkpoint_sha256=digest(checkpoint),
        actual_paths=len(mapping), padding_paths=initial.shape[1]*initial.shape[-1]//count-len(mapping),
        wall_seconds=wall+time.perf_counter()-started)
    write_json(report, r); write_json(folder/'status.json', dict(state='completed', elapsed=p['duration']))
    return r


def run(root, workers=4):
    root = Path(root).resolve()
    if not 1 <= workers <= 4:
        raise ValueError('Use one to four single-thread CPU workers')
    with (root/'coordinator.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        p=read(root/'protocol.json');verify(p);frozen.runtime_check(p);torch.set_num_threads(1)
        for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
            os.environ[name]='1'
        started=time.perf_counter()
        try:
            with ProcessPoolExecutor(max_workers=workers,mp_context=get_context('spawn')) as pool:
                write_json(root/'status.json',dict(state='running',stage='independent_references',completed=0,total=len(p['trials'])))
                pending=[pool.submit(frozen.reference_worker,(root,t)) for t in p['trials']]
                for i,future in enumerate(as_completed(pending),1):
                    if not future.result()['passed']:
                        raise RuntimeError('Raw dual-reference gate failed')
                    write_json(root/'status.json',dict(state='running',stage='independent_references',completed=i,total=len(p['trials'])))
                if not base.reference_qualification(root,p)['passed']:
                    raise RuntimeError('Normalized reference gate failed; production blocked')
                # CPU verification and GPU production have disjoint output paths.
                pending=[pool.submit(frozen.cpu_worker,(root,t,l)) for l in p['levels'] for t in p['trials']]
                for level in p['levels']:
                    done=[f for f in pending if f.done()]
                    for future in done:future.result()
                    write_json(root/'status.json',dict(state='running',stage='production_and_cpu_verification',
                        level=level['level'],cpu_completed=len(done),cpu_total=len(pending),total_production_paths=p['total_production_paths']))
                    production(root,p,level)
                for i,future in enumerate(as_completed(pending),1):
                    future.result()
                    write_json(root/'status.json',dict(state='running',stage='cpu_verification',completed=i,total=len(pending)))
            result=base.assess(root)
            write_json(root/'timing.json',dict(coordinator_wall_seconds=time.perf_counter()-started,
                cpu_verification_overlaps_gpu=True,workers=workers))
            write_json(root/'status.json',dict(state='completed' if result['passed'] else 'completed_with_unresolved_checks',
                passed=result['passed'],primary_paths=p['primary_paths'],diagnostic_paths=p['diagnostic_paths'],
                total_production_paths=p['total_production_paths'],independent_histories=3,new_histories=0))
            return result
        except Exception as error:
            write_json(root/'status.json',dict(state='failed',error=str(error)));raise


if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('action',choices=('prepare','run','assess'))
    a.add_argument('--output',type=Path,default=Path('output/studies/phase-carry-recipient-state-response-grouped'))
    a.add_argument('--parent',type=Path,default=Path('output/studies/phase-carry-shared-context-response'))
    a.add_argument('--review',type=Path,default=Path('docs/phase_carry_shared_context_response_assessment.json'))
    a.add_argument('--protocol-doc',type=Path,default=Path('docs/phase_carry_recipient_state_response_grouped.md'))
    a.add_argument('--workers',type=int,default=4);args=a.parse_args()
    result=prepare(args.output,args.parent,args.review,args.protocol_doc) if args.action=='prepare' else run(args.output,args.workers) if args.action=='run' else base.assess(args.output)
    print({k:v for k,v in result.items() if k in ('passed','primary_paths','diagnostic_paths','total_production_paths')})
