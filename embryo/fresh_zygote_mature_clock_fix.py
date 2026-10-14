"""Zero-elapsed checkpoint scheduling fix; retain the failed frozen attempt.

Only the mature coordinator's integer scheduler accepts elapsed zero. No solver,
kernel, chemical clock, checkpoint array or acceptance threshold changes.
"""
import argparse
from copy import deepcopy
from pathlib import Path
import os

from . import fresh_zygote_mature_qualification as f
from .qualify_development_gpu import read, write, digest, pin

OLD=Path('output/studies/fresh-zygote-mature-qualification')
OLD_SHA='f6c099718c3d368646912660fb91fc2ce1280ce11130fb105cf07790f9811131'
AUDIT=Path('archive/study-launches/2026-10-12-fresh-zygote-mature-qualification/prelaunch-audit.json')
LAUNCH=Path('archive/study-launches/2026-10-12-fresh-zygote-mature-clock-fix')
DESIGN=Path('docs/fresh_zygote_mature_clock_fix.md')
ORIGINAL_STEPS=f._steps
ORIGINAL_PREFIX=f.pilot.prefix
FIRST='h1101_off_coarse'


def elapsed_steps(duration,dt):
    if duration==0:
        ORIGINAL_STEPS(dt,dt)  # Retain positive, finite dt validation.
        return 0
    return ORIGINAL_STEPS(duration,dt)


def verify(p):
    f.q.verify(p)
    old=read(OLD/'protocol.json');audit=read(AUDIT)
    if digest(OLD/'protocol.json')!=OLD_SHA or not audit['passed'] or audit['protocol_sha256']!=OLD_SHA:
        raise ValueError('Immutable, independently audited starts required')
    changed={'jobs','study_root','source_sha256','input_sha256','stop_policy'}
    if any(p[key]!=value for key,value in old.items() if key not in changed):
        raise ValueError('Only scheduling/paths/provenance may change')
    expected=[dict(job,source=str(Path(p['study_root'])/job['key']/'source.npz')) for job in old['jobs']]
    if p['jobs']!=expected:raise ValueError('Clock-fix endpoint design changed')
    for mapping in ('source_sha256','input_sha256'):
        if any(p[mapping].get(k)!=v for k,v in old[mapping].items()):raise ValueError('Inherited pin changed')
    for job,original in zip(p['jobs'],old['jobs']):
        source=Path(job['source'])
        if not source.is_symlink() or source.resolve()!=Path(original['source']).resolve():
            raise ValueError('Continuation starts must alias the unchanged audited bytes')
    status=read(OLD/'status.json');result=read(OLD/FIRST/'gpu-result.json')
    if status['state']!='failed' or result['final_time']!=18. or result['failure']['message']!='duration and dt must be positive and finite':
        raise ValueError('Original scheduling failure must remain explicit')
    prefix=read(OLD/FIRST/'actual-start-prefix.json')
    if (not prefix['passed'] or not prefix['exact_gpu_restart'] or prefix['steps']!=8 or
            prefix['protocol_sha256']!=OLD_SHA or any(v>p['prefix_limits'][k] for k,v in prefix['errors'].items())):
        raise ValueError('Reused actual-start prefix must remain accepted')
    if p['reused_prefix_case']!=FIRST or p['scheduling_fix']!='Elapsed zero maps to integer step zero; positive durations retain the original resolver.':
        raise ValueError('Unregistered scheduler/prefix change')


def prepare(root):
    root=Path(root).resolve()
    if root.exists():raise FileExistsError('Fresh corrected-study directory required')
    old=read(OLD/'protocol.json');f.q.verify(old)
    if digest(OLD/'protocol.json')!=OLD_SHA:raise ValueError('Original frozen protocol changed')
    root.mkdir(parents=True);p=deepcopy(old);p['study_root']=str(root)
    p['jobs']=[dict(job,source=str(root/job['key']/'source.npz')) for job in old['jobs']]
    for job,original in zip(p['jobs'],old['jobs']):
        folder=root/job['key'];folder.mkdir();Path(job['source']).symlink_to(Path(original['source']))
    p['source_sha256'].update(pin([Path(__file__),DESIGN,Path('tests/test_fresh_zygote_mature_clock_fix.py'),LAUNCH/'audit_prelaunch.py',LAUNCH/'launch.py']))
    inputs=[OLD/'protocol.json',OLD/'status.json',OLD/FIRST/'gpu-result.json',AUDIT,
        OLD/FIRST/'actual-start-prefix.json',OLD/FIRST/'failure-state.npz']
    prefix=read(OLD/FIRST/'actual-start-prefix.json');inputs+=list(prefix['evidence_sha256'])
    p['input_sha256'].update(pin(inputs));p['reused_prefix_case']=FIRST
    p['scheduling_fix']='Elapsed zero maps to integer step zero; positive durations retain the original resolver.'
    p['stop_policy']='Preserve the original zero-time scheduling failure. Stop on any corrected-study prefix, physical, restart or pair failure; no threshold relaxation, fallback or long extension.'
    write(root/'protocol.json',p);verify(p)
    write(root/'status.json',dict(state='prepared',protocol_sha256=digest(root/'protocol.json')))
    return p


def prefix(root,job,p):
    if job['key']!=FIRST:return ORIGINAL_PREFIX(root,job,p)
    report=read(OLD/FIRST/'actual-start-prefix.json');output=Path(root)/FIRST/'actual-start-prefix.json'
    if output.exists():raise FileExistsError('Never overwrite a reused prefix record')
    reused=dict(report,reused=True,original_protocol_sha256=report['protocol_sha256'],
        original_report=str((OLD/FIRST/'actual-start-prefix.json').resolve()),
        original_report_sha256=digest(OLD/FIRST/'actual-start-prefix.json'),
        protocol_sha256=digest(Path(root)/'protocol.json'))
    write(output,reused);return reused


def run(root):
    # Process-local adapter for coordinator scheduling and audited provenance.
    # Neither the resolution module nor the original source files are changed.
    f._steps=elapsed_steps;f.verify=verify;f.pilot.prefix=prefix
    return f.run(root)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=('prepare','run'));parser.add_argument('root',type=Path)
    args=parser.parse_args();(prepare if args.action=='prepare' else run)(args.root)


if __name__=='__main__':main()
