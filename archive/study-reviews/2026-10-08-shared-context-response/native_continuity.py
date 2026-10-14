"""Compare saved native frozen-parent trajectories; never integrate or restore."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib, json, sys
sys.dont_write_bytecode = True
import numpy as np
ROOT = Path('/home/sji/Documents/Programming/Python/embryogenesis-simulation')
OLD = ROOT/'output/studies/phase-carry-frozen-response'
NEW = ROOT/'output/studies/phase-carry-shared-context-response'
REPORT = Path('/tmp/shared-context-native-continuity.json')
pins = {}; stats = {}
read = lambda p:json.loads(Path(p).read_text())

def metadata(p):
    s=Path(p).stat()
    return (s.st_dev,s.st_ino,s.st_size,s.st_mtime_ns,s.st_mode)

def sha(p,expected=None):
    p=Path(p); before=metadata(p)
    with p.open('rb') as f: value=hashlib.file_digest(f,'sha256').hexdigest()
    assert metadata(p)==before and (expected is None or value==expected), str(p)
    assert str(p) not in pins or pins[str(p)]==value
    pins[str(p)]=value; stats[str(p)]=before
    return value

def signature(*values):
    h=hashlib.sha256()
    for a in values:
        a=np.ascontiguousarray(a)
        h.update(a.dtype.str.encode()); h.update(str(a.shape).encode()); h.update(a.tobytes())
    return h.hexdigest()

def problems(p):
    contexts={}
    for c in p['contexts']:
        sha(c['source'],p['input_sha256'][str(Path(c['source']).resolve())])
        with np.load(c['source'],allow_pickle=False) as z:
            contexts[c['key']]={k:z[k].copy() for k in ('delta','volumes','ids','initial')}
    values={}
    for t in p['trials']:
        sha(t['source'],p['input_sha256'][str(Path(t['source']).resolve())])
        with np.load(t['source'],allow_pickle=False) as z: initial=z['initial'].copy()
        x=contexts[t['context']]
        values[t['key']]=dict(initial=initial,signature=signature(x['delta'],x['volumes'],x['ids'],initial))
    return values

report=dict(passed=False,reviewed_utc=datetime.now(timezone.utc).isoformat(),
    scope='Descriptive saved-trajectory continuity audit only; no integration or scientific outcome test. A common-amount case is compared only after exact full-state/operator/volume/IDs equality.',
    integrations_performed=0,phase_arrays_loaded_or_copied=False,all_mismatches=[])
try:
    po=read(OLD/'protocol.json'); pn=read(NEW/'protocol.json')
    ho=sha(OLD/'protocol.json'); hn=sha(NEW/'protocol.json')
    compact=read(ROOT/'docs/phase_carry_frozen_response_assessment.json')
    assert compact['full_study_accepted'] and compact['protocol_sha256']==ho
    sha(OLD/'status.json',compact['original_status_sha256'])
    sha(OLD/'summary.json',compact['original_summary_sha256'])
    for k in ('beta','da','db','start','duration','interval','levels','device','runtime'):
        assert po[k]==pn[k],k
    vo,vn=problems(po),problems(pn)
    native=[]
    for t in pn['trials']:
        matches=[o['key'] for o in po['trials'] if vo[o['key']]['signature']==vn[t['key']]['signature']]
        if matches:native.append(dict(new_trial=t['key'],target=t['target'],dose=t['dose'],parent_matches=matches))
    assert len(native)==27
    assert sum(t['target'] is None for t in native)==9
    assert sum(t['dose']=='fractional' for t in native)==12
    assert sum(t['dose']=='amount' for t in native)==6
    levels=[]; waiting=[]
    for level in pn['levels']:
        l=level['level']; a=OLD/'production'/l; b=NEW/'production'/l
        if not all((b/f).exists() for f in ('paths.npz','result.json')):
            waiting.append(l);continue
        ra,rb=read(a/'result.json'),read(b/'result.json')
        sha(a/'result.json');sha(b/'result.json')
        assert ra['protocol_sha256']==ho and rb['protocol_sha256']==hn
        assert ra['level']==rb['level']==level
        sha(a/'paths.npz',ra['paths_sha256']);sha(b/'paths.npz',rb['paths_sha256'])
        with np.load(a/'paths.npz',allow_pickle=False) as z, np.load(b/'paths.npz',allow_pickle=False) as q:
            at=z['times'].copy();bt=q['times'].copy();av=z['trajectory'].copy();bv=q['trajectory'].copy()
            ak=list(z['trial_keys']);bk=list(q['trial_keys'])
        assert np.array_equal(at,bt) and len(at)==401 and at[0]==0. and at[-1]==60.
        assert av.shape==(401,33,2,16) and bv.shape==(401,69,2,16)
        assert all(np.isfinite(v).all() and np.min(v)>0 for v in (av,bv))
        assert set(bk)==set(vn) and len(bk)==len(set(bk))
        for t in pn['trials']:assert np.array_equal(bv[0,bk.index(t['key'])],vn[t['key']]['initial'])
        cases=[]
        for match in native:
            j=bk.index(match['new_trial'])
            for oldkey in match['parent_matches']:
                i=ak.index(oldkey)
                assert np.array_equal(av[0,i],bv[0,j])
                cases.append(dict(**match,parent_trial=oldkey,
                    byte_identical_trajectory=np.array_equal(av[:,i],bv[:,j]),
                    max_absolute_log_difference=float(np.abs(np.log(av[:,i]/bv[:,j])).max()),
                    max_absolute_concentration_difference=float(np.abs(av[:,i]-bv[:,j]).max())))
        levels.append(dict(level=level,logical_native_cases=len(native),exact_comparisons=len(cases),
            cases=cases,max_absolute_log_difference=max(c['max_absolute_log_difference'] for c in cases),
            max_absolute_concentration_difference=max(c['max_absolute_concentration_difference'] for c in cases)))
    before=dict(pins)
    for name,value in list(before.items()):sha(name,value)
    for name,value in stats.items():assert metadata(name)==value
    report.update(passed=not waiting,prepared_inputs_verified=True,old_protocol_sha256=ho,
        new_protocol_sha256=hn,native_exact_initial_problem_matches=native,
        native_count=27,native_controls=9,native_fractional=12,native_common_amount=6,
        completed_levels=levels,waiting_for_saved_levels=waiting,
        immutable_before_sha256=before,immutable_after_sha256=dict(pins),
        immutable_before_after_equal=before==pins,
        script_path=str(Path(__file__)),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
except BaseException as exc:
    import traceback
    report.update(error=str(exc),traceback=traceback.format_exc(),all_mismatches=[str(exc)])
    raise
finally:
    REPORT.write_text(json.dumps(report,indent=2)+'\n')
    print('Native saved-path continuity',report.get('prepared_inputs_verified'),report.get('waiting_for_saved_levels'),REPORT)
