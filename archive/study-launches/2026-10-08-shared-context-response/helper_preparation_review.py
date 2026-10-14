"""Literal independent equality checks of read-only shared-context helpers."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, inspect, json, sys
sys.dont_write_bytecode=True
import numpy as np
from embryo import phase_carry_shared_context_response as s
ROOT=Path('/home/sji/Documents/Programming/Python/embryogenesis-simulation')
REPORT=Path('/tmp/shared-context-helper-preparation-review.json')
digest=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
source=Path(s.__file__);before=digest(source)
old,proof,pins=s.accepted_parent(ROOT/'output/studies/phase-carry-frozen-response',ROOT/'docs/phase_carry_frozen_response_assessment.json')
contexts,trials,factorials=s.design(old['selections'],Path('/tmp/shared-read-only-nonexistent-root'))
assert not Path('/tmp/shared-read-only-nonexistent-root').exists()
saved={}
for context in old['contexts']:
    with np.load(context['source'],allow_pickle=False) as z:saved[context['key']]={k:z[k].copy() for k in z.files}
unique=set();maxerr=0.;checks=[]
for context in contexts:
    chemical=saved[context['chemical_context']];physical=saved[context['physical_context']]
    values,accounting=s.transferred_context(chemical,physical)
    for name in ('initial','ids'):assert np.array_equal(values[name],chemical[name])
    for name in ('delta','volumes'):assert np.array_equal(values[name],physical[name])
    for trial in [t for t in trials if t['context']==context['key']]:
        changed=values['initial'].copy();factor=1.;error=0.
        if trial['target'] is not None:
            i=values['ids'].tolist().index(trial['target'])
            pair=[saved[f'seed-{trial["seed"]}_{name}']['initial'][0,i] for name in ('sham',f'redistributed-cell-{trial["target"]}')]
            changed,record=s.shared_pulse(values['initial'],values['volumes'],values['ids'],trial['target'],trial['dose'],pair)
            expected=values['initial'].copy()
            factor=.9 if trial['dose']=='fractional' else 1-.1*min(pair)/float(expected[0,i])
            expected[0,i]*=factor
            assert np.array_equal(changed,expected) and record['factor']==factor
            error=record['pulse_amount_relative_error'];maxerr=max(maxerr,error)
        unique.add((context['key'],*[(a.dtype.str,a.shape,a.tobytes()) for a in (values['delta'],values['volumes'],values['ids'],changed)]))
        checks.append(dict(trial=trial['key'],factor=float(factor),literal_array_equality=True,pulse_amount_relative_error=error))
assert (len(saved),len(contexts),len(trials),len(factorials),len(unique))==(9,21,69,6,57)
assert maxerr<=2e-13 and digest(source)==before
functions={name:hashlib.sha256(inspect.getsource(getattr(s,name)).encode()).hexdigest()
    for name in ('accepted_parent','design','transferred_context','shared_pulse')}
report=dict(passed=True,reviewed_utc=datetime.now(timezone.utc).isoformat(),parent_protocol_sha256=digest(ROOT/'output/studies/phase-carry-frozen-response/protocol.json'),
    source_path=str(source),source_sha256=before,helper_source_sha256=functions,
    actual_parent_contexts=9,shared_contexts=21,logical_trials=69,factorials=6,
    unique_initial_value_problems=57,max_pulse_amount_relative_error=maxerr,
    all_arrays_match_literal_source_or_recipient_only_pulse=True,checks=checks,
    phases_loaded_or_copied=False,integrations_performed=0,output_files_written=0,
    script_path=str(Path(__file__)),script_sha256=digest(__file__),
    scope='Calls read-only accepted-parent/design/assignment/pulse helpers; validates against literal independent NumPy source selection and recipient-only canonical dosing. No prepare/run/production or new scientific evidence.')
REPORT.write_text(json.dumps(report,indent=2)+'\n')
print('Preparation helpers passed',REPORT,digest(REPORT),'runnerSHA',before)
