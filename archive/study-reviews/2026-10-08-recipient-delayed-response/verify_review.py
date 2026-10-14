"""Verify the completed delayed assay and unchanged accepted parents."""
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path

REPO=Path(__file__).resolve().parents[3]
ARCHIVE=Path(__file__).resolve().parent
STUDY=REPO/'output/studies/phase-carry-recipient-delayed-response'
SNAPSHOT=REPO/'archive/study-document-snapshots/2026-10-08-before-recipient-delayed-response-launch/manifest.json'


def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def verify(pins):
    for name,expected in pins.items():
        path=Path(name)
        if not path.is_absolute():path=REPO/path
        assert sha(path)==expected,f'Changed immutable evidence: {name}'


p=read(STUDY/'protocol.json');r=read(ARCHIVE/'independent_review.json')
a=read(REPO/'docs/phase_carry_recipient_delayed_response_assessment.json');snapshot=read(SNAPSHOT)
assert r['passed'] and r['inputs_unchanged'] and a['full_study_accepted']
assert read(STUDY/'status.json')['state']=='completed' and read(STUDY/'summary.json')['passed']
for pins in (p['source_sha256'],p['input_sha256'],r['inputs_sha256'],snapshot['unchanged_guards']):verify(pins)
parent_pins={name:value for name,value in snapshot['files_sha256'].items() if name not in ('README.md','docs/plan.md')}
verify(parent_pins)
for key,file in [('protocol_sha256','protocol.json'),('original_summary_sha256','summary.json'),('original_status_sha256','status.json')]:
    assert a[key]==sha(STUDY/file)
assert a['independent_review_sha256']==sha(ARCHIVE/'independent_review.json')
assert a['counts']==r['counts'] and len(r['scientific_comparisons'])==48
assert (len(p['source_sha256']),len(p['input_sha256']))==(64,369)
verify(read(ARCHIVE/'figure-provenance.json')['files_sha256'])
assets=[REPO/'README.md',REPO/'docs/plan.md',SNAPSHOT]
assets += [REPO/'docs'/name for name in ('phase_carry_recipient_delayed_response.md',
    'phase_carry_recipient_delayed_response_launch.md','phase_carry_recipient_delayed_response_launch_verification.json',
    'phase_carry_recipient_delayed_response_assessment.md','phase_carry_recipient_delayed_response_assessment.json')]
assets += [REPO/name for name in ('embryo/phase_carry_recipient_delayed_response.py','tests/test_phase_carry_recipient_delayed_response.py')]
assets += [x for x in ARCHIVE.iterdir() if x.is_file() and x.name!='review-provenance.json']
assets += [x for x in (REPO/'archive/study-launches/2026-10-08-recipient-delayed-response').iterdir() if x.is_file()]
result=dict(verified_utc=datetime.now(timezone.utc).isoformat(),passed=True,
    protocol_sha256=a['protocol_sha256'],original_summary_sha256=a['original_summary_sha256'],original_status_sha256=a['original_status_sha256'],
    independent_review_sha256=a['independent_review_sha256'],independent_review_script_sha256=sha(ARCHIVE/'independent_review.py'),
    source_and_input_pins_unchanged=True,source_pin_count=64,input_pin_count=369,immutable_consumed_artifacts=len(r['inputs_sha256']),
    ledger_manuscript_and_parent_unchanged=True,unchanged_guards={**snapshot['unchanged_guards'],**parent_pins},
    files_sha256={str(x.relative_to(REPO)):sha(x) for x in sorted(set(assets))},figure_visually_reviewed=True,
    scope='Completed fixed-delay response assay; reporting updated, accepted parents/ledger/manuscript unchanged, next amount-control experiment not launched.')
(ARCHIVE/'review-provenance.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('files_sha256','unchanged_guards')},indent=2))
