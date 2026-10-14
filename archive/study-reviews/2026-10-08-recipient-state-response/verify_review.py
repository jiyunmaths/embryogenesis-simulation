"""Verify immutable assay evidence and reporting assets after curation."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
ARCHIVE = Path(__file__).resolve().parent
STUDY = REPO/'output/studies/phase-carry-recipient-state-response-grouped'
SNAPSHOT = REPO/'archive/study-document-snapshots/2026-10-08-before-recipient-state-response-launch/manifest.json'


def load(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify(pins):
    for name, expected in pins.items():
        path = Path(name)
        if not path.is_absolute():
            path = REPO/path
        assert sha(path) == expected, f'Changed immutable evidence: {name}'


p = load(STUDY/'protocol.json')
r = load(ARCHIVE/'independent_review.json')
a = load(REPO/'docs/phase_carry_recipient_state_response_assessment.json')
snapshot = load(SNAPSHOT)
assert r['passed'] and r['inputs_unchanged'] and a['full_study_accepted']
assert load(STUDY/'status.json')['state'] == 'completed' and load(STUDY/'summary.json')['passed']
for pins in (p['source_sha256'], p['input_sha256'], r['inputs_sha256'], snapshot['unchanged_guards']):
    verify(pins)
parent_pins = {k:v for k,v in snapshot['files_sha256'].items() if k not in ('README.md','docs/plan.md')}
verify(parent_pins)
pilot = load(REPO/'output/studies/phase-carry-recipient-state-response/protocol.json')
verify(pilot['source_sha256']); verify(pilot['input_sha256'])
assert not (REPO/'output/studies/phase-carry-recipient-state-response/summary.json').exists()
for field, filename in [('protocol_sha256','protocol.json'),('original_summary_sha256','summary.json'),
                        ('original_status_sha256','status.json')]:
    assert a[field] == sha(STUDY/filename)
assert a['independent_review_sha256'] == sha(ARCHIVE/'independent_review.json')
assert a['counts'] == r['counts']
assert (len(p['source_sha256']),len(p['input_sha256']),len(r['inputs_sha256'])) == (62,250,1162)
assert len(r['scientific_comparisons']) == 24
for outcome in r['history_outcomes']:
    assert outcome['above_screen'] == {7:0,8:2,9:2}[outcome['seed']]
for row in r['baseline_descriptions']:
    assert row['initial_species_ratios'] == [1.,1.]
verify(load(ARCHIVE/'figure-provenance.json')['files_sha256'])
diagnostic = load(ARCHIVE/'reference-endpoint-diagnostic.json')
verify(diagnostic['inputs_sha256'])
assert diagnostic['inputs_unchanged'] and diagnostic['post_run_descriptive']

assets = [REPO/'README.md',REPO/'docs/plan.md',SNAPSHOT]
assets += [REPO/'docs'/name for name in (
    'phase_carry_recipient_state_response.md','phase_carry_recipient_state_response_grouped.md',
    'phase_carry_recipient_state_response_launch.md','phase_carry_recipient_state_response_launch_verification.json',
    'phase_carry_recipient_state_response_assessment.md','phase_carry_recipient_state_response_assessment.json')]
assets += [REPO/name for name in ('embryo/phase_carry_recipient_state_response.py',
    'embryo/phase_carry_recipient_state_response_grouped.py',
    'tests/test_phase_carry_recipient_state_response.py','tests/test_phase_carry_recipient_state_response_grouped.py')]
assets += [x for x in ARCHIVE.iterdir() if x.is_file() and x.name != 'review-provenance.json']
for folder in ('2026-10-08-recipient-state-response','2026-10-08-recipient-state-response-grouped'):
    assets += [x for x in (REPO/'archive/study-launches'/folder).iterdir() if x.is_file()]
files = {str(path.relative_to(REPO)):sha(path) for path in sorted(set(assets))}
report = dict(verified_utc=datetime.now(timezone.utc).isoformat(),passed=True,
    protocol_sha256=a['protocol_sha256'],original_summary_sha256=a['original_summary_sha256'],
    original_status_sha256=a['original_status_sha256'],
    independent_review_sha256=a['independent_review_sha256'],
    independent_review_script_sha256=sha(ARCHIVE/'independent_review.py'),
    source_and_input_pins_unchanged=True,source_pin_count=62,input_pin_count=250,
    immutable_consumed_artifacts=1162,ledger_manuscript_and_parent_unchanged=True,
    unchanged_guards={**snapshot['unchanged_guards'],**parent_pins},
    unlaunched_pilot_unchanged=True,files_sha256=files,
    figure_visually_reviewed=True,
    scope='Completed matched-recipient assay only; reporting updated, accepted parent and ledger/manuscript preserved; next delayed challenge not launched.')
(ARCHIVE/'review-provenance.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k not in ('files_sha256','unchanged_guards')},indent=2))
