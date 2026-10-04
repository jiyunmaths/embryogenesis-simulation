import json
import shutil

import pytest

from embryo import precision_claim_audit as audit


@pytest.fixture(autouse=True)
def isolated_evidence_repository(tmp_path, monkeypatch):
    """Exercise provenance checks without Git-ignored scientific run outputs.

    Preserve the real archived claim registry, but give its evidence references
    isolated synthetic bytes. Actual scientific evidence is checked separately
    by the ledger verification commands in the study workspace.
    """
    baseline = json.loads((audit.REPO/audit.ORIGINAL).read_text())
    for field in ('evidence_sha256', 'history_checkpoint_sha256'):
        for name in baseline[field]:
            target = tmp_path/name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text('Synthetic unit-test evidence for '+name+'\n')
            baseline[field][name] = audit.digest(target)
    for name in ('embryo/results_ledger.py', 'embryo/precision_claim_audit.py'):
        target = tmp_path/name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(audit.REPO/name, target)
    baseline['source_sha256']['embryo/results_ledger.py'] = audit.digest(tmp_path/'embryo/results_ledger.py')
    target = tmp_path/audit.ORIGINAL
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(baseline))
    fixtures = [
        dict(passed=True, independent_histories=1,
             same_dt_interventions=dict(fine=dict(errors=dict(chemical_log_max=.067)),
                                        finer=dict(errors=dict(chemical_log_max=.104)))),
        dict(passed=True, independent_histories=1), {}, {},
    ]
    for name, payload in zip(audit.AUDIT_EVIDENCE, fixtures):
        target = tmp_path/name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload))
    monkeypatch.setattr(audit, 'REPO', tmp_path)


def original():
    return json.loads((audit.REPO/audit.ORIGINAL).read_text())


def test_annotation_preserves_every_original_decision_count_and_evidence_hash():
    baseline = original()
    result = audit.annotate(baseline)
    assert [{k:v for k,v in r.items() if k != 'precision_audit'}
            for r in result['entries']] == baseline['entries']
    assert result['evidence_sha256'] == baseline['evidence_sha256']
    assert result['history_checkpoint_sha256'] == baseline['history_checkpoint_sha256']
    assert result['unique_primary_histories'] == 3
    assert next(r for r in result['entries'] if r['id']=='historical-development-refinement')['validation'] == 'fail'
    audit.verify(result)


def test_an_unclassified_new_study_blocks_annotation_instead_of_getting_accepted():
    ledger = original()
    ledger['entries'].append(dict(id='unreviewed-new-study'))
    with pytest.raises(ValueError, match='every curated'):
        audit.annotate(ledger)


def test_changed_original_acceptance_is_rejected_even_if_audit_metadata_is_unchanged():
    ledger = audit.annotate(original())
    row = next(r for r in ledger['entries'] if r['id']=='historical-development-refinement')
    row['validation'] = 'pass'
    with pytest.raises(ValueError, match='original claims'):
        audit.verify(ledger)


def test_pending_claim_cannot_be_promoted_by_changing_its_annotation():
    ledger = audit.annotate(original())
    row = next(r for r in ledger['entries'] if r['id']=='frozen-exchange-response')
    assert row['precision_audit']['dependency'] == 'frozen_inherited_no_carry'
    row['precision_audit']['carry_status'] = 'validated'
    with pytest.raises(ValueError, match='Entry annotation differs'):
        audit.verify(ledger)


def test_standalone_reference_is_not_invalidated_by_moving_phase_rounding():
    ledger = audit.annotate(original())
    row = next(r for r in ledger['entries'] if r['id']=='conservative-reference')
    assert row['precision_audit']['carry_status'] == 'not_applicable_to_reference_test'
    assert row['validation'] == next(r for r in original()['entries'] if r['id']==row['id'])['validation']


def test_single_history_carry_support_cannot_be_promoted_to_blanket_validation():
    ledger = audit.annotate(original())
    ledger['precision_audit']['qualified_carry_scope'] = 'All histories and assays validated.'
    with pytest.raises(ValueError, match='blanket validation'):
        audit.verify(ledger)
