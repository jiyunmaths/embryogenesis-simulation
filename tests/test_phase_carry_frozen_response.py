"""Qualify exact fixed operators, actual pulse normalization and resumable batches."""
import json
from pathlib import Path

import numpy as np
import pytest
import torch

from embryo import phase_carry_frozen_response as study
from embryo.feedback_long import digest
from embryo.gpu_backend import gm_step


def fixture(tmp_path, duration=.6):
    contexts = []
    initial = np.array([[.8, 1.1], [.9, 1.2]], dtype=np.float64)
    delta = np.array([[-.1, .1], [.05, -.05]])
    volumes, ids = np.array([1., 2.]), np.array([11, 17])
    for seed in study.HISTORIES:
        for recipient in (None, 11, 17):
            name = 'sham' if recipient is None else f'redistributed-cell-{recipient}'
            source = tmp_path/'contexts'/f'{seed}-{name}.npz'
            values = initial.copy()
            if recipient is not None:
                values[0, ids.tolist().index(recipient)] *= .8
            study.save_npz(source, delta=delta, volumes=volumes, initial=values, ids=ids)
            contexts.append(dict(key=f'seed-{seed}_{name}', seed=seed, background=name,
                targets=[11, 17] if recipient is None else [recipient], source=str(source)))
    trials = study.trial_design(contexts)
    for trial in trials:
        context = next(c for c in contexts if c['key'] == trial['context'])
        with np.load(context['source']) as z:
            values = z['initial'].copy()
        factor = 1. if trial['target'] is None else .9
        if trial['target'] is not None:
            values[0, ids.tolist().index(trial['target'])] *= factor
        trial.update(source=str(tmp_path/'starts'/(trial['key']+'.npz')), factor=factor,
            initial_log_amplitude=float(abs(np.log(factor))))
        study.save_npz(trial['source'], initial=values)
    p = dict(contexts=contexts, trials=trials, histories=list(study.HISTORIES),
        selections={str(s):dict(selected_ids=[11, 17]) for s in study.HISTORIES},
        levels=[dict(level=l, dt=dt) for l, dt in study.LEVELS], beta=2., da=.02, db=.55,
        duration=duration, interval=.15, checkpoint_interval=.3, effect_min=.01,
        late_window=24., criteria=study.CRITERIA, device='cpu',
        reference_settings=dict(rtol=1e-13, atol=1e-15, max_step=.15),
        runtime=dict(torch=str(torch.__version__), torch_cuda=torch.version.cuda, expected_gpu_name='unused'),
        source_sha256={}, input_sha256=study.pin([*[c['source'] for c in contexts], *[t['source'] for t in trials]]), scope='CPU test fixture')
    study.write_json(tmp_path/'protocol.json', p)
    return p


def test_exact_design_has_nine_contexts_and_33_logical_starts(tmp_path):
    p = fixture(tmp_path)
    assert len(p['contexts']) == 9
    assert len(p['trials']) == len({t['key'] for t in p['trials']}) == 33
    assert sum(t['target'] is None for t in p['trials']) == 9
    assert sum(t['dose'] == 'fractional' for t in p['trials']) == 12
    assert sum(t['dose'] == 'amount' for t in p['trials']) == 12
    assert len(p['levels'])*len(p['trials']) == 99
    assert study.unique_problem_count(p) == 21


def test_batch_has_exact_operators_and_positive_control_padding(tmp_path):
    p = fixture(tmp_path)
    initial, matrix, mapping, count = study.batch_layout(p)
    assert initial.shape == (2, 18, 5) and matrix.shape == (18, 18)
    assert len(mapping) == 33 and count == 2
    assert np.all(initial > 0)
    for ci, context in enumerate(p['contexts']):
        with np.load(context['source']) as z:
            np.testing.assert_array_equal(matrix[2*ci:2*ci+2, 2*ci:2*ci+2], z['delta'])
            if context['background'] != 'sham':
                np.testing.assert_array_equal(initial[:, 2*ci:2*ci+2, 3], z['initial'])
                np.testing.assert_array_equal(initial[:, 2*ci:2*ci+2, 4], z['initial'])
    off = matrix.copy()
    for ci in range(9):
        off[2*ci:2*ci+2, 2*ci:2*ci+2] = 0
    assert np.count_nonzero(off) == 0


def test_batched_unchanged_gm_matches_every_independent_trial(tmp_path):
    p = fixture(tmp_path)
    initial, delta, mapping, count = study.batch_layout(p)
    values = [torch.from_numpy(x.copy()) for x in (initial[0], initial[1], delta)]
    a, b = gm_step(*values, .00375, 2., .02, .55)
    actual = study.unpack(np.stack([a.numpy(), b.numpy()]), mapping, count)
    for trial, result in zip(p['trials'], actual):
        _, matrix, _, state, _ = study.load_context(p, trial)
        aa, bb = gm_step(torch.from_numpy(state[0]), torch.from_numpy(state[1]),
            torch.from_numpy(matrix), .00375, 2., .02, .55)
        np.testing.assert_allclose(result, np.stack([aa.numpy(), bb.numpy()]), rtol=0, atol=1e-15)


def test_one_context_cannot_communicate_with_another(tmp_path):
    p = fixture(tmp_path)
    initial, delta, mapping, count = study.batch_layout(p)
    first = [torch.from_numpy(x.copy()) for x in (initial[0], initial[1], delta)]
    other = initial.copy(); other[:, :2] *= 3.
    second = [torch.from_numpy(x.copy()) for x in (other[0], other[1], delta)]
    a = gm_step(*first, .00375, 2., .02, .55)
    b = gm_step(*second, .00375, 2., .02, .55)
    for x, y in zip(a, b):
        torch.testing.assert_close(x[2:], y[2:], rtol=0, atol=0)


@pytest.mark.parametrize('invalid', [np.array([[0., 1.], [1., 1.]]),
    np.array([[np.nan, 1.], [1., 1.]]), np.ones((2, 2), dtype=np.float32)])
def test_frozen_context_rejects_invalid_or_nonfloat64_chemistry(invalid):
    with pytest.raises(ValueError):
        study.validate_context([[-1., 1.], [.5, -.5]], [1., 2.], invalid, [11, 17])


def test_fixed_operator_is_preserved_without_symmetrizing():
    delta = np.array([[-.1, .1], [.05, -.05]])
    result = study.fixed_transport(delta, np.array([1., 2.]))
    np.testing.assert_array_equal(result.delta.toarray(), delta)
    np.testing.assert_array_equal(result.volumes, [1., 2.])
    changed = delta.copy(); changed[0, 1] += .0001
    with pytest.raises(ValueError):
        study.fixed_transport(changed, np.array([1., 2.]))


def test_modified_prepared_operator_is_rejected_by_input_hash(tmp_path):
    p = fixture(tmp_path); file = Path(p['contexts'][0]['source'])
    with np.load(file) as z:
        changed = {k:z[k].copy() for k in z.files}
    changed['delta'] *= 2.
    study.save_npz(file, **changed)
    with pytest.raises(ValueError, match='preparation'):
        study.load_context(p, p['trials'][0])


def test_stage_checks_use_actual_unchanged_gm_and_reject_invalid(tmp_path):
    p = fixture(tmp_path); state, delta, _, _ = study.batch_layout(p)
    a, b, matrix = [torch.from_numpy(x.copy()) for x in (state[0], state[1], delta)]
    expected = gm_step(a, b, matrix, .00375, 2., .02, .55)
    actual = study.checked_torch_step(a, b, matrix, .00375, p, count=1)
    for first, second in zip(expected, actual):
        torch.testing.assert_close(first, second, rtol=0, atol=0)
    b[0, 0] = 0.
    with pytest.raises(FloatingPointError):
        study.checked_torch_step(a, b, matrix, .00375, p, count=1)


def test_pulse_is_exactly_recipient_only_and_amount_accounted():
    initial = np.array([[.1, 2.], [.3, 4.]])
    changed = initial.copy(); changed[0, 0] *= .9
    removed = float((initial[0, 0]-changed[0, 0])*2.)
    record = dict(factor=.9, dose='amount', requested_removed_amount=.02,
        removed_activator_amount=removed, initial_target_activator_amount=.2,
        initial_log_amplitude=float(abs(np.log(.9))))
    assert study.pulse_accounting(initial, changed, np.array([2., 3.]), np.array([11, 17]), 11, record, .02) < 2e-13
    changed[1, 1] *= .9
    with pytest.raises(ValueError, match='recipient-only'):
        study.pulse_accounting(initial, changed, np.array([2., 3.]), np.array([11, 17]), 11, record, .02)


def test_full_cpu_pipeline_preserves_all_paths_and_reference_response_gates(tmp_path):
    p = fixture(tmp_path)
    for trial in p['trials']:
        assert study.reference_worker((tmp_path, trial))['passed']
    assert study.reference_qualification(tmp_path, p)['passed']
    for level in p['levels']:
        study.production(tmp_path, p, level)
        for trial in p['trials']:
            study.cpu_worker((tmp_path, trial, level))
    result = study.assess(tmp_path)
    assert result['passed']
    assert len(result['production_accuracy']) == 99
    assert len(result['orders']) == 33
    assert len(result['raw_timestep_comparisons']) == 66
    assert all(r['passed'] for r in result['response_comparisons'])
    assert result['primary_paths'] == 66 and result['diagnostic_paths'] == 33


def test_interrupted_batch_resumes_exact_committed_state(tmp_path, monkeypatch):
    p = fixture(tmp_path); original = study.checked_torch_step; calls = 0
    def interrupted(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 110:
            raise RuntimeError('Deliberate interruption')
        return original(*args, **kwargs)
    monkeypatch.setattr(study, 'checked_torch_step', interrupted)
    with pytest.raises(RuntimeError, match='interruption'):
        study.production(tmp_path, p, p['levels'][0])
    with np.load(tmp_path/'production/coarse/checkpoint.npz') as z:
        assert json.loads(str(z['metadata']))['step'] == 80
    monkeypatch.setattr(study, 'checked_torch_step', original)
    study.production(tmp_path, p, p['levels'][0])
    with np.load(tmp_path/'production/coarse/paths.npz') as z:
        resumed = z['trajectory'].copy()
    clean = tmp_path/'clean'; clean.mkdir(); other = fixture(clean)
    study.production(clean, other, other['levels'][0])
    with np.load(clean/'production/coarse/paths.npz') as z:
        np.testing.assert_array_equal(resumed, z['trajectory'])


def test_changed_checkpoint_is_not_silently_restarted(tmp_path, monkeypatch):
    p = fixture(tmp_path); original = study.checked_torch_step
    monkeypatch.setattr(study, 'checked_torch_step', lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError('interrupt')))
    with pytest.raises(RuntimeError):
        study.production(tmp_path, p, p['levels'][0])
    cp = tmp_path/'production/coarse/checkpoint.npz'
    with np.load(cp) as z:
        values = {k:z[k].copy() for k in z.files}
    values['state'].flat[0] *= 1.01
    study.save_npz(cp, **values)
    monkeypatch.setattr(study, 'checked_torch_step', original)
    with pytest.raises(ValueError, match='checkpoint'):
        study.production(tmp_path, p, p['levels'][0])


def test_changed_completed_production_cannot_authorize_cache_reuse(tmp_path):
    p = fixture(tmp_path); level = p['levels'][0]
    study.production(tmp_path, p, level)
    cp = tmp_path/'production/coarse/checkpoint.npz'
    with np.load(cp) as z:
        values = {k:z[k].copy() for k in z.files}
    values['state'].flat[0] *= 1.01
    study.save_npz(cp, **values)
    with pytest.raises(ValueError, match='evidence'):
        study.production(tmp_path, p, level)


def test_wrong_runtime_never_falls_back_to_cpu(tmp_path):
    p = fixture(tmp_path); p['runtime']['torch'] = 'changed'
    with pytest.raises(ValueError, match='runtime'):
        study.runtime_check(p)
    p['runtime']['torch'] = str(torch.__version__); p['device'] = 'invalid:0'
    with pytest.raises((ValueError, RuntimeError)):
        study.runtime_check(p)


def test_small_pulse_raw_reference_agreement_does_not_hide_normalized_error(tmp_path):
    p = fixture(tmp_path, duration=60.)
    times = np.arange(401)*.15; paths = {}
    for trial in p['trials']:
        _, _, _, initial, _ = study.load_context(p, trial)
        path = np.repeat(initial[None], len(times), axis=0)
        paths['dop853', trial['key']] = path
        paths['radau', trial['key']] = path.copy()
    trial = next(t for t in p['trials'] if t['target'] == 11 and t['dose'] == 'amount')
    trial['factor'] = .99999
    # Only this pulse is deliberately weakened, identically for both solvers.
    control = next(t for t in p['trials'] if t['context'] == trial['context'] and t['target'] is None)
    for method in ('dop853', 'radau'):
        paths[method, trial['key']] = paths[method, control['key']].copy()
        paths[method, trial['key']][:, 0, 0] *= trial['factor']
    paths['radau', trial['key']][4, 0, 0] *= np.exp(1e-10)
    assert abs(np.log(paths['radau', trial['key']]/paths['dop853', trial['key']])).max() < 1e-8
    summaries = {m:study.summarize_responses(paths, m, p, times) for m in ('dop853', 'radau')}
    report = study.compare_responses(paths, 'dop853', 'radau', p, times, summaries, limit=1e-8)
    assert not report['passed']
    failed = next(r for r in report['pulse_rows'] if r['trial'] == trial['key'])
    assert failed['normalized_response_error'] > 1e-8


def accepted_fixture(tmp_path):
    parent = tmp_path/'parent'; parent.mkdir()
    archive = tmp_path/'archive'; archive.mkdir()
    source = tmp_path/'model.py'; source.write_text('unchanged model\n')
    scientific_plan = tmp_path/'science/plan.md'; scientific_plan.parent.mkdir(); scientific_plan.write_text('scientific plan\n')
    study.write_json(parent/'protocol.json', dict(source_sha256=study.pin([source])))
    study.write_json(parent/'status.json', dict(state='completed'))
    study.write_json(parent/'summary.json', dict(passed=True))
    keys = dict(protocol_sha256=digest(parent/'protocol.json'),
        original_status_sha256=digest(parent/'status.json'), original_summary_sha256=digest(parent/'summary.json'))
    evidence_file = archive/'provenance_verification.json'
    study.write_json(evidence_file, dict(passed=True, protocol_sha256=keys['protocol_sha256'],
        all_mismatches=[], evidence_sha256={}))
    study.write_json(archive/'review-provenance.json', dict(passed=True, **keys,
        files_sha256={**study.pin([evidence_file, scientific_plan]),
            str(Path('README.md').resolve()):'historical reporting hash',
            str(Path('docs/plan.md').resolve()):'historical reporting hash'}))
    public = tmp_path/'assessment.json'
    study.write_json(public, dict(full_study_accepted=True, **keys, review_archive=str(archive)))
    return parent, public, source, scientific_plan


def test_compact_parent_review_skips_only_exact_reporting_paths(tmp_path):
    parent, public, source, scientific_plan = accepted_fixture(tmp_path)
    old, review, evidence, inputs = study.accepted_parent(parent, public)
    assert review['full_study_accepted'] and scientific_plan in inputs
    assert source.as_posix() in old['source_sha256']
    scientific_plan.write_text('changed model plan\n')
    with pytest.raises(ValueError, match='review artifact'):
        study.accepted_parent(parent, public)


@pytest.mark.parametrize('which', ['summary', 'source', 'public'])
def test_changed_or_unaccepted_parent_blocks_new_preparations(tmp_path, which):
    parent, public, source, _ = accepted_fixture(tmp_path)
    if which == 'summary':
        study.write_json(parent/'summary.json', dict(passed=False))
    elif which == 'source':
        source.write_text('changed model\n')
    else:
        doc = json.loads(public.read_text()); doc['full_study_accepted'] = False
        study.write_json(public, doc)
    with pytest.raises(ValueError):
        study.accepted_parent(parent, public)
