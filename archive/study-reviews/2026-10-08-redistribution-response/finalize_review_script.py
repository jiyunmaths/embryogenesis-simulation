"""Save a completed review after all independent audits; no scientific run."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import shutil

root = Path('/home/sji/Documents/Programming/Python/embryogenesis-simulation')
review = root / 'archive/study-reviews/2026-10-08-redistribution-response'
snapshot = root / 'archive/study-document-snapshots/2026-10-08-before-completed-redistribution-response-review'
read = lambda p: json.loads(Path(p).read_text())
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
p = read(root / 'outputs/phase-carry-redistribution-response/protocol.json')
s = read(root / 'outputs/phase-carry-redistribution-response/summary.json')
status = read(root / 'outputs/phase-carry-redistribution-response/status.json')
r = read(review / 'response_reconstruction.json')
n = read(review / 'numerical_verification.json')
provenance = read('/tmp/redistribution-response-provenance-verification.json')
assert provenance['passed'] and len(provenance['completed_paths']) == 66
assert provenance['prepared_paths_verified'] == len(provenance['context_gates']) == 66
assert len(provenance['paired_reports']) == 6 and len(provenance['pilot_quality']) == 3
assert provenance['source_hashes_verified'] == len(p['source_sha256']) == 53
assert provenance['input_hashes_verified'] == len(p['input_sha256']) == 2412
assert s['passed'] and status['state'] == 'completed' and status['passed']
assert r['independent_response_reconstruction_passed'] and n['summary']['all_reconstructed_gates_pass']
assert not n['failures'] and not n['recorded_metric_disagreements']
assert r['protocol_sha256'] == provenance['protocol_sha256'] == sha(root / 'outputs/phase-carry-redistribution-response/protocol.json')
guard = read(snapshot / 'manifest.json')
for file, digest in {**guard['files_sha256'], **guard['unchanged_ledger_and_manuscript_guards']}.items():
    assert sha(root / file) == digest, file

for src, name in [('/tmp/redistribution-response-provenance-verification.json', 'provenance_verification.json'),
                  ('/tmp/verify_redistribution_response_provenance.py', 'provenance_initial_script.py'),
                  ('/tmp/verify_redistribution_response_stages.py', 'provenance_verification_script.py'),
                  (__file__, 'finalize_review_script.py')]:
    dst = review / name
    assert not dst.exists(), dst
    shutil.copy2(src, dst)

compact = dict(reviewed_utc=datetime.now(timezone.utc).isoformat(),
    protocol_sha256=r['protocol_sha256'], full_study_accepted=True,
    original_summary_sha256=sha(root / 'outputs/phase-carry-redistribution-response/summary.json'),
    original_status_sha256=sha(root / 'outputs/phase-carry-redistribution-response/status.json'),
    independent_histories=3, new_histories=0, paths=66, new_paths=48, reused_paths=18,
    source_hashes_verified=53, input_hashes_verified=2412,
    runtime_binaries_verified=2, prepared_states_verified=66, context_gates_verified=66,
    immutable_paths_verified=provenance['immutable_files_rechecked'],
    pilot_reports=3, full_reports=3, full_mechanical_pairs=33, full_response_pairs=24,
    quality_extrema={k:n['summary'][k] for k in ('max_volume_error','min_radius','max_clipping','max_boundary','max_dilution_error')},
    minimum_whole_window_contrast=n['summary']['minimum_sampled_contrast'],
    full_timestep_maxima=n['summary']['mechanical_max'],
    max_normalized_response_error=r['max_full_network_normalized_response_error'],
    max_relative_auc_error=r['max_relative_target_auc_error'],
    max_recovery_time_error=r['max_recovery_time_discrepancy'],
    max_response_effect_error=r['max_coarse_fine_primary_effect_discrepancy'],
    max_pulse_amount_relative_error=r['max_relative_pulse_amount_error'],
    effect_screen=.01, history_outcomes=r['history_outcomes'], recipients=r['recipients'],
    all_target_recovered=True, all_network_recovered=True,
    target_recovery_interval=[2.25,4.35], network_recovery_interval=[2.25,5.4],
    numerical_audit_limitations=n['limitations'],
    historical_document_guard_scope=provenance['guard_mismatch_summary'],
    initial_amount_matching_time=450., challenge_time=510., horizon=570.,
    review_archive=str(review.relative_to(root)),
    pre_review_reporting_snapshot=str(snapshot.relative_to(root)),
    scope='Selective delayed responses after redistribution at matched initial t450 amounts in three existing mature histories. '
          'Challenge-time quantities/geometry differ. No new histories, divisions, autonomous/inherited identity, '
          'nearest-neighbor isolation, unique mechanical mediation, full behavioral dosage factorial, '
          'fresh carry zygote development, general spatial/transport validation, or wholesale ledger/manuscript revalidation.')
assessment_json = root / 'docs/phase_carry_redistribution_response_assessment.json'
assessment_md = root / 'docs/phase_carry_redistribution_response_assessment.md'
assert not assessment_json.exists() and not assessment_md.exists()
assessment_json.write_text(json.dumps(compact, indent=2) + '\n')
draft = Path('/tmp/redistribution-response-assessment-draft.md').read_text()
draft = draft.replace('The independent response and numerical reconstructions pass; final provenance acceptance remains pending the separate immutable-evidence audit.',
    'All 66 logical paths pass the unchanged numerical, physical and independently verified provenance gates.')
draft = draft.replace('Model time and chemical quantities remain dimensionless.',
    'The tested diffusivity ratio is D_b/D_a=27.5 and polarity-tension contrast chi=0.35 on a 72-cubed grid. '
    'Phase fields retain float64 update-residual carry and chemistry remains float64. Model time and chemical quantities remain dimensionless.')
placeholder = '**[Final audit placeholder: expected 53 scientific sources and 2,412 frozen inputs; insert final immutable-evidence path count and acceptance result after archive audit.]**'
assert placeholder in draft
draft = draft.replace(placeholder,
    f"The [independent provenance audit](../archive/study-reviews/2026-10-08-redistribution-response/provenance_verification.json) "
    f"verifies 53 scientific source hashes, 2,412 input hashes, two pinned runtime binaries, all 66 exact preparations/context gates, "
    f"and {provenance['immutable_files_rechecked']:,} immutable paths before and after review. "
    'The user-authorized experimental-literature manuscript revision after launch is recorded separately from the frozen scientific pins; '
    'it does not change this assay\'s equations or outcomes. No manuscript file is a frozen source/input pin for this assay.')
draft += ('\n\nThe [compact assessment](phase_carry_redistribution_response_assessment.json), '
          '[six-recipient CSV](../archive/study-reviews/2026-10-08-redistribution-response/recipient-summary.csv), '
          'and [pre-review reporting snapshot](../archive/study-document-snapshots/2026-10-08-before-completed-redistribution-response-review/manifest.json) '
          'preserve the decisions and reporting provenance. Scientific sources, kernels, tests, protocols, preparations, histories, '
          'checkpoints, prior assessments, ledger and manuscript remain unchanged by this review. README and the working plan point to '
          'the completed assessment; no new simulation is launched.\n')
assert 'placeholder' not in draft.lower() and 'acceptance remains pending' not in draft
assessment_md.write_text(draft)

def replace_once(text, old, new):
    assert text.count(old) == 1, old
    return text.replace(old, new)

readme = (root / 'README.md').read_text()
readme = replace_once(readme,
    'delayed responses after initially amount-preserving redistribution are now under test; acceptance remains pending.',
    'the completed redistribution-response assay now qualifies selective later response changes after an equal-initial-amount intervention at this point.')
readme = replace_once(readme,
    'responses from those new redistribution backgrounds are now being tested by the [active delayed-response assay](docs/phase_carry_redistribution_response_launch.md), with outcomes pending.',
    'the [completed redistribution-response assay](docs/phase_carry_redistribution_response_assessment.md) now detects the same selective history-level response outcomes after controlling original amounts.')
readme = replace_once(readme,
    'Delayed pulses from the initially amount-preserving backgrounds are now [underway](docs/phase_carry_redistribution_response_launch.md); their outcomes remain pending.',
    'The following [completed delayed-pulse assay](docs/phase_carry_redistribution_response_assessment.md) qualifies selective later response changes from these initially amount-preserving backgrounds.')
old = next(line for line in readme.splitlines() if line.startswith('**Active delayed redistribution-response test:**'))
new = ('**Completed delayed redistribution-response test:** the [independent assessment](docs/phase_carry_redistribution_response_assessment.md) '
       'accepts all 66 logical paths: 48 new and 18 exact reused paths within three existing mature histories, t=510–570. '
       'Response effects exceed the prespecified full-window 0.01 screen in 0/2 recipients in history 7 and 1/2 in histories 8 and 9 '
       '(cells 19 and 29), at both fractional and common-amount doses and both timesteps. All 48 pulse paths recover toward their own moving controls. '
       'Cells 20 and 30 retain about 16.3-fold and 17.0-fold lower baseline activator despite below-screen normalized response effects; '
       'chemical levels and response behavior are distinct attributes. The intervention matched original amounts at t=450; '
       'challenge-time amounts and physical states can differ. Maximum raw log discrepancy is 0.000378 and normalized response error 1.76e-5. '
       'This qualifies selective context-sensitive behavior at matched initial amounts, with zero new histories and no autonomous/inherited identity claim. '
       '[Frozen protocol](docs/phase_carry_redistribution_response.md), [historical launch](docs/phase_carry_redistribution_response_launch.md), '
       'and [result figure](archive/study-reviews/2026-10-08-redistribution-response/redistribution-response.png). '
       'Next isolate the contribution of ongoing mechanics with frozen t=510 response controls.')
readme = replace_once(readme, old, new)
(root / 'README.md').write_text(readme)

plan = (root / 'docs/plan.md').read_text()
old = next(line for line in plan.splitlines() if line.startswith('**Active redistribution-response assay:**'))
new = ('**Completed redistribution-response assay:** the [independent assessment](phase_carry_redistribution_response_assessment.md) '
       'accepts all 66 logical paths (48 new, 18 exact reused) and all context, pilot/full and physical gates. '
       'Above-screen delayed responses occur in 0/2, 1/2 and 1/2 recipients in histories 7, 8 and 9, respectively, at both doses and timesteps; '
       'all 48 pulse paths recover toward their own continuing controls. Initial species amounts were matched at t=450, not at the t=510 challenge. '
       'Cells 20/30 retain about 16.3/17.0-fold lower activator with below-screen normalized delayed effects. '
       'The supported result is selective context-sensitive chemical behavior in three existing mature histories, not new developmental replication '
       'or autonomous/inherited identity. [Frozen protocol](phase_carry_redistribution_response.md) and '
       '[historical launch](phase_carry_redistribution_response_launch.md).')
plan = replace_once(plan, old, new)
plan = replace_once(plan,
    'delayed responses from those initially amount-preserving backgrounds are now [under test](phase_carry_redistribution_response_launch.md), with outcomes pending.',
    'the [completed initially amount-preserving follow-up](phase_carry_redistribution_response_assessment.md) reproduces these selective delayed outcomes.')
old = next(line for line in plan.splitlines() if line.startswith('**Current priority: assess delayed responses'))
new = ('**Current priority: separate ongoing mechanical mediation of delayed responses.** '
       'The [accepted redistribution-response study](phase_carry_redistribution_response_assessment.md) shows that changing the original bulk amounts '
       'is not required for the delayed effects in cells 19 and 29 at this tested point. Next freeze each actual t=510 transport operator and measured volumes '
       'and repeat chemistry-only pulse/own-control comparisons using both original doses and independent chemical-solver checks. '
       'This tests whether continuing mechanical evolution during t=510–570 is necessary; each frozen background still inherits distinct chemistry and physical context. '
       'A shared operator and recipient-state matching are separate controls for those initial differences. Geometry replay additionally requires independently qualified '
       'source-time sampling, so the existing source-spacing failure must not be bypassed. Retain the [phase-update claim audit](results_ledger.md#phase-update-precision-audit), '
       'original failed checks, and the no-carry ledger/manuscript snapshot. These mature continuations do not qualify fresh carry-corrected zygote development '
       'or broadly revalidate older parameters, amplitudes, positive pulses, reservoir studies, or the 262.1-fold ratio. This review launches no new experiment.')
plan = replace_once(plan, old, new)
plan = replace_once(plan,
    'Preserve earlier failures and defer parameter-map expansion until the active redistribution-response test is assessed;',
    'Preserve earlier failures; the redistribution-response test is now accepted at its tested point, while wider parameter-map expansion remains separate;')
plan = replace_once(plan,
    'The [delayed redistribution-response test](phase_carry_redistribution_response_launch.md) is active with own controls and both pulse doses; scientific outcomes and mechanical mediation remain separate.',
    'The [completed redistribution-response test](phase_carry_redistribution_response_assessment.md) detects selective later effects in 0/2, 1/2 and 1/2 recipients across histories at both doses, '
    'with all pulses recovering. Original amount matching is controlled; challenge-time state and mechanical mediation remain separate questions.')
(root / 'docs/plan.md').write_text(plan)

for file, digest in guard['unchanged_ledger_and_manuscript_guards'].items():
    assert sha(root / file) == digest, file
files = {str(path.relative_to(root)): sha(path) for path in review.iterdir() if path.is_file()}
files.update({str(path.relative_to(root)):sha(path) for path in [assessment_json,assessment_md,root/'README.md',root/'docs/plan.md']})
manifest = dict(reviewed_utc=datetime.now(timezone.utc).isoformat(), passed=True,
    protocol_sha256=r['protocol_sha256'], original_summary_sha256=compact['original_summary_sha256'],
    original_status_sha256=compact['original_status_sha256'], files_sha256=files,
    unchanged_ledger_and_manuscript_guards=guard['unchanged_ledger_and_manuscript_guards'],
    scope='Independent completed-study review and reporting update only. Source/input/runtime/physical/response audits pass. '
          'No simulations, threshold changes, scientific-source or output writes, ledger refresh, or manuscript results promotion.')
(review / 'review-provenance.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps({'accepted':True,'histories':3,'paths':66,'counts':['0/2','1/2','1/2'],
                  'assessment':str(assessment_md),'review_archive':str(review)},indent=2))
