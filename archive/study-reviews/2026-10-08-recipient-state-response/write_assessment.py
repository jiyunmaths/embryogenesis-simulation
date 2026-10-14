"""Curate the accepted recipient-matching assay, preserving parent evidence."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
ARCHIVE = Path(__file__).resolve().parent
STUDY = REPO / 'output/studies/phase-carry-recipient-state-response-grouped'


def load(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


r = load(ARCHIVE/'independent_review.json')
p = load(STUDY/'protocol.json')
status = load(STUDY/'status.json')
summary = load(STUDY/'summary.json')
diagnostic = load(ARCHIVE/'reference-endpoint-diagnostic.json')
timing = load(STUDY/'timing.json')
assert r['passed'] and r['inputs_unchanged'] and summary['passed'] and status['passed']
assert status['state'] == 'completed' and sha(STUDY/'protocol.json') == r['protocol_sha256']
effects = {(x['seed'], x['cell'], x['matching_reference'], x['physics']): x
           for x in r['scientific_comparisons']}
baselines = {(x['seed'], x['cell'], x['matching_reference'], x['physics']): x
             for x in r['baseline_descriptions']}
cells = [(7, 22), (7, 23), (8, 19), (8, 20), (9, 29), (9, 30)]
effect_rows = '\n'.join('| '+f'{s} / {c}'+' | '+' | '.join(
    f"{effects[s,c,l,h]['full_effects']['fine']:.6f}" for l,h in [('S','S'),('S','R'),('R','S'),('R','R')])
    +' | '+('Yes' if all(effects[s,c,l,h]['above_screen'] for l in ('S','R') for h in ('S','R')) else 'No')+' |'
    for s,c in cells)
baseline_rows = '\n'.join('| '+f'{s} / {c}'+' | '+' | '.join(
    f"{baselines[s,c,'S',h]['final_species_ratios'][0]:.6f}" for h in ('S','R'))+' | '+' | '.join(
    f"{baselines[s,c,'S',h]['late_network_log_rms']:.6f}" for h in ('S','R'))+' |' for s,c in cells)
fine_metrics = [x for y in r['scientific_comparisons'] for x in y['fine_pulse_metrics'].values()]
all_recovered = all(x[k] is not None for x in fine_metrics for k in ('target_recovery_time', 'network_recovery_time'))
assert all_recovered
assessment = dict(
    reviewed_utc=datetime.now(timezone.utc).isoformat(), full_study_accepted=True,
    protocol_sha256=r['protocol_sha256'], original_summary_sha256=sha(STUDY/'summary.json'),
    original_status_sha256=sha(STUDY/'status.json'), counts=r['counts'],
    source_hashes_verified=len(p['source_sha256']), input_hashes_verified=len(p['input_sha256']),
    immutable_consumed_artifacts=len(r['inputs_sha256']), gate_counts=r['gate_counts'], maxima=r['maxima'],
    maximum_dose_relative_error=r['input_audit']['maximum_pulse_relative_amount_error'],
    fine_all_pulses_recovered=all_recovered,
    recovery_ranges={k:[min(x[k] for x in fine_metrics),max(x[k] for x in fine_metrics)]
                     for k in ('target_recovery_time','network_recovery_time')},
    history_outcomes=r['history_outcomes'], scientific_comparisons=r['scientific_comparisons'],
    descriptive_contrasts=r['descriptive_contrasts'], baseline_descriptions=r['baseline_descriptions'],
    endpoint_reference_diagnostic=dict(post_run_descriptive=True, no_new_acceptance_threshold=True,
        pairs=len(diagnostic['rows']),
        maximum_endpoint_log_discrepancy=diagnostic['maximum_endpoint_raw_log_reference_difference'],
        maximum_late_log_discrepancy=diagnostic['maximum_late_raw_log_reference_difference'],
        path=str(ARCHIVE/'reference-endpoint-diagnostic.json'),
        sha256=sha(ARCHIVE/'reference-endpoint-diagnostic.json')),
    recorded_metrics_compared=r['recorded_metrics_compared'],
    max_recorded_metric_difference=r['maximum_recorded_metric_difference'],
    native_continuity_pairs=len(r['native_continuity']), native_continuity_max_log_error=r['maxima']['native_log'],
    timing=dict(**timing, gpu_production_seconds=sum(x['wall_seconds'] for x in r['stage_checks'])),
    independent_review=str(ARCHIVE/'independent_review.json'), independent_review_sha256=sha(ARCHIVE/'independent_review.json'),
    review_archive=str(ARCHIVE.relative_to(REPO)),
    prelaunch_snapshot='archive/study-document-snapshots/2026-10-08-before-recipient-state-response-launch/manifest.json',
    scope=('Whole-surrounding-network chemical sensitivity despite identical initial recipient chemistry and shared physics. '
           'Immediate matching and pulse, freely evolving target and own control; unmatched surrounding amounts. '
           'Three existing mature histories, no new histories. Endpoint reference agreement is post-run descriptive. '
           'No settled response, autonomous identity, neighbor specificity, new development, ledger/manuscript promotion or next launch.'))
(REPO/'docs/phase_carry_recipient_state_response_assessment.json').write_text(json.dumps(assessment,indent=2)+'\n')
text = f'''# Recipient-state matching — accepted assessment

**Matching the recipient's initial chemistry does not remove surrounding-network effects.** Cells **19, 20, 29 and 30** clear the prespecified response screen using both measured recipient reference states and both shared physical contexts. The counts are **0/2, 2/2 and 2/2** recipients in histories 7, 8 and 9, respectively, in every reference/context combination. All numerical gates pass. These are immediate post-reset responses in three existing mature histories, not settled response traits or autonomous cell identities.

## What was controlled

The [prespecified protocol](phase_carry_recipient_state_response.md) asks whether surrounding chemistry changes a recipient's response when its own initial activator and inhibitor are equal. For each selected recipient, both surrounding starts use the same measured pair: either its sham pair or its redistributed pair. Every other cell's concentrations remain exactly as measured. The two chemical starts share the same saved conservative transport operator, matching volumes and cell-ID order. Repeat this on both measured physical contexts.

The matching and 10% activator depletion occur together at t=510. The target then evolves freely; each pulse is compared with its own matched, unpulsed continuation through t=570. There is no continuing clamp or pre-pulse equilibration. The target concentration, volume and removed activator amount are equal between surrounding starts within each reference/physical-context comparison. Fractional and common-amount depletion therefore coincide; there is one challenge design, not two independent doses. Global chemical amounts remain unmatched. Preparation records separately preserve physical-volume reassignment and local-matching amount changes for both species.

There are **45 distinct backgrounds and 93 distinct initial-value problems**: 45 controls and 48 target pulses. The sham/sham native control is shared across the two recipients within each history; all starts are calculated anew. Three timesteps produce **279 GPU paths**, **279 CPU checks**, and **186 DOP853/Radau reference trajectories**. Forty-five unchanged cases reproduce accepted parent paths at five numerical settings, providing 225 continuity comparisons. No parent path is counted as a new calculation or new history. Six selected recipients are nested within **three histories**, with **zero new developmental histories**.

Only Gierer–Meinhardt chemistry and conservative exchange advance. Transport and volumes are frozen: mechanics, polarity, geometric dilution and division are not advanced. The [qualified grouped GPU layout](phase_carry_recipient_state_response_grouped.md) uses independent chemical columns on nine measured physical blocks, with 42 duplicated positive-control padding columns excluded from reported paths. It changes no equation or shared scientific kernel.

## Response results

The primary effect is the full 60-unit time RMS of the difference between the target's signed, two-species pulse/own-control log responses, normalized by the actual absolute initial log pulse amplitude. The screen is strictly **greater than 0.01**. Each table entry holds initial target chemistry and physics equal between the two surrounding starts. Below-screen is not equivalence.

| History / cell | Sham match, sham physics | Sham match, redistributed physics | Redistributed match, sham physics | Redistributed match, redistributed physics | Above screen in all four |
| --- | ---: | ---: | ---: | ---: | --- |
{effect_rows}

Cells 19/29 retain detected differences after target matching. Cells 20/30, previously below screen with unmatched target chemistry, now have effects **0.0474–0.1006** under all four matching/context combinations. Holding the recipient's starting concentrations equal therefore does not eliminate the influence of the rest of the chemical network. Matching changes the subsequent trajectory; comparing old and new effects does not yield an additive percentage attributable to local versus surrounding chemistry.

All **48 fine pulse paths recover** toward their own controls. Target and volume-weighted network recovery each range from **0.45 to 5.40 units**, using the prespecified sustained-recovery definition and remaining observation requirement. Pulse recovery does not erase differences between the two chemical backgrounds.

The largest physical-context and signed interaction RMS values are **0.011955** and **0.009881**, respectively. They were prespecified as descriptive continuous measures with no scientific screen; the primary chemical screen must not be applied to them. Signed waveform identities hold to `1.11e-16`, but their scalar RMS magnitudes are not additive.

## Target states reappear, but immediate responses depend on the reset

Every matched control starts with redistributed/sham target ratios exactly **one** for both species. Later, differences reappear. The following fine results use the sham matching reference; the alternative reference returns to the same endpoints to numerical precision in the observed window.

| History / cell | Final activator R/S: sham physics | Final activator R/S: redistributed physics | Late network log RMS: sham physics | Late network log RMS: redistributed physics |
| --- | ---: | ---: | ---: | ---: |
{baseline_rows}

Cells 20/30 recover approximately **16.6–20.9-fold lower activator** under redistributed surrounding chemistry despite initially equal target states. The surrounding network can restore local differences in this model. Whole-network log separation during elapsed 36–60 units remains approximately **0.6398–0.9962**. These are descriptive finite-horizon state differences, not an equilibrium or biological identity acceptance test.

A separate [post-run endpoint diagnostic](../archive/study-reviews/2026-10-08-recipient-state-response/reference-endpoint-diagnostic.json) compares the two matching references at fixed surrounding chemistry and physics. Across 24 unpulsed whole-network pairs, maximum final log discrepancy is `2.61e-13`, and maximum discrepancy in the last 24 units is `6.53e-11`. This diagnostic integrates no new trajectories and introduces no acceptance threshold. It shows agreement for these measured references and this horizon; it establishes neither unique attractors nor general equivalence.

The immediate pulse effects differ between the two matching references even though their unpulsed endpoints agree. Cells 20/30 show their large response differences during early reset relaxation, as the [figure](../archive/study-reviews/2026-10-08-recipient-state-response/recipient-state-response.png) illustrates. Thus initial matching reveals surrounding-network sensitivity, but it does not yet establish different response behavior after the controls have settled. Delaying the pulse is the next necessary distinction.

## Numerical qualification

| Saved-array check | Maximum / range | Prespecified limit |
| --- | ---: | ---: |
| DOP853/Radau raw log discrepancy | `6.42e-12` | `1e-8` |
| DOP853/Radau normalized pulse discrepancy | `5.91e-11` | `1e-8` |
| DOP853/Radau signed contrast discrepancy | `2.79e-11` | `1e-8` |
| Production/reference raw log discrepancy | `7.25e-6` | `0.01` |
| CPU/PyTorch raw log discrepancy | `4.96e-14` | `1e-8` |
| Production normalized pulse discrepancy | `1.14e-5` | `0.01` |
| Production signed contrast discrepancy | `4.07e-6` | `0.01` |
| RMS contrast discrepancy | `1.66e-7` | `0.01` |
| Relative activator-response area discrepancy | `7.47e-6` | `0.02` |
| Recovery-time discrepancy | `0` | `0.3` |
| Observed order, all 93 resolved cases | `2.00056–2.00699` | `1.7–2.3` |
| Relative requested/measured pulse amount error | `9.99e-16` | `2e-13` |

All 93 raw dual-reference checks, 279 production accuracy checks, 186 raw timestep comparisons, 432 pulse comparison rows and 540 contrast comparison rows pass. Primary decisions agree at all three timesteps and both reference solvers. Saved states are positive; padding and checkpoint histories agree exactly. Internal Euler-stage checks are unchanged runtime/source evidence, not independently reconstructed unsaved observations.

A separate [raw-array audit](../archive/study-reviews/2026-10-08-recipient-state-response/independent_review.json), importing no project analysis helpers, reconstructs **6,357 recorded quantities** to maximum discrepancy `2.22e-16`. All **62 source pins, 250 input pins and 1,162 consumed artifacts** verify unchanged. The 225 native case/method comparisons have zero log discrepancy. The [launch evidence](phase_carry_recipient_state_response_launch.md) retains 13 design tests, three grouping tests, literal prepared-input checks and all actual CUDA prefixes. The initial grouping-test missing-import failure is retained; it was fixed before preparation, without changing scientific criteria.

GPU production took **165.82 seconds**. Full coordinator elapsed time was **287.80 seconds** (about 4.8 minutes), with four single-thread CPU verification workers overlapping GPU production. The grouped layout was **4.12 times faster in the short prefix benchmark**; no full ungrouped timing experiment was run. Worker-time sums are not elapsed wall time.

The [PNG](../archive/study-reviews/2026-10-08-recipient-state-response/recipient-state-response.png), [SVG](../archive/study-reviews/2026-10-08-recipient-state-response/recipient-state-response.svg), [24-comparison CSV](../archive/study-reviews/2026-10-08-recipient-state-response/recipient-summary.csv), [compact assessment](phase_carry_recipient_state_response_assessment.json) and review provenance preserve the reported result. Raw chemical trajectories remain in `output/studies/phase-carry-recipient-state-response-grouped`; no 3D fields were copied. The earlier ungrouped preparation remains unlaunched and unchanged.

## Next priority and limits

**Delay the challenge after target matching**, keeping each evolving control and the same crossed chemical/physical contexts. Compare immediate and delayed responses to determine which differences describe reset relaxation and which remain after local chemical states have re-established. This follow-up has not been launched. A subsequent initial-amount-matching control would separate surrounding allocation from total network dosage.

The present intervention changes the chemistry of the other fifteen cells as a whole and leaves their global amounts unmatched. It establishes neither nearest-neighbor specificity nor dosage-independent allocation effects. Retain the limits of three selected mature histories, historical no-carry developmental inheritance, sixteen-cell resolution, frozen mechanics and the inherited conductance approximation. No autonomous or inherited biological identity, fresh carry-corrected zygote differentiation, population prevalence or universal polarity mechanism is established. The results ledger, manuscript, accepted parent, baseline simulator and shared kernels remain unchanged; no subsequent experiment is launched here.
'''
(REPO/'docs/phase_carry_recipient_state_response_assessment.md').write_text(text)

readme = REPO/'README.md'
lines = readme.read_text().splitlines()
for i, line in enumerate(lines):
    if line.startswith('**Completed shared physical-context test:**'):
        lines[i] = line.replace('with recipient-state matching needed next to test surrounding-network effects.',
                               'with the recipient-state matching follow-up below testing surrounding-network effects.')
    if line.startswith('**Recipient-state matching is running:**'):
        lines[i] = ('**Completed recipient-state matching:** the [assessment](docs/phase_carry_recipient_state_response_assessment.md) '
            'accepts 279 GPU paths, 279 CPU checks and 186 reference trajectories. With initial target activator/inhibitor equal '
            'and physics shared, cells 19/20/29/30 retain above-screen response differences under both measured matching references '
            'and both physical contexts. Counts are 0/2, 2/2 and 2/2 across the three existing histories; all 48 fine pulses recover. '
            'Surrounding chemistry restores local concentration differences after matching, including about 16.6–20.9-fold lower '
            'activator in cells 20/30. The two local matching references return to the same contextual endpoints, but their immediate '
            'responses differ: this is surrounding-network sensitivity during a reset, not yet settled response behavior or autonomous '
            'identity. Surrounding amounts remain unmatched. All numerical gates, 62 source pins, 250 input pins and 1,162 consumed '
            'artifacts verify. The [grouped GPU layout](docs/phase_carry_recipient_state_response_grouped.md) changes no equations; '
            'the full study took about 4.8 minutes with concurrent CPU checks. **Next: delay the pulse after matching** to distinguish '
            'reset relaxation from settled behavior; the follow-up is not launched. [Figure](archive/study-reviews/2026-10-08-recipient-state-response/recipient-state-response.png).')
readme.write_text('\n'.join(lines)+'\n')

plan = REPO/'docs/plan.md'
lines = plan.read_text().splitlines()
for i, line in enumerate(lines):
    if line.startswith('**Current priority: qualify recipient-state matching.**'):
        lines[i] = ('**Completed recipient-state matching:** the [accepted assessment](phase_carry_recipient_state_response_assessment.md) '
            'qualifies all 279 GPU paths, 279 CPU checks and 186 references. Cells 19/20/29/30 clear the primary screen under both '
            'matching references and both physical contexts; counts are 0/2, 2/2 and 2/2 across three existing mature histories. '
            'Every fine pulse recovers and target-state differences reappear after initially equal local chemistry. The two matching '
            'references return to the same contextual endpoints, while immediate response amplitudes differ. Global surrounding '
            'amounts remain unmatched. All source/input pins and separate saved-array reconstructions pass; no ledger or manuscript promotion.')
        lines.insert(i+1, '')
        lines.insert(i+2, ('**Current priority: distinguish reset responses from settled responses.** Plan a delayed challenge after '
            'recipient matching, with each arm’s own continuing control and the same shared physical contexts. Retain both matching '
            'references and all selected recipients; compare immediate and delayed effects without choosing a favorable delay or '
            'dropping failed numerical gates. Freeze a separate protocol before launch. This is not yet launched. Follow it with '
            'a surrounding-amount control to separate allocation from global dosage; moving mechanics, new developmental histories '
            'and autonomy require later experiments.'))
        break
content = '\n'.join(lines)+'\n'
content = content.replace('while recipient-state matching remains necessary to separate local chemical state from surrounding-network sensitivity.',
    'and the [accepted recipient-state matching](phase_carry_recipient_state_response_assessment.md) detects surrounding-network '
    'response effects despite equal initial local chemistry, in 0/2, 2/2 and 2/2 recipients under both references/contexts. '
    'Those immediate reset responses still require delayed challenges before being interpreted as settled behavior.')
plan.write_text(content)
print('Wrote assessment Markdown/JSON and updated README/plan reporting only.')
