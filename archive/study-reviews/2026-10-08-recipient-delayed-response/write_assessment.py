"""Curate the accepted fixed-delay assay and reporting documents."""
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path

REPO=Path(__file__).resolve().parents[3]
ARCHIVE=Path(__file__).resolve().parent
STUDY=REPO/'output/studies/phase-carry-recipient-delayed-response'


def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


r=read(ARCHIVE/'independent_review.json');p=read(STUDY/'protocol.json');summary=read(STUDY/'summary.json')
assert r['passed'] and r['inputs_unchanged'] and summary['passed'] and read(STUDY/'status.json')['state']=='completed'
assert r['protocol_sha256']==sha(STUDY/'protocol.json')
metrics=[m for x in r['scientific_comparisons'] for m in x['fine_pulse_metrics'].values()]
all_recovered=all(x[k] is not None for x in metrics for k in ('target_recovery_time','network_recovery_time'))
a=dict(reviewed_utc=datetime.now(timezone.utc).isoformat(),full_study_accepted=True,
    protocol_sha256=r['protocol_sha256'],original_summary_sha256=sha(STUDY/'summary.json'),original_status_sha256=sha(STUDY/'status.json'),
    counts=r['counts'],source_hashes_verified=len(p['source_sha256']),input_hashes_verified=len(p['input_sha256']),
    immutable_consumed_artifacts=len(r['inputs_sha256']),gate_counts=r['gate_counts'],maxima=r['maxima'],
    source_and_settling_checks=r['input_audit']['source_checks'],
    fine_all_pulses_recovered=all_recovered,recovery_ranges={k:[min(x[k] for x in metrics),max(x[k] for x in metrics)] for k in ('target_recovery_time','network_recovery_time')},
    history_outcomes=r['history_outcomes'],scientific_comparisons=r['scientific_comparisons'],
    immediate_vs_delayed=r['immediate_vs_delayed'],matching_reference_response_contrasts=r['matching_reference_response_contrasts'],
    baseline_descriptions=r['baseline_descriptions'],descriptive_contrasts=r['descriptive_contrasts'],
    maximum_matching_reference_normalized_difference=max(x['maximum_normalized_difference'] for x in r['matching_reference_response_contrasts']),
    maximum_matching_reference_response_rms=max(x['rms'] for x in r['matching_reference_response_contrasts']),
    maximum_saved_fine_control_raw_log_drift=max(x['maximum_raw_log_drift'] for x in r['control_drift']),
    maximum_dose_relative_error=r['input_audit']['maximum_pulse_relative_amount_error'],
    recorded_metrics_compared=r['recorded_metrics_compared'],max_recorded_metric_difference=r['maximum_recorded_metric_difference'],
    duplicate_continuity_comparisons=len(r['duplicate_continuity']),duplicate_continuity_max_log_error=r['maxima']['duplicate_log'],
    timing={**read(STUDY/'timing.json'),'gpu_production_seconds':sum(x['wall_seconds'] for x in r['stage_checks'])},
    independent_review=str(ARCHIVE/'independent_review.json'),independent_review_sha256=sha(ARCHIVE/'independent_review.json'),
    review_archive=str(ARCHIVE.relative_to(REPO)),prelaunch_snapshot='archive/study-document-snapshots/2026-10-08-before-recipient-delayed-response-launch/manifest.json',
    scope='Fixed60-unit delayed challenge after original matching, no new reset. Three existing histories; source settling verified, frozen physics. Unequal local chemistry and unmatched surrounding totals at challenge. No autonomous identity, new development, equivalence below screen or ledger/manuscript promotion.')
(REPO/'docs/phase_carry_recipient_delayed_response_assessment.json').write_text(json.dumps(a,indent=2)+'\n')
cells=[(7,22),(7,23),(8,19),(8,20),(9,29),(9,30)]
effects={(x['seed'],x['cell'],x['matching_reference'],x['physics'],x['dose']):x for x in r['scientific_comparisons']}
tables=[]
for dose,label in [('fractional','Fractional depletion'),('amount','Equal removed amount')]:
    rows=['| History / cell | Sham match, sham physics | Sham match, redistributed physics | Redistributed match, sham physics | Redistributed match, redistributed physics |',
          '| --- | ---: | ---: | ---: | ---: |']
    for s,c in cells:
        rows.append('| '+f'{s} / {c}'+' | '+' | '.join(f"{effects[s,c,L,P,dose]['full_effects']['fine']:.6f}" for L,P in [('S','S'),('S','R'),('R','S'),('R','R')])+' |')
    tables.append('### '+label+'\n\n'+'\n'.join(rows))
text='''# Delayed response after recipient matching — accepted assessment

**The large immediate-reset effects in cells 20/30 fall below the response screen after waiting, while cells 19/29 retain detected differences.** After the fixed 60-unit wait, history counts are **0/2, 1/2 and 1/2** in histories 7, 8 and 9 under both doses, both original matching references and both shared physical contexts. All numerical gates pass. Persistent chemical levels and persistent response behavior are distinct attributes in this model.

## What this experiment changes

The [prespecified protocol](phase_carry_recipient_delayed_response.md) challenges at t=570, after the matched targets and their entire chemical networks evolve freely for 60 units from t=510. Every new method starts from the same accepted fine unpulsed endpoint. There is **no second reset or continuing clamp**. Retain both original target-matching references, all selected recipients and both shared measured operator/volume contexts. Controls continue unperturbed while paired target pulses are observed through t=630.

All 45 source backgrounds agree with their parent reference and other-timestep endpoints within maximum raw log discrepancy **2.93e-13**, against the prespecified `1e-8` source gate. All pass the predefined descriptive settling screen. Maximum preceding 24-unit log drift is **6.71e-11**, and maximum absolute challenge-time log-concentration rate is **5.73e-14** per model time unit. This is measured state recovery, not an assumed equilibration or a proof of exact mathematical equilibrium.

Targets have regained different concentrations by challenge time, so the experiment includes both **10% fractional depletion** and **equal removed amount**. Common amounts are calculated separately for each original-reference/shared-physics pair from its actual delayed target levels. Each pulse is normalized by its own actual initial log amplitude and compared with its own continuing control. Equal removed amount controls the challenge dose; it does not equalize pre-existing target concentrations or network species totals.

There are **45 controls and 96 logical pulses**, giving **141 logical starts and 117 distinct initial-value problems**. Twenty-four dose preparations coincide; both labels are retained, without adding independent evidence. All new starts are calculated anew: 282 primary plus 141 diagnostic GPU paths, **423 GPU paths**, **423 CPU checks**, and **282 DOP853/Radau reference solves**. The independent units remain **three existing mature developmental histories**, with two selected recipients nested in each and zero new histories. The earlier wait trajectories and immediate challenges are accepted input/comparison evidence, not new replication.

## Delayed response results

The primary effect is the two-species time RMS of the difference between signed pulse/own-control log responses, normalized by each pulse's absolute initial log amplitude, over the full 60-unit observation window. The strict operational screen remains **greater than 0.01**. Each comparison holds its physical operator and volume vector common between chemical origins. All 48 decisions agree across the three timesteps and both reference solvers. Below-screen is not equivalence or absence of an effect.

'''+ '\n\n'.join(tables)+'''

Cells 19/29 clear the screen under all eight reference/dose/physical-context combinations: delayed effects range **0.0727–0.0770**. Cells 20/30 fall below it in every combination, with effects **0.00649–0.00870**. Their immediate fractional effects were approximately **0.0474–0.1006** after local matching. This supports a distinction between their large reset-dependent responses and their weaker responses after the chemical states have recovered. It does not imply that their delayed responses are zero.

Cell 22 also stays below screen, but its largest delayed effect is **0.009869**, close to the operational cutoff. Its weak detected classification should not be interpreted as a universal biological absence. Cell 23 has much smaller effects, approximately `5.12e-5–1.06e-4`.

All **96 logical fine pulse paths recover** toward their own controls, with target recovery **2.25–4.35 units** and volume-weighted network recovery **2.25–5.40 units**. Repeated dose labels do not make these independent recovery experiments. Recovery uses the same sustained-displacement criterion and remaining observation requirement as the parent study.

The largest descriptive physical-context contrast RMS is **0.004153**, and the largest signed interaction RMS is **0.003638**. These contrasts have no scientific positive-result screen. Their signed waveform identities hold to `5.55e-17`; scalar RMS values are not additive causal fractions.

## Concentration differences persist after the reset response subsides

The target states that reappeared during the waiting interval persist through the delayed observation. Cells 20/30 retain approximately **16.6–20.9-fold lower activator** under redistributed surrounding chemistry than under sham surrounding chemistry on either shared physical context. Their large chemical-level differences therefore coexist with below-screen delayed normalized response contrasts. A strong difference in concentration does not, by itself, establish a different response phenotype under this particular pulse assay.

All saved fine unpulsed controls remain unchanged at float64 observation precision during the new window. The source rate and dual-reference checks are the evidence for the near-stationary state; saved equality is not a general claim of exact equilibrium or unique attractors.

At fixed chemical origin, physical context and dose, the two original matching references give maximum normalized delayed response difference **3.60e-11**, and maximum full-window response RMS difference **2.55e-11**. Their agreement is descriptive, without a new equivalence screen. This contrasts with their different immediate response amplitudes: the initial local reset choice affects the transient response, but these two measured choices return to the same contextual state and delayed behavior in this window. It does not establish forgetting of every perturbation or absence of other basins.

At the delayed pulse, target chemistry is no longer equal between chemical origins. These delayed differences consequently support **conditional chemical-state response sensitivity**; they do not isolate an ongoing surrounding-network effect at fixed target chemistry. The accepted immediate-matching assay supplies the separate evidence that surrounding chemistry affects initially identical targets. Taken together, the assays support network-maintained local states and selective later behavior, without proving autonomous cell identity.

## Numerical qualification and evidence

| Saved-array check | Maximum / range | Prespecified limit |
| --- | ---: | ---: |
| DOP853/Radau raw log discrepancy | `4.47e-12` | `1e-8` |
| DOP853/Radau normalized pulse discrepancy | `4.24e-11` | `1e-8` |
| DOP853/Radau signed contrast discrepancy | `3.53e-11` | `1e-8` |
| Production/reference raw log discrepancy | `1.20e-6` | `0.01` |
| CPU/PyTorch raw log discrepancy | `2.45e-13` | `1e-8` |
| Production normalized pulse discrepancy | `1.15e-5` | `0.01` |
| Production signed contrast discrepancy | `4.47e-6` | `0.01` |
| RMS contrast discrepancy | `2.35e-7` | `0.01` |
| Relative activator-response area discrepancy | `7.71e-6` | `0.02` |
| Recovery-time discrepancy | `0` | `0.3` |
| Observed order, 96 resolved cases | `2.00044–2.00268` | `1.7–2.3` |
| Relative measured/requested pulse amount error | `1.71e-14` | `2e-13` |

All 141 dual-reference checks, 423 production accuracy checks, 282 raw timestep comparisons, 864 normalized pulse comparison rows and 1,080 signed contrast comparison rows qualify. All 96 pulse cases resolve second order. The 45 near-stationary controls have errors below the prespecified reference-noise resolution criterion; their order is reported as unresolved, not assigned an order of two. Saved positivity, padding and checkpoint consistency pass; unsaved internal Euler-stage checks remain unchanged runtime/source evidence.

A separate [saved-array reconstruction](../archive/study-reviews/2026-10-08-recipient-delayed-response/independent_review.json), importing no project analysis helpers, checks **12,084 recorded quantities** to maximum discrepancy `2.22e-16`. All **64 source pins, 369 input pins and 1,573 consumed artifacts** verify unchanged. All 120 duplicated dose/method comparisons have zero log discrepancy. The [launch evidence](phase_carry_recipient_delayed_response_launch.md) preserves eight tests, the literal input audit, all actual CUDA prefixes and the retained pre-preparation duplicate-keyword failure. No criterion changed.

The accepted grouped layout uses nine physical blocks and 23 independent chemical columns, with 66 positive-control padding columns excluded from reported counts. It uses unchanged float64 PyTorch matrix operations and the unchanged scientific SSP-RK2 kernel. No mechanics/polarity or geometry is advanced. GPU production took **168.58 seconds**; full elapsed coordinator time was **432.73 seconds**, about **7.2 minutes**, with four single-thread CPU verification workers overlapping the GPU. Worker-time sums are not elapsed time.

The [figure](../archive/study-reviews/2026-10-08-recipient-delayed-response/recipient-delayed-response.png), [SVG](../archive/study-reviews/2026-10-08-recipient-delayed-response/recipient-delayed-response.svg), [48-comparison CSV](../archive/study-reviews/2026-10-08-recipient-delayed-response/recipient-summary.csv), [compact assessment](phase_carry_recipient_delayed_response_assessment.json) and review provenance retain the evidence. Raw chemical paths remain under `output/studies/phase-carry-recipient-delayed-response`; no 3D fields were copied.

## Next priority and limits

**Repeat initial local matching with the surrounding network's initial activator and inhibitor amounts also controlled.** This would separate initial surrounding allocation from initial total dosage while preserving the distinction between immediate and delayed responses. Specify the adjustment and positivity/amount checks before launch; later reaction-driven amount differences remain outcomes rather than an assumed continuing clamp. This follow-up is not launched here.

This test does not match global network amounts, isolate direct neighbors, keep target chemistry equal at the delayed pulse, or advance mechanics. It uses three selected mature histories, historical no-carry developmental inheritance, sixteen-cell compartments and the inherited conductance approximation. It establishes neither autonomous/inherited biological identity, fresh carry-corrected zygote differentiation, population prevalence nor a universal polarity mechanism. Accepted parents, the results ledger, manuscript, baseline simulator and shared kernels remain unchanged.
'''
(REPO/'docs/phase_carry_recipient_delayed_response_assessment.md').write_text(text)

path=REPO/'README.md';lines=path.read_text().splitlines()
for i,line in enumerate(lines):
    if line.startswith('**Completed recipient-state matching:**'):
        lines[i]=line.replace('**Next: delay the pulse after matching** to distinguish reset relaxation from settled behavior; the follow-up is not launched.',
            'The completed delayed assay below distinguishes reset relaxation from behavior after state recovery.')
    if line.startswith('**Delayed recipient challenge is running:**'):
        lines[i]=('**Completed delayed recipient challenge:** the [assessment](docs/phase_carry_recipient_delayed_response_assessment.md) '
            'accepts 423 GPU paths, 423 CPU checks and 282 reference solves across 141 logical starts (117 distinct problems). '
            'After a fixed 60-unit wait, cells 19/29 retain above-screen responses, while cells 20/30 fall below screen: counts '
            'are 0/2, 1/2 and 1/2 across three existing histories under both doses, both original matching references and both '
            'physical contexts. All tested pulses recover. Cells 20/30 still retain about 16.6–20.9-fold activator differences, '
            'separating persistent chemical levels from large immediate reset responses. The two original matching references '
            'give nearly identical delayed behavior (maximum normalized difference 3.60e-11). All 45 source backgrounds pass '
            'the predefined settling screens. The target is not rematched at the delayed pulse, and surrounding totals remain '
            'unmatched; this supports conditional state sensitivity, not autonomous identity. All numerical gates and immutable '
            'evidence checks pass. **Next: control initial surrounding species amounts alongside local matching.** The follow-up '
            'is not launched. [Figure](archive/study-reviews/2026-10-08-recipient-delayed-response/recipient-delayed-response.png).')
path.write_text('\n'.join(lines)+'\n')
path=REPO/'docs/plan.md';lines=path.read_text().splitlines()
for i,line in enumerate(lines):
    if line.startswith('**Current priority: qualify the delayed recipient challenge.**'):
        lines[i]=('**Completed delayed recipient challenge:** the [assessment](phase_carry_recipient_delayed_response_assessment.md) '
            'qualifies 423 GPU paths, 423 CPU checks and 282 references, with all 48 decisions agreeing across methods. After '
            'the fixed 60-unit wait, cells 19/29 remain above screen and cells 20/30 fall below it under both doses/references/physics: '
            'history counts 0/2, 1/2 and 1/2. All pulses recover. Cells 20/30 retain large concentration differences, and the two '
            'initial matching references return to nearly identical delayed responses. All source/settling checks and separate '
            'saved-array reconstructions pass. These are three existing mature histories, not new developmental replication.')
        lines.insert(i+1,'')
        lines.insert(i+2,('**Current priority: control initial surrounding chemical amounts.** Define a separate intervention that '
            'matches the recipient state and the surrounding network’s initial species totals on shared physics, while retaining '
            'distinct allocation. Specify positivity, per-cell/total amount accounting and immediate/delayed outcomes before '
            'launch. Later reaction-driven amount changes are measured outcomes. This is not launched; no continuing clamp, '
            'ledger/manuscript promotion or autonomous identity claim is implied.'))
        break
s='\n'.join(lines)+'\n'
s=s.replace('Those immediate reset responses still require delayed challenges before being interpreted as settled behavior.',
    'The [accepted delayed recipient assay](phase_carry_recipient_delayed_response_assessment.md) returns to selective '
    '0/2, 1/2 and 1/2 outcomes after a verified settling interval, with large concentration differences retained in cells 20/30 '
    'but below-screen delayed responses. Target chemistry differs again at the delayed pulse, and surrounding totals remain unmatched.')
path.write_text(s)
print('Saved accepted assessment and updated README/plan reporting.')
