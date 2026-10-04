"""Annotate historical claims without rerunning or changing scientific decisions.

Outcome agreement is a candidate for revalidation, not evidence of carry
insensitivity. Frozen chemistry inherits its geometric/chemical preparations.
"""
import argparse
from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
ORIGINAL = 'archive/study-document-snapshots/2026-10-04-before-claim-precision-audit/docs/results_ledger.json'


def record(key, dependency, outcome, trajectory, action, *, diagnostic='', supplemental=False):
    return dict(id=key, dependency=dependency, outcome_claims=outcome,
                trajectory_claims=trajectory, diagnostic_claim=diagnostic,
                required_revalidation=action, supplemental=supplemental,
                carry_status=('not_applicable_to_reference_test' if dependency == 'reference'
                              else 'not_revalidated_with_carry'))


RECORDS = [
    record('initiation', 'moving_no_carry',
           'Formation versus suppression in the paired zygote history.',
           'Developmental contrast curves, onset, and late contrast magnitudes.',
           'Repeat paired development with carry from the zygote; mature restarts do not replace this test.'),
    record('polarity-components', 'moving_no_carry',
           'Which mechanical components permit or suppress formation.',
           'Component response amplitudes, formation times, and geometric trajectories.',
           'Repeat matched component controls with carry, preserving chemical polarity modulation and initial state.'),
    record('homogeneous-spectrum', 'moving_no_carry', '',
           'Developmental spectral-band exits, modal growth, and counterfactual geometric attribution.',
           'Recompute spectra on carry-generated snapshots; unchanged kinetics alone does not fix the moving graph.'),
    record('moving-maintenance', 'moving_no_carry',
           'Developed contrast survives the switch to full coupling in histories 7–9.',
           'Contrast minima, concentration paths, cell-association correlations, and response timing.',
           'Repeat paired t=90–150 continuations with carry in all three histories; separately audit upstream development.'),
    record('frozen-coexistence', 'frozen_inherited_no_carry',
           'Uniform and patterned states coexist locally on each tested endpoint graph.',
           'Equilibrium concentrations, eigenvalues, and recovery errors/times.',
           'Retain validity conditional on the old graphs; rerun the same basin/stability assay on actual carry endpoints.'),
    record('original-frozen-recovery', 'frozen_inherited_no_carry',
           'Perturbed states recover; strong contrast is retained on the original no-feedback graph.',
           'Endpoint contrast, relaxation paths, and recovery times.',
           'Reassay carry-generated endpoints with matched chemical starts and unchanged perturbation criteria.'),
    record('all-pair-exchange', 'frozen_inherited_no_carry',
           'Return/reorganization/likeness classifications for chemical transplants.',
           'Distances to transferred states, response amplitudes, and the 48/48 split of 96 informative cases.',
           'Reprepare graphs and chemical states with carry; repeat pair selection, report classification margins and changed counts.'),
    record('common-environment', 'frozen_inherited_no_carry',
           'Tested isolation releases erase differences; reservoir conditions admit multiple stable states.',
           'Concentrations, response sizes, release times, and regime boundaries.',
           'Repeat sampled-state preparations with carry; distinguish intrinsic reservoir-ODE behavior from upstream state sampling.'),
    record('finite-reservoir', 'frozen_inherited_no_carry',
           'Some finite-reservoir regimes sustain differences; others lose them.',
           'Reservoir/chemical trajectories, transition times, and quantitative regime boundaries.',
           'Reassay carry-derived initial states with identical amounts, reservoir size and exchange parameters.'),
    record('frozen-exchange-response', 'frozen_inherited_no_carry',
           'Exchanged-cell response is nearer the donor reference on the tested graphs.',
           'Waveform distances, donor-nearer margins, pulse amplitudes, and settling/recovery times.',
           'Reprepare carry-derived endpoints and states; rerun donor/recipient/sham responses and check sign margins.'),
    record('moving-exchange-response', 'moving_no_carry',
           'Donor-nearer response sign and reorganization after exchange in histories 7–9.',
           'Response amplitudes, waveform distances, recovery times, and cell concentration paths.',
           'Repeat carry formation/preparation and paired moving exchange, sham and pulse controls; retain original dose definitions.'),
    record('neighbor-context', 'frozen_inherited_no_carry',
           'Conservative neighbor averaging shifts the recipient; state/response classifications change.',
           'Mixing ratios, late-state shifts, waveform distances, and 11/36 and 46/72 intervention counts.',
           'Reprepare carry endpoints and matched challenges; recompute classifications, counts and threshold margins by history.'),
    record('delayed-pulses', 'frozen_inherited_no_carry',
           'The two selected immediate nonrecoveries recover after settling.',
           'Delay-dependent recovery times, pulse amplitudes and amount-matched dose factors.',
           'Reassess selection on carry-derived challenges before repeating fractional and matched-amount pulses.'),
    record('gpu-validation', 'old_method_numerical', '', '',
           'Keep original backend comparisons; require context gates for carry and do not infer long-horizon equivalence from a short prefix.',
           diagnostic='CPU/GPU agreement qualifies the tested no-carry implementation, not carry insensitivity.'),
    record('mature-timestep', 'old_method_numerical', '', '',
           'Repeat the selected mature response/retention refinement with carry after matched preparations.',
           diagnostic='Old-method timestep passes do not bound carry/no-carry method-change errors.'),
    record('survival-timestep', 'old_method_numerical', '', '',
           'Repeat carry maintenance refinement and separate continuation accuracy from upstream developmental accuracy.',
           diagnostic='The recorded t=90–150 pass applies only to the old accumulation method.'),
    record('conservative-reference', 'reference', '', '',
           'No phase-carry rerun is required for this controlled operator convergence test.',
           diagnostic='Manufactured regular-grid transport convergence; no moving phase accumulation.'),
    record('geometric-closure', 'reference', '', '',
           'Retain the declared flat passes and general-geometry failures; carry does not repair conductance closure.',
           diagnostic='Controlled geometric closure tests, not a carry-dependent moving trajectory.'),
    record('positive-skew-prototype', 'reference', '', '',
           'Retain prototype scope; a reference pass does not promote the prototype into live simulations.',
           diagnostic='Separate prescribed conservative-positive flux reference patches.'),
    record('boundary-skew-prototype', 'reference', '', '',
           'Retain prototype scope; no live-backend substitution.',
           diagnostic='Separate prescribed boundary-corrected flux reference patches.'),
    record('historical-development-refinement', 'old_method_numerical', '', '',
           'Retain the original failure. A new carry developmental gate would be a separate result.',
           diagnostic='The historical full-development gate failed; later carry continuation passes do not accept it.'),
    record('cleavage-measurement', 'calibrated_diagnostic', '', '',
           'Retain analytic measurement calibration; full moving cleavage/development still needs its own carry validation.',
           diagnostic='Diagnostic calibration resolves a measurement problem, not phase accumulation or developmental convergence.'),
    record('parameter-robustness', 'moving_and_frozen_inherited_no_carry',
           'Formation/maintenance and endpoint coexistence classifications at the three tested parameter points.',
           'Spectral crossing times, conductance/area decline percentages, correlations, shape changes, and history-9 frozen/moving contrast ratio 262.1.',
           'Repeat matched moving/frozen controls with carry, then reassay actual endpoints; rederive all quantitative comparisons.', supplemental=True),
    record('polarity-robustness', 'moving_and_frozen_inherited_no_carry',
           'History-9 formation at zero directional tension; developed-pattern persistence across sampled contrasts.',
           'Onsets, spectral crossings, conductance declines, concentration paths and contrast magnitudes.',
           'Running carry controls address near-uniform initiation at chi=0/0.35 in histories 7–9; chi=0.7 and developed-start maintenance remain untested.', supplemental=True),
]

AUDIT_EVIDENCE = [
    'outputs/phase-carry-formation/summary.json',
    'outputs/phase-carry-convergence/summary.json',
    'outputs/parameter-robustness-moving/assessment.json',
    'outputs/polarity-robustness/assessment.json',
]
CARRY_SCOPE = ('One existing history (9), mature t=150–390, near-uniform start, '
               'chi=0, D_b/D_a=27.5; three timestep paths, persistent formation '
               'and actual frozen-endpoint coexistence. Initial development still used no carry.')
DISPOSITION = ('No curated biological claim is broadly revalidated with carry. '
               'Outcome-level claims are candidates, not certified robust. '
               'Trajectory-level quantities remain original no-carry results until rederived.')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate_coverage(ledger):
    ids = [r['id'] for r in RECORDS if not r['supplemental']]
    actual = [r['id'] for r in ledger['entries']]
    if len(ids) != len(set(ids)) or len(actual) != len(set(actual)) or set(ids) != set(actual):
        raise ValueError('Claim audit must explicitly cover every curated ledger entry.')


def annotate(ledger):
    validate_coverage(ledger)
    result = deepcopy(ledger)
    result.pop('precision_audit', None)
    # Annotate the original scientific snapshot rather than silently replacing it.
    result['captured_utc'] = json.loads((REPO/ORIGINAL).read_text())['captured_utc']
    indexed = {r['id']: r for r in RECORDS}
    for entry in result['entries']:
        entry['precision_audit'] = deepcopy(indexed[entry['id']])
    formation = json.loads((REPO/AUDIT_EVIDENCE[0]).read_text())
    qualification = json.loads((REPO/AUDIT_EVIDENCE[1]).read_text())
    if (not formation['passed'] or not qualification['passed']
            or formation['independent_histories'] != 1 or qualification['independent_histories'] != 1):
        raise ValueError('The limited completed carry support must match its actual recorded decision.')
    result['precision_audit'] = dict(
        schema=1, assessed_utc=datetime.now(timezone.utc).isoformat(),
        original_snapshot=ORIGINAL, original_snapshot_sha256=digest(REPO/ORIGINAL),
        original_generator_sha256=json.loads((REPO/ORIGINAL).read_text())['source_sha256'],
        records=deepcopy(RECORDS),
        evidence_sha256={name: digest(REPO/name) for name in AUDIT_EVIDENCE},
        same_dt_chemical_log_max={k:v['errors']['chemical_log_max']
                                 for k,v in formation['same_dt_interventions'].items()},
        qualified_carry_scope=CARRY_SCOPE, disposition=DISPOSITION,
        source_sha256={'embryo/precision_claim_audit.py':digest(__file__)})
    result['source_sha256']['embryo/results_ledger.py'] = digest(REPO/'embryo/results_ledger.py')
    verify_preserved(result)
    return result


def verify_preserved(ledger):
    original = json.loads((REPO/ORIGINAL).read_text())
    stripped = [{k:v for k,v in e.items() if k != 'precision_audit'} for e in ledger['entries']]
    if stripped != original['entries']:
        raise ValueError('Audit must not rewrite original claims, decisions, counts or limitations.')
    for key in ('captured_utc', 'evidence_sha256', 'history_checkpoint_sha256',
                'history_registry', 'unique_primary_histories'):
        if ledger[key] != original[key]:
            raise ValueError('Original evidence snapshot changed: '+key)


def verify(ledger):
    validate_coverage(ledger)
    verify_preserved(ledger)
    audit = ledger['precision_audit']
    if audit['qualified_carry_scope'] != CARRY_SCOPE or audit['disposition'] != DISPOSITION:
        raise ValueError('Limited carry support cannot be promoted to blanket validation.')
    if audit['records'] != RECORDS:
        raise ValueError('Audit registry and rendered annotations differ.')
    indexed = {r['id']:r for r in RECORDS}
    for entry in ledger['entries']:
        if entry['precision_audit'] != indexed[entry['id']]:
            raise ValueError('Entry annotation differs: '+entry['id'])
    checks = {audit['original_snapshot']:audit['original_snapshot_sha256'],
              **audit['source_sha256'], **audit['evidence_sha256'],
              **ledger['source_sha256'], **ledger['evidence_sha256'],
              **ledger['history_checkpoint_sha256']}
    for name, expected in checks.items():
        if digest(REPO/name) != expected:
            raise ValueError('Precision audit evidence changed: '+name)
    formation = json.loads((REPO/AUDIT_EVIDENCE[0]).read_text())
    observed = {k:v['errors']['chemical_log_max'] for k,v in formation['same_dt_interventions'].items()}
    if audit['same_dt_chemical_log_max'] != observed:
        raise ValueError('Carry/no-carry comparison changed.')


def markdown_overlay(ledger, original_markdown):
    audit = ledger['precision_audit']
    differences = audit['same_dt_chemical_log_max']
    lines = [
        '## Phase-update precision audit', '',
        f"Claim audit: {audit['assessed_utc']}. The evidence snapshot below retains its original date and decisions. **All curated biological results use no-carry moving paths or inherit their graphs/states. They are not broadly carry-validated.**", '',
        f"At fixed timestep, carry versus no carry changes maximum log concentration by **{differences['fine']:.8f}** (dt=0.001875) and **{differences['finer']:.8f}** (dt=0.0009375). Old-method timestep passes do not bound this change. Matching final contrast in one case does not establish matching trajectories or robustness of other claims.", '',
        '**Outcome-level claims:** retained as original-method observations and candidates for revalidation. None is certified insensitive to carry merely because it is categorical. Signs and threshold classifications require margins as well as matched reruns.', '',
        '| Study | Outcome-level claim | Carry disposition |',
        '|---|---|---|',
    ]
    for row in audit['records']:
        if row['outcome_claims']:
            suffix=' (later study; outside original 22)' if row['supplemental'] else ''
            lines.append(f"| `{row['id']}`{suffix} | {row['outcome_claims']} | Revalidation pending. |")
    lines += ['', '**Trajectory-level claims:** the quantities below are computed with the original method, including those in reports whose own numerical gates passed. Recompute on carry trajectories before using them as current-method quantitative evidence.', '',
              '| Study | Original-method quantities requiring rederivation |', '|---|---|']
    for row in audit['records']:
        if row['trajectory_claims']:
            suffix=' (later study)' if row['supplemental'] else ''
            lines.append(f"| `{row['id']}`{suffix} | {row['trajectory_claims']} |")
    lines += ['', 'Counts derived from classifications (for example 48/48, 11/36 or 46/72) also need recalculation; they are not protected by their categorical source. The history-9 **262.1-fold** frozen/moving contrast ratio remains an **old-method** result. Its small moving denominator makes quantitative revalidation particularly important.', '',
              '### Frozen assays and upstream dependence', '',
              'Frozen chemical ODE integration has no moving phase update. Independent-solver agreement and local stability remain evidence **conditional on the original graphs and chemical preparations**. To generalize these findings to carry mechanics, generate carry endpoints/states and repeat the original assays there. An additional frozen-ODE tolerance check on the same old inputs cannot test this dependence.', '',
              'A carry continuation from a stored t=90 or t=150 state still inherits no-carry development. It tests the continuation conditional on that state. Fresh carry development and cleavage qualification are required for zygote-to-organization claims.', '',
              '### Numerical and reference claims', '',
              '| Study | What the recorded decision establishes | Carry action |', '|---|---|---|']
    for row in audit['records']:
        if row['diagnostic_claim']:
            lines.append(f"| `{row['id']}` | {row['diagnostic_claim']} | {row['required_revalidation']} |")
    lines += ['', '### Completed carry evidence and revalidation order', '',
              audit['qualified_carry_scope'], '',
              'The [completed formation comparison](phase_carry_formation.md) and [three-timestep qualification](phase_carry_convergence.md) support that limited case. Neither replaces the original ledger rows, adds a history, or validates the 262.1-fold comparison at its different contrast/window. Three existing carry endpoints supporting coexistence do not revalidate all six older graphs or the reservoir/exchange studies.', '',
              '1. Finish and assess the already-running [matched carry polarity controls](phase_carry_polarity.md): three histories, near-uniform starts, chi=0/0.35, two timesteps. Keep numerical acceptance separate from initiation/suppression outcomes.',
              '2. Repeat three-history carry maintenance and actual-endpoint coexistence; then matched moving exchange/sham/pulse assays, reporting donor-nearer margins and response amplitudes.',
              '3. Rebuild frozen neighbor, exchange and reservoir preparations from carry endpoints; retain perturbations and report outcomes, margins, counts and trajectories separately.',
              '4. Recompute parameter/polarity quantitative mechanisms with carry: matched moving/frozen ratios, spectra, crossing/onset times and conductance changes. Remaining contrasts, developed starts and fresh developmental histories need separate gates.', '',
              'No new science is launched by this audit. Do not expand the parameter map or treat original-method passes as carry-insensitivity evidence. The [machine-readable audit](results_ledger.json) records the dependency and required rerun for every curated entry plus the two later robustness studies.', '']
    body = original_markdown.replace('| Validation |', '| Validation in original record |')
    body = body.replace('## Reproduce the consolidation', '## Reproduce and verify the audited consolidation')
    body = body.replace('No new scientific experiment was launched during consolidation. The next scientific decision remains moving-geometry confirmation of the context interventions, after reviewing this ledger and manuscript.',
        'No scientific output is changed by the audit. Before expansion, revalidate the claims listed above with carry. Run `python -m embryo.precision_claim_audit --verify` for the additional claim/dependency and original-record preservation checks. A ledger refresh that changes curated claims requires a new audit baseline; it cannot silently reuse these labels.')
    return body.replace('## Initiation\n', '\n'.join(lines)+'\n## Initiation\n', 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    target = REPO/'docs/results_ledger.json'
    ledger = json.loads(target.read_text())
    from .results_ledger import markdown
    if args.verify:
        verify(ledger)
        if (REPO/'docs/results_ledger.md').read_text() != markdown(ledger):
            raise ValueError('Audited Markdown does not match the machine-readable ledger.')
        print('All 22 entries and two supplemental studies audited; original decisions and evidence preserved.')
    else:
        ledger = annotate(ledger)
        target.write_text(json.dumps(ledger, indent=2, allow_nan=False)+'\n')
        (REPO/'docs/results_ledger.md').write_text(markdown(ledger))
        print('Annotated 22 original entries and two later robustness studies; no scientific output changed.')


if __name__ == '__main__':
    main()
