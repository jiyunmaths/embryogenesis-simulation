"""Snapshot the consolidated evidence without running or changing simulations.

Histories are explicitly registered and reused across assays. Directory status,
numerical acceptance, and scientific interpretation are separate fields.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from .feedback_long import digest

REPO = Path(__file__).resolve().parent.parent
HISTORIES = {str(seed): dict(id=f'attribute-history-{seed}', seed=seed,
    model='Conservative AttributeSimulation; no supplied fate switch',
    lineage='No-feedback zygote development to t=90; mature descendants reused in later assays',
    source=('outputs/attribute-development/no_feedback/final_state.npz' if seed == 7 else
            f'outputs/feedback-survival-validation/seed-{seed}/development/latest_state.npz'))
    for seed in (7, 8, 9)}


def snapshot():
    reports, entries = {}, []

    def read(name):
        relative = 'outputs/'+name
        if relative not in reports:
            reports[relative] = json.loads((REPO/relative).read_text())
        return reports[relative]

    def add(key, theme, seeds, claim, validation, limits, files, details=None):
        for name in files: read(name)
        entries.append(dict(id=key, theme=theme,
            history_ids=[HISTORIES[str(seed)]['id'] for seed in seeds],
            histories=len(seeds),
            claim=claim, validation=validation, limitations=limits,
            evidence=['outputs/'+f for f in files], details=details or {}))

    development = read('attribute-development/comparison.json')
    add('initiation', 'Initiation', [7],
        'No-feedback development produces strong chemical contrast; full coupling suppresses it in the paired zygote comparison.',
        'quality_pass' if development['all_quality_pass'] else 'quality_fail',
        'One paired history; formation contrast is not a discrete identity or shape-causality result.',
        ['attribute-development/protocol.json', 'attribute-development/comparison.json'],
        {k:dict(log_activator_spread=v['late_mean_attribute_spread'][0],
                mean_axis_ratio=v['late_mean_axis_ratio']) for k,v in development['arms'].items()})
    components = read('feedback-long/comparison.json')
    ablation = read('feedback-polarity-ablation/comparison.json')
    component_pass = components['all_quality_pass'] and all(v['quality_pass'] for v in ablation['new_arms'].values())
    add('polarity-components', 'Initiation', [7],
        'Polar mechanics alone suppresses initiation; removing its mechanical contribution from tension+adhesion restores formation on the matched geometry.',
        'quality_pass' if component_pass else 'quality_fail',
        'One mature geometry, two prepared chemical families, twelve continuations; chemical modulation of polarity remains active.',
        ['feedback-long/protocol.json', 'feedback-long/comparison.json',
         'feedback-polarity-ablation/protocol.json', 'feedback-polarity-ablation/comparison.json'],
        dict(component_continuations=len(components['results'])+len(ablation['new_arms']),
             ablation_verdict=ablation['verdict']))
    add('homogeneous-spectrum', 'Initiation', [7],
        'The direct branch loses graph-supported homogeneous instability during development.',
        'descriptive_spectral_analysis',
        'Instantaneous frozen-graph linearization; not a nonlinear-maintenance or moving-system stability proof.',
        ['feedback-spectrum/protocol.json', 'feedback-spectrum/spectra.json', 'feedback-spectrum/counterfactuals.json'])
    maintenance_files = ['feedback-survival/protocol.json',
                         'feedback-survival-validation/summary.json']
    maintenance = {}
    for seed in (7,8,9):
        name = ('feedback-survival' if seed == 7 else f'feedback-survival-validation/seed-{seed}/survival')
        data = read(name+'/comparison.json'); maintenance[str(seed)] = data['arms']
        maintenance_files += [name+'/comparison.json']
    maintenance_pass = all(v['quality_pass'] for arms in maintenance.values() for v in arms.values())
    add('moving-maintenance', 'Maintenance', [7,8,9],
        'All three histories retain developed chemical contrast and cell association after switching full coupling on; paired keep-off controls also retain it.',
        'quality_pass' if maintenance_pass else 'quality_fail',
        'Finite horizon t=90-150 after the 16-cell cap; not direct-feedback formation, continued division, or inherited identity.',
        maintenance_files, dict(by_history=maintenance, paired_branches=6))
    endpoints, endpoint_files = {}, []
    for seed in (7,8,9):
        root = 'feedback-endpoint-bistability'+('' if seed == 7 else f'-seed-{seed}')
        data = read(root+'/results.json'); endpoints[str(seed)] = {
            k:dict(bistability_supported=v['chemical_bistability_supported'], trials=len(v['trials']))
            for k,v in data.items()}
        endpoint_files += [root+'/protocol.json', root+'/results.json']
    endpoint_pass = all(v['bistability_supported'] for d in endpoints.values() for v in d.values())
    add('frozen-coexistence', 'Maintenance', [7,8,9],
        'Both endpoint graphs in each of the three histories support locally stable uniform and patterned chemical states.',
        'chemical_acceptance_pass' if endpoint_pass else 'chemical_acceptance_fail',
        'Six graphs and 252 nested starts are not independent histories; no full moving-system stability or bifurcation classification.',
        endpoint_files, dict(by_history=endpoints, graphs=6, chemical_starts=sum(v['trials'] for d in endpoints.values() for v in d.values())))
    persistence = read('attribute-persistence/results.json')
    add('original-frozen-recovery', 'Maintenance', [7],
        'Small perturbations recover on both original graphs; only the no-feedback graph retains strong chemical contrast.',
        'tolerance_and_recovery_pass' if all(v['reference_pass'] and all(t['recovered'] for t in v['trials']) for v in persistence.values()) else 'check_fail',
        'Two graphs from one paired history; recovery of uniform chemistry is not differentiation.',
        ['attribute-persistence/protocol.json', 'attribute-persistence/results.json'])
    exchange = read('attribute-exchange/results.json')
    add('all-pair-exchange', 'Context', [7],
        'Among 96 informative chemical transplants, half return and half reorganize; none meets transferred-concentration likeness.',
        'tolerance_pass' if exchange['reference_pass'] else 'tolerance_fail',
        '120 pairs within one fixed graph; post-run endpoint stability is exploratory, not independent replication.',
        ['attribute-exchange/protocol.json', 'attribute-exchange/results.json', 'attribute-exchange/assessment.json'],
        dict(pair_counts=exchange['all_pairs']))
    common = read('attribute-common-environment/verification.json')
    add('common-environment', 'Context', [7],
        'Isolation erases differences in the tested releases; common reservoir exchange can support multiple stable chemical states.',
        'independent_checks_pass' if all(v['passed'] for v in common.values()) else 'independent_checks_fail',
        '2256 sampled cell states descend from one history; fixed reservoirs supply external support. An initial strength-4 tolerance failure is retained; further refinement passes, with one selected independent Radau check per arm.',
        ['attribute-common-environment/protocol.json', 'attribute-common-environment/results.json',
         'attribute-common-environment/verification.json'])
    finite = read('attribute-finite-reservoir/summary.json')
    add('finite-reservoir', 'Context', [7],
        'An evolving finite reservoir can support differences; formation and release depend on exchange regime.',
        'qualified_with_retained_failures',
        'Original strong-release numerical failure, unstable symmetric endpoints, and an unsettled follow-up remain exceptions; not new developmental histories.',
        ['attribute-finite-reservoir/protocol.json', 'attribute-finite-reservoir/results.json',
         'attribute-finite-reservoir/verification.json', 'attribute-finite-reservoir/summary.json'], finite)
    frozen_response = read('cell-response-exchange/assessment.json')
    add('frozen-exchange-response', 'Context', [7,8,9],
        'Exchanged-cell responses favor the donor reference on both endpoint graphs in each history.',
        'independent_solver_and_settling_pass' if frozen_response['all_relaxations_settled'] and frozen_response['max_solver_log_error'] < 1e-5 else 'check_fail',
        '48 reference comparisons nested in three histories; donor-nearer is not unchanged transfer or autonomous identity.',
        ['cell-response-exchange/protocol.json', 'cell-response-exchange/results.json', 'cell-response-exchange/assessment.json'], frozen_response)
    moving = read('exchange-response-histories/comparison.json')
    add('moving-exchange-response', 'Context', [7,8,9],
        'Each history has eight donor-nearer moving response comparisons; exchanged chemistry reorganizes rather than retaining transferred values unchanged.',
        'completed_quality_and_backend_checks',
        'Prepared mature patterned basins; 24 comparisons reuse the same three histories, not new formation or population-frequency evidence.',
        ['cell-exchange-response-moving/protocol.json', 'cell-exchange-response-moving/comparison.json',
         'exchange-response-histories/protocol.json', 'exchange-response-histories/status.json', 'exchange-response-histories/comparison.json'],
        dict(by_history={'7':dict(donor=sum(x['nearest_reference']=='donor' for x in moving['historical_seed7_comparisons']))}
             | {str(h['seed']):h['descriptive_counts'] for h in moving['new_histories']}))
    neighbors = read('neighbor-context/results.json'); summary = read('neighbor-context/summary.json')
    by_history = {}
    for seed in (7,8,9):
        contexts = [c for c in neighbors['contexts'] if c['seed']==seed]
        rows = [r for c in contexts for r in c['comparisons']]
        low = {7:20,8:20,9:30}[seed]
        mix = next(r for c in contexts if c['background']=='fresh_exchange' for r in c['comparisons']
                   if r['target']==low and r['arm']=='neighbor_mix')
        by_history[str(seed)] = dict(challenges=len(rows), state_shifts=sum(r['state_outcome']=='shifted' for r in rows),
            pulse_comparisons=2*len(rows), response_shifts=sum(t['response_shift'] for r in rows for t in r['responses']),
            conservative_recipient_state_ratio=mix['late_state_ratio'])
    add('neighbor-context', 'Context', [7,8,9],
        'Conservative neighbor averaging shifts the fresh low-state recipient in all three histories; broader challenges alter late states and pulse responses.',
        'independent_solver_and_sham_pass' if summary['numerical_pass'] else 'check_fail',
        '36 challenges and 72 comparisons nested in three histories; immediate pulses occur during reorganization; geometry is frozen.',
        ['neighbor-context/protocol.json', 'neighbor-context/results.json', 'neighbor-context/summary.json', 'neighbor-context/status.json'],
        dict(by_history=by_history, nested_counts=summary['counts']))
    delayed = read('neighbor-context-delayed/results.json')
    add('delayed-pulses', 'Context', [8],
        'Both selected seed-8 immediate-pulse nonrecoveries recover when challenged after settling, under percentage and original-amount controls; alternative endpoints also resist these pulses.',
        'independent_solver_and_control_pass' if delayed['numerical_pass'] else 'check_fail',
        'Post-selected follow-up of two cases in one reused history; 24 pulses are not a replication rate, a measured delay window, or moving-context evidence.',
        ['neighbor-context-delayed/protocol.json', 'neighbor-context-delayed/results.json', 'neighbor-context-delayed/status.json',
         'neighbor-context-delayed/verification.json'],
        dict(counts=delayed['counts'], max_solver_log_error=delayed['max_solver_log_error'],
             max_control_drift_log=delayed['max_control_drift_log']))
    gpu = read('gpu-backend-validation/comparison.json')
    add('gpu-validation', 'Numerical', [7],
        'Four full accepted mature CPU/GPU control/pulse replays pass; additional contexts are checked before scientific use.',
        'backend_pass' if gpu['passed'] and gpu['scientific_ready'] else 'backend_fail',
        'Backend agreement in the tested mature regime, not an extra history, developmental convergence, or GPU division validation.',
        ['gpu-backend-validation/protocol.json', 'gpu-backend-validation/comparison.json', 'gpu-backend-validation/status.json'])
    refined = read('exchange-response-histories-refined/refinement.json')
    add('mature-timestep', 'Numerical', [8,9],
        'Both existing histories pass mature exchange/retention and selected same-state response timestep halving; classifications and sampled recoveries agree.',
        'refinement_pass' if refined['passed'] else 'refinement_fail',
        '18 repeated continuations, four full backend replays, ten context checks; no new history, spatial convergence, or refined development.',
        ['exchange-response-histories-refined/protocol.json', 'exchange-response-histories-refined/refinement.json',
         'exchange-response-histories-refined/status.json'])
    for key, seeds, files, claim, limits in (
        ('survival-timestep', [7], ['feedback-survival-validation/refinement.json'],
         'The paired t=90-150 seed-7 maintenance continuation passes timestep halving.', 'Same physical starts; no validation of preceding cleavage.'),
        ('conservative-reference', [], ['live-transport-validation-complete/report.json'],
         'Conservative volume-weighted transport converges on the controlled regular reference; random-walk scaling does not give fixed physical diffusivity.', 'Regular reference does not validate diffuse-contact conductance closure.'),
        ('geometric-closure', [], ['geometry-transport-validation/comparison.json', 'geometry-transport-validation/status.json'],
         'Flat orthogonal reference checks pass; all three declared general-geometry closure checks fail.', 'Curvature, diffuse gaps, and nonorthogonal flux remain physical approximation limits.'),
        ('positive-skew-prototype', [], ['positive-skew-flux/comparison.json'],
         'The separate conservative-positive flux prototype passes its specified patches.', 'Not coupled to live mechanics; does not repair the current conductance approximation.'),
        ('boundary-skew-prototype', [], ['skew-boundary-correction/comparison.json'],
         'The separate corrected-boundary prototype passes its prescribed patches.', 'Not live-integrated; no general moving-tissue acceptance.'),
        ('historical-development-refinement', [], ['development-refinement/comparison.json', 'development-refinement/status.json'],
         'The historical full-development refinement completes five cases but fails the combined predeclared gate.', 'Historical fate model; not three new attribute-model histories, and no full-development acceptance transferred.'),
        ('cleavage-measurement', [], ['cleavage-measurement/comparison.json'],
         'A separate calibrated diagnostic resolves the cleavage measurement failure without rewriting the original failed screen.', 'Controlled analytic cleavage, not full developmental or sharp-interface convergence.'),
    ):
        data = read(files[0]); decision = data.get('passed')
        if key == 'survival-timestep': decision = all(v['passed'] for v in data.values())
        add(key, 'Numerical', seeds, claim, 'pass' if decision else 'fail', limits, files)

    inventory = []
    curated_roots = {Path(name).parts[1] for name in reports}
    for root in sorted((REPO/'outputs').iterdir()):
        if not root.is_dir(): continue
        state_file = root/'status.json'
        status = json.loads(state_file.read_text()) if state_file.exists() else {}
        inventory.append(dict(directory=str(root.relative_to(REPO)),
            disposition=('cache' if 'cache' in root.name else 'curated' if root.name in curated_roots else 'historical_or_unassessed'),
            runner_state=status.get('state', 'no_root_status'),
            acceptance_flags={k:status[k] for k in ('passed','numerical_pass','all_checks_pass','all_quality_pass') if k in status},
            root_reports=[p.name for p in sorted(root.glob('*.json')) if p.name not in ('protocol.json','launch.json')]))
    return dict(schema=1, captured_utc=datetime.now(timezone.utc).isoformat(),
        replication_unit='Developmental history; cells, branches, graphs, pulses, preparation choices and numerical repeats are nested.',
        unique_primary_histories=3, histories_added_by_consolidation=0, history_registry=HISTORIES,
        history_count_note='Three unique attribute histories globally, reused across rows. The paired seed-7 direct branch is an intervention in that registered random-stream history. Historical fate-model and analytic tests are not pooled with these histories.',
        scope='Consolidation snapshot; no new scientific simulations. Curated decisions and evidence are separate from the complete directory inventory.',
        entries=entries, evidence_sha256={name:digest(REPO/name) for name in reports},
        source_sha256={str(Path(__file__).relative_to(REPO)):digest(__file__)},
        history_checkpoint_sha256={v['source']:digest(REPO/v['source']) for v in HISTORIES.values()},
        directory_inventory=inventory)


def markdown(ledger):
    if 'precision_audit' in ledger:
        from .precision_claim_audit import markdown_overlay
        original = dict(ledger)
        original.pop('precision_audit')
        return markdown_overlay(ledger, markdown(original))
    lines = ['# Results ledger', '', 'Evidence snapshot: '+ledger['captured_utc']+'.', '',
        '**Replication unit: developmental history.** The core attribute-model evidence represents **three distinct histories (7, 8, 9)**, reused across experiments. Cells, graphs, branches, pulses, perturbation seeds, and timestep/backend repeats are nested measurements. Do not sum histories across rows or report intervention fractions as population success probabilities.', '',
        'Initiation, maintenance, and context sensitivity are separate claims. Numerical acceptance is separate from a scientific outcome; failed gates and unsettled states remain evidence. The history column counts registered attribute-model histories. The historical fate-model developmental-refinement row uses one separate history with five nested numerical repeats; analytic reference problems have no developmental history.', '']
    for theme in ('Initiation','Maintenance','Context','Numerical'):
        lines += ['## '+theme, '', '| Study | Attribute histories | Evidence and current conclusion | Validation |', '|---|---:|---|---|']
        for row in ledger['entries']:
            if row['theme'] != theme: continue
            links=', '.join(f"[{Path(p).parent.name}/{Path(p).name}](../{p})" for p in row['evidence'] if Path(p).name not in ('protocol.json','status.json','verification.json'))
            lines.append(f"| `{row['id']}` | {row['histories']} | {row['claim']} {links} | {row['validation']} |")
        lines += ['']
        for row in ledger['entries']:
            if row['theme'] == theme: lines += [f"- **{row['id']}:** {row['limitations']}"]
        lines += ['']
    lines += ['## History-level context assessment', '',
        '| History | State shifts / challenges | Response shifts / comparisons | Fresh recipient conservative-mixing ratio |',
        '|---|---:|---:|---:|']
    neighbor = next(r for r in ledger['entries'] if r['id']=='neighbor-context')
    for seed,row in neighbor['details']['by_history'].items():
        lines += [f"| {seed} | {row['state_shifts']} / {row['challenges']} | {row['response_shifts']} / {row['pulse_comparisons']} | {row['conservative_recipient_state_ratio']:.3f} |"]
    lines += ['', 'The counts describe nested interventions within each history. The conservative recipient shift repeats across all three histories. Both immediate-pulse nonrecoveries are in history 8; its delayed follow-up is selected and does not add a fourth history.', '',
        '## Unsupported claims and decisions', '',
        '- Chemical differences and donor-nearer responses do not establish autonomous cell types, a prescribed type count, or biological function.',
        '- Frozen chemical local stability is not full moving-system stability. A finite-horizon recovery is not permanent memory or inheritance.',
        '- Mature continuation timestep/backend passes do not replace the failed historical developmental gate or establish attribute-development spatial convergence.',
        '- Conservative exchange does not resolve the failed general geometric conductance closure. Passed flux prototypes are not silently substituted into live results.',
        '- Polarity-dependent initiation suppression has one matched geometry; necessity is conditional on the tested coefficient set. The chemical multiplier remains active.', '',
        '## Reproduce the consolidation', '',
        'Run `python -m embryo.results_ledger` to deliberately refresh the snapshot and `python -m embryo.results_ledger --verify` to check its current evidence/source hashes. Refreshing records report values; it does not rerun assays or strengthen their inference. [Machine-readable ledger](results_ledger.json) contains study IDs, explicit history IDs, details, SHA-256 evidence hashes, and the complete output-directory inventory. Uncurated folders and caches are not promoted to accepted studies by their names or completion status.', '',
        'No new scientific experiment was launched during consolidation. The next scientific decision remains moving-geometry confirmation of the context interventions, after reviewing this ledger and manuscript.', '']
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args(); target = REPO/'docs/results_ledger.json'
    if args.verify:
        ledger = json.loads(target.read_text())
        if 'precision_audit' in ledger:
            from .precision_claim_audit import verify
            verify(ledger)
        for name, expected in {**ledger['evidence_sha256'], **ledger['source_sha256'], **ledger['history_checkpoint_sha256']}.items():
            if digest(REPO/name) != expected: raise ValueError('Ledger evidence changed: '+name)
        if (REPO/'docs/results_ledger.md').read_text() != markdown(ledger):
            raise ValueError('Human and machine-readable ledgers differ')
        print('Ledger evidence, source hashes, and rendered Markdown verified.')
    else:
        ledger = snapshot()
        if target.exists() and 'precision_audit' in json.loads(target.read_text()):
            from .precision_claim_audit import annotate
            ledger = annotate(ledger)
        target.write_text(json.dumps(ledger,indent=2,allow_nan=False)+'\n')
        (REPO/'docs/results_ledger.md').write_text(markdown(ledger))
        print(f"Recorded {len(ledger['entries'])} curated studies; three shared primary histories.")


if __name__ == '__main__': main()
