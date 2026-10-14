"""Record a separate assessment only after successful saved-evidence review."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = Path('output/studies/fresh-zygote-mature-clock-fix')
PREFIX = Path('docs/fresh_zygote_mature_clock_fix_assessment')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    review_path = HERE / 'independent_review.json'
    review = json.loads(review_path.read_text())
    summary = json.loads((ROOT / 'summary.json').read_text())
    status = json.loads((ROOT / 'status.json').read_text())
    assert review['review_completed'] and review['scoped_mature_timestep_accepted']
    assert len(review['pairs']) == 12 and len(review['checkpoints']) == 72
    assert all(x['independently_passed'] for x in review['pairs'])
    assert status['state'] == 'completed_provisional' and not status['qualification_accepted']
    assert not PREFIX.with_suffix('.json').exists() and not PREFIX.with_suffix('.md').exists()
    pairs = review['pairs']
    maxima = {k: max(x['errors'][k] for x in pairs) for k in pairs[0]['errors']}
    homologous = {k: max(x['homologous_errors'][k] for x in pairs) for k in pairs[0]['homologous_errors']}
    field_max = max(x['field_max'] for x in pairs)
    dense_max = max(x['dense_peak']['relative_spectrum'] for x in pairs)
    audit = {k: (min if k == 'min_radius' else max)(x['audit'][k] for x in review['paths'])
             for k in review['paths'][0]['audit']}
    endpoints = review['fine_endpoint_observations']
    result = dict(review_completed=True, scoped_mature_timestep_accepted=True,
        utc=datetime.now(timezone.utc).isoformat(), protocol_sha256=review['protocol_sha256'],
        independent_review=dict(path=str(review_path.relative_to(Path.cwd())), sha256=sha(review_path)),
        histories=review['histories'], independent_histories=3, new_scientific_histories=0,
        numerical_paths=24, time_pairs=12, physical_time_window=[18, 24], grid=120, extent=2.24,
        timesteps=[.00375, .001875], observations_per_path=41, dense_graphs_per_path=401,
        complete_checkpoint_states_reconstructed=72, maximum_regular_errors=maxima,
        maximum_homologous_errors=homologous, maximum_density_field_error=field_max,
        maximum_dense_spectrum_error=dense_max, recorded_every_step_extrema=audit,
        paired_topology_mismatch_times=sum(len(x['topology_mismatch_times']) for x in pairs),
        fine_endpoint_observations=endpoints,
        original_scheduling_failure_retained=True, original_developmental_failure_retained=True,
        production_acceptance_fields_unchanged=True, manuscript_and_ledger_unchanged=True,
        long_window_qualification_accepted=False, biological_identity_qualified=False,
        general_convergence_accepted=False, new_simulation_launched=False,
        scope=review['scope'], coordinator_seconds=summary['seconds'])
    PREFIX.with_suffix('.json').write_text(json.dumps(result, indent=2) + '\n')
    labels = {'off': 'Feedback off', 'tension_adhesion': 'Tension + adhesion',
              'polarity': 'Polarity', 'full': 'Full combination'}
    numerical_table = '\n'.join(
        f"| {x['history']} | {labels[x['branch']]} | {100*x['field_max']:.6f}% | "
        f"{100*x['dense_peak']['relative_spectrum']:.6f}% | Pass |" for x in pairs)
    scientific_table = '\n'.join(
        f"| {x['history']} | {labels[x['branch']]} | {x['start']['log_activator_sd']:.6f} → "
        f"{x['end']['log_activator_sd']:.6f} | {x['log_contrast_growth']:.3f} | "
        f"{100*(-x['conductance_relative_change']):.3f}% | "
        f"{x['start']['frozen_uniform_unstable_modes']} → {x['end']['frozen_uniform_unstable_modes']} |"
        for x in endpoints)
    text = f'''# Mature-timestep qualification: completed assessment

The independent saved-array audit **accepts all twelve registered timestep comparisons**. The actual developmental endpoints of histories 1101–1103 can be continued from physical t=18 to24 at **dt=0.00375 versus 0.001875**, at grid120 and extent2.24, under all four tested mechanical contexts. This qualifies the short mature clock comparison; long-window accuracy and persistent organization remain untested.

[Registered design](fresh_zygote_mature_qualification.md), [scheduler correction](fresh_zygote_mature_clock_fix.md), [assessment data](fresh_zygote_mature_clock_fix_assessment.json), [independent saved-array review](../archive/study-reviews/2026-10-12-fresh-zygote-mature/independent_review.json), [audit scope](../archive/study-reviews/2026-10-12-fresh-zygote-mature/scope.md), [audit log](../archive/study-reviews/2026-10-12-fresh-zygote-mature/review.log). Protocol SHA-256: `{review['protocol_sha256']}`. Original production status/summary, sources, thresholds and previous failures remain unchanged. Acceptance is recorded here separately from the coordinator's provisional fields.

## Numerical result

All **24 paths**, **24 actual-start prefixes** and **12 paired comparisons** pass. Twenty-three prefixes are newly executed; the identical, already passing first prefix is reused with explicit provenance from the retained scheduling failure. These paths are nested within **three existing independent histories**, each with four mechanical branches and two timesteps; there are zero new histories. Actual chemistry, geometry, polarity, numerical carry/correction, lineage and random states are inherited without a chemical reset, supplied pattern or mechanical switch.

The largest independently reconstructed differences are:

- Native occupancy/chemical/polarity density fields: **{100*field_max:.6f}%**, against the unchanged 1% gate.
- Dense graph spectrum: **{100*dense_max:.6f}%**, against the unchanged 1% gate.
- Homologous transport-matrix relative norm: **{100*homologous['relative_transport_max']:.6f}%**, against 1%.
- Homologous chemical log concentration: **{homologous['chemical_log_max']:.6g}**, against 0.01.

| History | Branch | Maximum density-field error | Maximum dense-spectrum error | Decision |
| --- | --- | ---: | ---: | --- |
{numerical_table}

No between-timestep contact-support mismatch occurs at the 401 dense observations in each path. Contacts still change **within** trajectories; matched support between clocks does not mean a frozen graph or remove the hard cutoff.

The recorded every-step extrema satisfy all physical limits: maximum volume error **{100*audit['max_volume_error']:.4f}%** (<5%), minimum equivalent radius **{audit['min_radius']:.4f} voxels** (≥4), zero clipping and maximum dilution-stage amount error **{audit['dilution_amount_error']:.3g}** (≤2e-14). Every path retains sixteen cells and fifteen completed cleavage events, with no active furrows. All recorded boundary/connectivity checks and correction-aware eight-step restarts pass. The final t=24 is reached in every case.

## Independent checks and limits

The NumPy/SciPy audit imports no embryo solver and uses no production comparison routines. It verifies **58 source pins and 494 input pins**, the runtime binary, original failures and all saved result evidence. It reconstructs **984 regular observations**, **9,624 dense conservative graph checks**, all weighted-distribution and homologous comparisons, and **72 complete checkpoint states/native density fields** at t=18,21,24. Generalized eigenproblems independently check conservation, weighted self-adjointness and operator sign. Initial checkpoints match their retimed sources byte for byte; targets, IDs/parents, cycle due times, lineage/events and random states remain unchanged. Every t=21 restart archive is byte-identical to its corresponding full checkpoint.

The eight-step CPU/GPU prefixes, every-step extrema and exact restart replays are production checks verified through immutable saved evidence, **not newly rerun CPU/GPU trajectories**. Sampled checkpoint volume/radius/connectedness are independently reconstructed. The audit-only rounding-diagnostic correction is retained in the review archive; no production result or registered gate changed.

This is a two-clock agreement screen at one spatial resolution, in twelve actual endpoint contexts, for six elapsed units. It does not establish an asymptotic convergence order, a full fine CPU horizon, many-cell continuum transport, all mature states or long-window accuracy. The original elapsed-zero scheduling failure and original history-1103 coarse developmental spectral failure remain failed. No ledger/manuscript claim is retroactively promoted.

## What the short continuation shows

Use fine-clock endpoints throughout the descriptive table. Chemical contrast is the across-cell SD of **natural-log activator**. Conductance is summed once over undirected edges; the percentage compares each branch with its own t=18 value.

| History | Branch | Log-activator SD: t=18 →24 | Growth factor | Conductance decrease | Frozen uniform growing modes |
| --- | --- | ---: | ---: | ---: | ---: |
{scientific_table}

![Inherited chemistry and contact conductances during the six-unit mature window](../archive/study-reviews/2026-10-12-fresh-zygote-mature/mature-short-window.png)

**Small chemical differences amplify in every tested branch**, including feedback off. Log contrast increases about **4.98–5.78-fold**, ending at **0.0107–0.0209**. This is early amplification, not established cell identities, sustained nonlinear patterning, commitment or memory.

**Contact conductances decline in every branch**, by **9.56–12.86%**. Passive mechanics therefore also remodels communication. Polarity branches show a somewhat larger decline and modestly smaller chemical contrast than their matched branches in this short interval, but their initial geometry and chemistry already differ because their mechanical contexts were active during development. These descriptive differences do not isolate a mature polarity switch or prove causal mediation.

The endpoint graphs retain **8–9 supported growing modes** at both ends. Conductance loss has not eliminated all unstable modes by t=24. These counts linearize uniform chemistry with fixed volumes/contacts; they are not a stability calculation for the entire moving system. Shape remains anisotropic in every context, including feedback off; no activator–inhibitor-specific tissue axis is established here.

## Decision and next test

Use dt=0.00375 as a **qualified short-window mature candidate**, retaining dt=0.001875 as its refinement comparison. No additional dt=0.0009375 diagnostic is required by this experiment.

Next, separately register the planned **60-elapsed-unit formation test** from the twelve inherited developmental endpoints, retaining each branch's chemistry, geometry, polarity and compensation and testing the full horizon at matched clocks. Measure sustained onset, continued growth or loss, and conductance/spectrum evolution by history. A numerically qualified uniform branch is a scientific outcome. Do not infer long-window accuracy from this six-unit pass or reset the mature chemistry. If the formation pilot and its own numerical gates pass, assess maintenance through the separately gated **240-unit total continuation**. Behavioral/context assays remain necessary before interpreting chemical patterns as candidate identities.

The coordinator finished **12 October 2026 at 10:42:55 America/New_York**, in **{summary['seconds']/3600:.3f} hours** (about 2 hours 8 minutes). This review launches no new simulation. README numerical scope and the lab notebook are updated separately; manuscript and results-ledger snapshots remain unchanged.
'''
    PREFIX.with_suffix('.md').write_text(text)
    print(json.dumps({'scoped_mature_timestep_accepted': True,
                      'assessment': str(PREFIX.with_suffix('.md'))}, indent=2))


if __name__ == '__main__':
    main()
