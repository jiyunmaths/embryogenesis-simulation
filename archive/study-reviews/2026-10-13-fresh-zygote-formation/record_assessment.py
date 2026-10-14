"""Record scoped acceptance only after the complete independent audit passes."""
from datetime import datetime, timezone
from pathlib import Path
from zoneinfo import ZoneInfo
import hashlib
import json
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = Path.cwd()
STUDY = ROOT / 'output/studies/fresh-zygote-formation'
PREFIX = ROOT / 'docs/fresh_zygote_formation_assessment'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    review_path = HERE / 'independent_review.json'
    review = json.loads(review_path.read_text())
    assert review['review_completed'] and review['scoped_formation_timestep_accepted']
    assert not PREFIX.with_suffix('.json').exists() and not PREFIX.with_suffix('.md').exists()
    summary = json.loads((STUDY / 'summary.json').read_text())
    status = json.loads((STUDY / 'status.json').read_text())
    assert status['state'] == 'completed_provisional' and not status['qualification_accepted']
    pairs = review['pairs']
    maxima = {key: max(x['errors'][key] for x in pairs) for key in pairs[0]['errors']}
    homologous = {key: max(x['homologous_errors'][key] for x in pairs) for key in pairs[0]['homologous_errors']}
    field_max = max(x['field_max'] for x in pairs)
    dense_max = max(x['dense_peak']['relative_spectrum'] for x in pairs)
    physical = {key: (min if key == 'min_radius' else max)(x['audit'][key] for x in review['paths'])
                for key in review['paths'][0]['audit']}
    endpoints = review['fine_endpoint_observations']
    assert len(endpoints) == 12 and len(review['checkpoints']) == 144
    assert all(f['sustained_formation'] and f['late_persistent'] for x in review['paths']
               for f in x['independently_recomputed_formation'])
    complete_local = datetime.fromtimestamp((STUDY / 'summary.json').stat().st_mtime, ZoneInfo('America/New_York'))
    protocol = json.loads((STUDY / 'protocol.json').read_text())
    cfg = protocol['jobs'][0]['config']
    band = sorted(float(v) for v in np.roots([cfg['signal_da'] * cfg['signal_dh'],
                  cfg['signal_beta'] * cfg['signal_da'] - cfg['signal_dh'], cfg['signal_beta']]))
    growth_ranges = []
    for endpoint in endpoints:
        path = STUDY / f"h{endpoint['history']}_{endpoint['branch']}_fine/gpu-history.json"
        # The audited relocation wrapper preserves canonical workspace paths.
        canonical = str(path.absolute())
        assert sha(path) == review['original_files_sha256'][canonical]
        rows = json.loads(path.read_text())
        growth_ranges.append(dict(history=endpoint['history'], branch=endpoint['branch'],
            minimum_growing_modes=min(row['frozen_uniform_unstable_modes'] for row in rows),
            minimum_max_growth=min(row['frozen_uniform_max_growth'] for row in rows)))
    result = dict(review_completed=True, scoped_formation_timestep_accepted=True,
        utc=datetime.now(timezone.utc).isoformat(), protocol_sha256=review['protocol_sha256'],
        independent_review=dict(path=str(review_path.relative_to(ROOT)), sha256=sha(review_path)),
        histories=review['histories'], independent_histories=3, new_scientific_histories=0,
        numerical_paths=24, time_pairs=12, physical_time_window=[18, 78], grid=120, extent=2.24,
        timesteps=[.00375, .001875], observations_per_path=401, dense_graphs_per_path=4001,
        complete_checkpoint_states_reconstructed=144, initial_arrays_losslessly_checked=432,
        maximum_regular_errors=maxima, maximum_homologous_errors=homologous,
        maximum_density_field_error=field_max, maximum_dense_spectrum_error=dense_max,
        recorded_every_step_extrema=physical,
        paired_topology_mismatch_times=sum(len(x['topology_mismatch_times']) for x in pairs),
        fine_endpoint_observations=endpoints,
        frozen_uniform_instability_band=band, fine_frozen_growth_ranges=growth_ranges,
        sustained_and_late_formation_history_counts_by_branch={branch: 3 for branch in ['off', 'tension_adhesion', 'polarity', 'full']},
        primary_threshold=.1, sensitivity_thresholds=[.05, .2], sustained_window=24., late_window=[36., 60.],
        original_scheduling_failure_retained=True, original_developmental_failure_retained=True,
        production_acceptance_fields_unchanged=True, manuscript_and_ledger_unchanged=True,
        maintenance_240_qualified=False, biological_identity_qualified=False,
        general_convergence_accepted=False, new_simulation_launched=False,
        scope=review['scope'], coordinator_seconds=summary['seconds'], completed_local=complete_local.isoformat())
    PREFIX.with_suffix('.json').write_text(json.dumps(result, indent=2) + '\n')
    labels = {'off': 'Feedback off', 'tension_adhesion': 'Tension + adhesion',
              'polarity': 'Polarity tension', 'full': 'Full coupling'}
    numerical_table = '\n'.join(
        f"| {x['history']} | {labels[x['branch']]} | {100*x['field_max']:.6f}% | "
        f"{100*x['dense_peak']['relative_spectrum']:.6f}% | {x['homologous_errors']['chemical_log_max']:.6g} | Pass |"
        for x in pairs)
    scientific_table = '\n'.join(
        f"| {x['history']} | {labels[x['branch']]} | {x['formation'][1]['first_sustained_onset_elapsed']:.2f} | "
        f"{x['end']['log_activator_sd']:.4f} | {x['formation'][1]['late_min_log_sd']:.4f} | "
        f"{100*(-x['conductance_relative_change']):.2f}% | {x['end']['frozen_uniform_unstable_modes']} |"
        for x in endpoints)
    min_late = min(x['formation'][1]['late_min_log_sd'] for x in endpoints)
    drops = [-100*x['conductance_relative_change'] for x in endpoints]
    max_shape_difference = max(abs(next(x for x in endpoints if x['history'] == h and x['branch'] == 'full')['end']['axis_ratio']
                                   - next(x for x in endpoints if x['history'] == h and x['branch'] == 'off')['end']['axis_ratio'])
                               for h in review['histories'])
    text = f'''# Inherited fresh-history formation: completed assessment

The independent saved-array audit **accepts all twelve registered full-horizon timestep comparisons**. In **all three independent histories (1101–1103)**, inherited developmental chemistry forms sustained cell-to-cell differences under every tested mechanical context, including feedback off. The primary threshold 0.1 and both sensitivity thresholds 0.05/0.2 give sustained and late-persistent outcomes at both clocks.

[Registered design](fresh_zygote_formation.md), [assessment data](fresh_zygote_formation_assessment.json), [independent saved-array review](../archive/study-reviews/2026-10-13-fresh-zygote-formation/independent_review.json), [audit scope](../archive/study-reviews/2026-10-13-fresh-zygote-formation/scope.md), [audit log](../archive/study-reviews/2026-10-13-fresh-zygote-formation/review.log). Protocol SHA-256: `{review['protocol_sha256']}`. Original production status/summary remain provisional records; scoped acceptance is recorded separately here. Original failures, sources, thresholds, manuscript and results-ledger snapshots remain unchanged.

## What was tested

Continue the twelve actual developmental preparations from physical **t=18 to 78**, without a chemistry reset, supplied pattern, feedback switch, further division or backend handoff. The t=18 chemistry already contains small inherited daughter fluctuations; it is not exactly uniform. Use dt=0.00375 and 0.001875 at grid 120³ and extent 2.24. These 24 numerical paths and twelve comparisons are nested within **three existing independent histories**, not 24 independent replicates. Each branch retains its own geometry, polarity and chemistry shaped during earlier development.

The mature continuation uses unchanged compensated resident mechanics: PyTorch arrays/matrices and custom CUDA mechanics/geometry/polarity, float32 visible phase with float64 residual carry and float64 compensated chemical integration. Qualified earlier cleavage preparations and this continuation use their separately tested clocks; this is not one uniform-timestep run from fertilization to t=78.

## Formation and mechanical context

Contrast is across-cell SD of **natural-log activator concentration**, a dimensionless chemical-state measurement. Primary onset is the beginning of a fully observed 24-unit interval above 0.1, measured as elapsed time since t=18. Late persistence requires every saved observation during elapsed 36–60 to exceed threshold. Both clocks agree on onset at 0.15-unit observation resolution; all sampled losses/recoveries are reported, and none occurs at the three thresholds in this cohort.

Use finer-clock values for the descriptive table:

| History | Context | Sustained onset, elapsed | Final contrast | Late minimum | Conductance decrease | Endpoint growing modes |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
{scientific_table}

The smallest late-window contrast is **{min_late:.4f}**, above even the 0.2 sensitivity threshold. Small inherited differences therefore grow into sustained nonlinear chemical heterogeneity in all twelve contexts. This supports formation along these qualified developmental histories without fixed A/B labels. It does not by itself establish biological cell identities or their functions.

![Chemical formation and contact-network remodeling in the three inherited developmental histories](../archive/study-reviews/2026-10-13-fresh-zygote-formation/formation.png)

**Polarity tension is not necessary for chemical formation in these preparations.** Feedback-off branches form in 3/3 histories. Polarity-only versus feedback off delays primary onset by 0.45, 0.60 and 0.30 model units; adding polarity to tension/adhesion gives the same three delays. Full coupling versus feedback off delays onset by 0.45, 0.45 and 0.15. These are modest context differences, not universal initiation suppression.

The effect on final contrast has **no consistent sign**: polarity-associated branches have greater contrast than feedback off in histories 1101/1102, but lower contrast in 1103. Every branch still passes sustained formation and late persistence. Do not report polarity as universally amplifying or suppressing differentiation. Baseline polarity evolution remains enabled even when its mechanical tension coupling is zero; this is not a complete polarity-removal experiment.

## Communication remodels without closing the formation window

Total undirected conductance decreases **{min(drops):.2f}–{max(drops):.2f}%**, including substantial loss under feedback off, which retains passive mechanics. Polarity-associated branches lose more conductance in all three histories. Each endpoint nevertheless retains **4–6 growing nonzero chemical modes**, with positive frozen uniform-state maximum growth. The observed remodeling has not removed all supported instability modes over this window.

For these parameters, the frozen uniform chemical instability band is approximately **{band[0]:.3f}<λ<{band[1]:.3f}**. Its boundaries are the roots of D_a D_b λ²+(βD_a−D_b)λ+β=0. Every finer-clock regular snapshot retains at least{min(x['minimum_growing_modes'] for x in growth_ranges)} growing modes; the smallest snapshot maximum growth is{min(x['minimum_max_growth'] for x in growth_ranges):.4f}. Network remodeling moves eigenvalues, but does not close the supported-mode window in the sampled formation trajectories. These are snapshot measurements, not guarantees between observations.

This qualifies the earlier mechanistic narrative. Moving mechanics can suppress initiation in the older mature-reset preparations, and conductance preservation restored formation there. **The present inherited developmental preparations all form patterns while contacts remodel.** Network remodeling is consequential, but its amount alone does not determine failure; the supported spectrum, state, preparation and timing matter. The cohorts and preparations differ, so this comparison does not isolate the older chemistry reset as the cause of suppression or invalidate those scoped mature results.

Frozen growth rates linearize uniform chemistry while holding current geometry fixed; they do not prove whole-system stability, causal mediation or endpoint bistability. Tissue axis ratios remain approximately1.40–1.42 in all contexts; the largest full-versus-off endpoint difference is only{max_shape_difference:.4f}. An activator–inhibitor-specific tissue axis or large-scale shape organization is not established.

## Numerical acceptance

All 24 paths,24 actual-start prefixes and twelve complete timestep comparisons pass unchanged gates. Largest independently reconstructed differences:

- Native occupancy/chemical/polarity density fields: **{100*field_max:.6f}%**, against 1%.
- Dense graph spectrum: **{100*dense_max:.6f}%**, against 1%.
- Homologous chemical log concentration: **{homologous['chemical_log_max']:.6g}**, against 0.01.
- Homologous transport-matrix norm: **{100*homologous['relative_transport_max']:.6f}%**, against 1%.

| History | Context | Maximum native-field error | Maximum dense-spectrum error | Maximum cellwise log error | Decision |
| --- | --- | ---: | ---: | ---: | --- |
{numerical_table}

Five isolated sampled contact-support mismatches are retained: history1101 feedback-off at t=40.35; history1102 polarity at43.515; history1103 feedback-off at36.54; history1103 polarity at28.47 and75.78. The largest field and chemical discrepancies occur in history1103 polarity and remain within the registered gates. Hard-cutoff discontinuities remain visible; matched classifications do not imply identical contact topology or exact trajectories.

Recorded every-step extrema pass: maximum volume error **{100*physical['max_volume_error']:.4f}%** against 5%, minimum equivalent radius **{physical['min_radius']:.4f} voxels** against 4, zero clipping and dilution amount error **{physical['dilution_amount_error']:.3g}** against 2e−14. All paths retain sixteen cells, fifteen inherited abscissions and no active furrows, with passing observed boundary/connectivity checks and exact eight-step compensation-aware replay at elapsed 12 and 36.

## Independent audit and limits

The NumPy/SciPy audit imports no embryo solver or production comparison routine and verifies **63 source pins and 523 input pins**, the accepted parent review, runtime binary and every saved result/prefix archive. It reconstructs **9,624 regular observations**, **96,024 dense conservative graphs**, all distribution/homologous comparisons and chemical growth rates, and **144 complete checkpoint states/native fields**. Generalized eigenproblems check conservation, weighted self-adjointness, operator sign and recorded spectra. Chemical outcomes are recalculated from concentrations rather than copied from the formation labels.

All 432 nonmetadata arrays in 24 retimed starts match the accepted source arrays byte for byte, and every initial full checkpoint exactly matches its retimed start. Only the declared Config step horizon changes. Targets, IDs/parents, cycle due times, lineage/events and random states remain unchanged through every saved checkpoint. Each restart-start archive is byte-identical to its corresponding full checkpoint.

CPU/GPU prefixes, every-step extrema and exact eight-step restarts are verified from immutable saved evidence, **not independently rerun trajectories**. Sampled volumes, radii, cell connectedness, boundary occupancy and numerical carry/correction are independently reconstructed. Acceptance is a two-clock comparison at one spatial resolution in these histories for 60 elapsed units; it does not establish asymptotic convergence order, a full fine CPU horizon, continuum transport accuracy or universal developmental behavior. Original scheduling and developmental failures remain failed.

## Decision and next priority

Accept dt=0.00375/0.001875 as the **scoped 60-unit formation comparison** in these twelve inherited contexts. No additional formation rerun or automatic timestep fallback is required by the registered criteria.

Next, separately register the **240-total-elapsed-unit maintenance test from t=18**, retaining the existing formation states and each branch's inherited chemistry/mechanics, with its own full-window timestep and physical gates. This means continuing the elapsed60 endpoints by 180 units to physical t=258, not adding 240 units. Formation does not qualify that longer clock window. Endpoint coexistence and perturbation/context response tests must follow before these differences are interpreted as candidate identities, commitment or memory.

The coordinator finished **{complete_local.strftime('%d %B %Y at %H:%M:%S')} America/New_York**, in **{summary['seconds']/3600:.3f}hours** (about 10 hours 5 minutes), consistent with the updated forecast. The audit took **{review['seconds']/60:.1f}minutes**. The scientific GPU is idle; this review launches no new simulation. README scope and the lab notebook are updated separately; manuscript and ledger remain unchanged.
'''
    PREFIX.with_suffix('.md').write_text(text)
    print(json.dumps({'scoped_formation_timestep_accepted': True, 'assessment': str(PREFIX.with_suffix('.md'))}, indent=2))


if __name__ == '__main__':
    main()
