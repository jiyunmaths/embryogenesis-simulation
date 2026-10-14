"""Write the stopped-pilot assessment and a diagnostic plot from the reviewed evidence."""
from datetime import datetime, timezone
from pathlib import Path
from zoneinfo import ZoneInfo
import hashlib
import json
import shutil

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent.relative_to(Path.cwd().resolve())
STUDY = Path('output/studies/fresh-zygote-pilot')


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    review_path = HERE/'independent_review.json'
    review = json.loads(review_path.read_text())
    assert review['review_completed'] and review['numerical_stop_confirmed'] and not review['qualification_accepted']
    out = Path('docs/fresh_zygote_pilot_assessment.json')
    assert not out.exists(), 'Never overwrite a completed assessment'
    status_path = STUDY/'status.json'; status = json.loads(status_path.read_text())
    stop = datetime.fromtimestamp(status_path.stat().st_mtime, timezone.utc)
    stop_local = stop.astimezone(ZoneInfo('America/New_York'))
    diagnostic = review['cutoff_diagnostic']; peak = diagnostic['original_peak']
    result = dict(assessed_utc=datetime.now(timezone.utc).isoformat(), execution_state='stopped_on_numerical_failure',
        review_completed=True, qualification_accepted=False, developmental_design_completed=False,
        numerical_stop_confirmed=True, protocol_sha256=review['protocol_sha256'],
        registered_independent_histories=3, fully_completed_matched_histories=[1101,1102],
        completed_paths=18, total_paths=24, completed_pairs=9, passed_pairs=8, total_pairs=12,
        stop_time_utc=stop.isoformat(), stop_time_local=stop_local.isoformat(), runtime_seconds=status['seconds'],
        failure=review['pairs'][-1], endpoints=review['endpoints'], unexecuted_cases=review['unexecuted_cases'],
        all_completed_individual_paths_pass_recorded_physical_gates=True,
        source_hashes_verified=review['source_count'], input_hashes_verified=review['input_count'],
        cutoff_diagnostic=diagnostic, review_scope=review['scope'],
        independent_saved_evidence_review=dict(path=str(review_path),sha256=sha(review_path)),
        numerical_quality=dict(max_volume_error=max(x['audit']['max_volume_error'] for x in review['paths']),
            min_effective_radius_voxels=min(x['audit']['min_radius'] for x in review['paths']),
            max_boundary_occupancy=max(x['boundary_max'] for x in review['paths']),
            max_transport_closure_error=max(x['graph_closure_max'] for x in review['paths']),
            max_dilution_amount_error=max(x['audit']['dilution_amount_error'] for x in review['paths'])),
        preserved_history=review['preserved_history'],
        post_cleavage_started=False, biological_identity_qualified=False, general_continuum_qualified=False,
        next_test=dict(history=1103,branch='off',grid=120,duration=18.,candidate_dt=.0001171875,
            reference_dt=.000234375,thresholds_unchanged=True,not_launched=True,
            policy='New separately registered full developmental timestep diagnostic. Retain the original failed pair and all unexecuted controls. No seed replacement or post hoc threshold relaxation.'))
    out.write_text(json.dumps(result,indent=2)+'\n')

    # Plot original measurements and explicitly labeled algebraic sensitivity.
    series = diagnostic['time_series']; x=np.array([v['time'] for v in series])
    y=100*np.array([v['relative_spectrum'] for v in series])
    z=100*np.array([v['posthoc_no_cutoff_relative_spectrum'] for v in series])
    fig,axes=plt.subplots(1,2,figsize=(11,4.2),layout='constrained')
    for ax in axes:
        ax.plot(x,y,color='#b23b32',linewidth=1.3,label='Recorded contact cutoff')
        ax.plot(x,z,color='#306a8b',linewidth=1.2,linestyle='--',label='No-cutoff reassembly (diagnostic)')
        ax.axhline(1.,color='#505050',linewidth=1.,linestyle=':',label='Registered 1% limit')
        ax.scatter([peak['time']],[100*peak['relative_spectrum']],s=30,color='#b23b32',zorder=3)
        ax.set_xlabel('Developmental time');ax.set_ylabel('Relative spectral discrepancy (%)')
        ax.grid(alpha=.2);ax.spines[['top','right']].set_visible(False)
    axes[0].set_xlim(0,18);axes[0].set_ylim(0,1.3);axes[0].set_title('History 1103: feedback off')
    axes[1].set_xlim(11.25,11.6);axes[1].set_ylim(.55,1.3);axes[1].set_title('Contact switches during cleavage')
    axes[1].annotate('15 vs 17 edges\n1.184% at t=11.43',xy=(11.43,100*peak['relative_spectrum']),
        xytext=(11.28,1.20),fontsize=9,arrowprops=dict(arrowstyle='-',color='#505050'))
    axes[0].legend(loc='upper left',fontsize=8)
    fig.suptitle('A localized transport-spectrum failure stops the fresh-zygote pilot',fontsize=13)
    for suffix in ('png','svg'):fig.savefig(HERE/f'cutoff-diagnostic.{suffix}',dpi=170)
    plt.close(fig)

    pair_rows=[]
    for row in review['pairs']:
        pair_rows.append(f"| {row['history']} | {row['branch']} | {100*row['recorded_field_max']:.4f}% | {100*row['independently_recomputed_regular_errors']['relative_spectrum_max']:.4f}% | {100*row['independently_recomputed_dense']['peak']['relative_spectrum']:.4f}% | {'Pass' if row['passed'] else '**Fail**'} |")
    endpoint_rows=[]
    for row in sorted((x for x in review['endpoints'] if x['level']=='development_half'),key=lambda x:(x['history'],x['branch'])):
        endpoint_rows.append(f"| {row['history']} | {row['branch']} | {row['log_activator_sd']:.6f} | {row['axis_ratio']:.5f} | {row['frozen_uniform_unstable_modes']} |")
    edge_rows=[]
    for row in diagnostic['differing_edges']:
        edge_rows.append(f"| {row['cells'][0]}–{row['cells'][1]} | {row['cutoff_margins'][0]:.6g} | {row['cutoff_margins'][1]:.6g} | {row['conductance'][0]:.6g} | {row['conductance'][1]:.6g} |")
    text=f'''# Fresh-zygote pilot: stopped-design assessment

The coordinator **stopped on a numerical failure at {stop_local:%H:%M EDT on %d %B %Y}**, after **18/24 trajectories and 9/12 timestep pairs**. Eight pairs pass; history **1103, feedback off**, fails the dense spectrum gate. Six registered history-1103 control trajectories remain unexecuted. No post-cleavage continuation was started. This is a stopped experiment, not a completed three-history qualification.

[Frozen design](fresh_zygote_pilot.md), [data assessment](fresh_zygote_pilot_assessment.json), [saved-evidence review](../{review_path}), [review scope and log](../{HERE}/review.log). The protocol remains `{review['protocol_sha256']}`. Original cases, failure records, thresholds and unused histories are preserved.

## What passed and what failed

All eighteen individual paths reached t=18 with sixteen cells, fifteen completed abscissions and no live furrows. Their recorded every-step physical, positivity and conservation gates pass. The saved-evidence audit verifies their event/lineage, regular histories, dense graphs, pressure records and evidence hashes; it also checks all actual endpoints and independently reconstructs every density checkpoint in the failed pair. Actual-start prefixes and eight-step live-cleavage restarts remain recorded runtime checks, not new full-horizon CPU/GPU replays.

| History | Branch | Density-field error | Regular spectrum error | Dense spectrum error | Pair |
| --- | --- | ---: | ---: | ---: | --- |
{chr(10).join(pair_rows)}

Both density-field and spectral limits are 1%. The failed pair passes every regular-history comparison and its density-field screen; its independently reconstructed field maximum is **0.2552%**. The separate dense equal-count spectrum maximum is **1.1836% at t=11.43**, against the unchanged **1%** limit. This is the only above-limit observation among that pair's comparable dense records. Dense unequal-count observations at t=9.075 and 11.88 remain reported and excluded from equal-dimension spectral comparison.

Across all completed paths, maximum volume error is **{100*result['numerical_quality']['max_volume_error']:.3f}%** (5% limit), minimum effective radius **{result['numerical_quality']['min_effective_radius_voxels']:.3f} voxels** (four-voxel limit), maximum boundary occupancy **{result['numerical_quality']['max_boundary_occupancy']:.3g}** (0.01 limit), and dilution amount error at most **{result['numerical_quality']['max_dilution_amount_error']:.3g}** (2e-14 limit). Weighted graph conservation/symmetry closure remains near machine precision. These checks do not override the dense spectral failure.

## Diagnosis: two contacts switch across the cutoff

At t=11.43, both trajectories have the same eight cell IDs, the same parent IDs and eight active furrows. Relative volume discrepancy is only **{100*diagnostic['peak_relative_volume_error']:.6f}%**. The coarser trajectory has **15 edges**, its finer partner **17**. Two contact overlaps lie on opposite sides of the hard cutoff (0.02 times the largest overlap):

| Cell pair | Coarser overlap minus cutoff | Finer overlap minus cutoff | Coarser conductance | Finer conductance |
| --- | ---: | ---: | ---: | ---: |
{chr(10).join(edge_rows)}

The finer 9–14 edge is just 2.44e-8 above its cutoff, but the hard selection turns it into a finite conductance rather than a gradual onset. Cleavage timing also differs slightly between timesteps; the diagnostic does not assign the contact changes to thresholding alone.

![Recorded failure and explicitly post hoc sensitivity reassembly](../{HERE}/cutoff-diagnostic.png)

Using the same saved contacts and volumes but an identical contact-edge support in both matrices reduces the peak discrepancy to **{100*diagnostic['matched_support_peak_errors']['intersection']:.4f}%** with the edge intersection, or **{100*diagnostic['matched_support_peak_errors']['union']:.4f}%** with the union. Reassembling conductances without cutoff gives **{100*peak['posthoc_no_cutoff_relative_spectrum']:.4f}%** at this time and **{100*diagnostic['posthoc_no_cutoff_peak']['posthoc_no_cutoff_relative_spectrum']:.4f}%** maximum over equal-count dense observations.

This supports a localized cutoff contribution to the numerical discrepancy rather than a visible volume collapse or broad physical failure. **These are algebraic sensitivity calculations on saved geometries.** No chemical trajectory was rerun with these matrices. They neither change the failed decision nor qualify a no-cutoff model, a smoother contact law or the full moving system.

## Scientific interpretation of the short endpoints

| History | Fine branch | SD of log activator at t=18 | Shape axis ratio | Frozen uniform unstable modes |
| --- | --- | ---: | ---: | ---: |
{chr(10).join(endpoint_rows)}

These values are descriptive endpoints of the stopped pilot. Chemistry remains close to uniform, with log-activator SD around 0.0019–0.0036. They do not establish persistent differentiation, response memory or cell identity. The short developmental window was designed to qualify cleavage before separately gated long chemical observations.

Shape anisotropy is present in feedback-off controls. Polarity makes a small matched change in the shape and transport spectrum in histories 1101 and 1102, but the stopped cohort provides no accepted demonstration of a large polarity-specific shape effect or identity formation. Nine versus eight positive-growth modes refer to linear chemistry on an actual frozen endpoint; evolving geometry, dilution, polarity dynamics and nonlinear maintenance are excluded from that calculation.

The replication unit is the newly seeded history: three histories were registered, two have their complete matched controls, and the third has only its failed feedback-off pair. Paths, cells, branches and repeats are nested. Eight of nine passing pairs is a numerical inventory, not a population success probability. History 1103 cannot be omitted or replaced to promote the cohort.

## Separate document preservation issue

All **38 source hashes and 382 input hashes** match. The review verifies **1,769 unchanged historical entries** and the original archived bytes for five authorized manuscript edits, using the latest [comment-resolution record](../{review['preserved_history']['authorization_record']}). Four current files match that record; `sections.md` has further user edits including five new inline comments and does not match the last resolved version. Its observed annotated bytes are preserved in the review folder; further editorial reconciliation remains pending. The review does not edit or resolve these new comments. These documents are not simulation inputs and did **not** cause the observed stop: the coordinator stopped at the numerical pair gate before reaching its final historical-location check. The original preservation manifest is unchanged.

## Next numerical test

Register a separate **history-1103 feedback-off full developmental repeat at dt=0.0001171875**, retaining grid 120, the same analytic zygote/random streams, chemistry, contact cutoff and every acceptance limit. Compare it with the existing dt=0.000234375 path through t=0–18, including the dense t=11.43 neighborhood, physical fields, events and full recorded carry state. A short restart at the offending time would not test the accumulated cleavage-timing error.

The purpose is to determine whether another timestep halving resolves the contact-spectrum failure. If it does not, diagnose the discontinuous contact-selection law under a separately specified model change and repeat its reference, conservation, spectrum and developmental checks. Keep this pilot failed. Complete the remaining matched histories/controls under an accepted prospective design before any 60/240-unit formation or identity follow-up. No new simulation is launched by this review, and the ledger/manuscript scientific claims are unchanged.
'''
    Path('docs/fresh_zygote_pilot_assessment.md').write_text(text)
    # Keep pre-edit scientific README/notebook bytes alongside this assessment.
    prior=HERE/'before-document-update';prior.mkdir(exist_ok=False)
    for name in ('README.md','docs/lab_notebook.md'):
        dest=prior/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(name,dest)
    (HERE/'assessment-provenance.json').write_text(json.dumps(dict(
        assessed_utc=datetime.now(timezone.utc).isoformat(),original_status_sha256=sha(status_path),
        independent_review_sha256=sha(review_path),assessment_json_sha256=sha(out),
        assessment_markdown_sha256=sha('docs/fresh_zygote_pilot_assessment.md'),
        plot_sha256={name:sha(HERE/name) for name in ('cutoff-diagnostic.png','cutoff-diagnostic.svg')},
        original_failed_decision_retained=True,new_simulation_launched=False),indent=2)+'\n')
    print('Wrote stopped-pilot assessment and cutoff diagnostic. Qualification remains false.')


if __name__=='__main__':main()
