"""Write the scoped completed assessment from the independently checked arrays."""
from datetime import datetime, timezone
from pathlib import Path
from zoneinfo import ZoneInfo
import hashlib
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent.relative_to(Path.cwd().resolve())
STUDY = Path('output/studies/fresh-zygote-time-refinement')
PARENT = Path('output/studies/fresh-zygote-pilot')


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def main():
    out=Path('docs/fresh_zygote_time_refinement_assessment.json')
    assert not out.exists(), 'Never overwrite a completed assessment'
    review_path=HERE/'independent_review.json';review=json.loads(review_path.read_text())
    assert review['review_completed'] and review['scoped_refinement_accepted'] and not review['original_cohort_accepted']
    # Both older dense traces still match the preceding saved-array audit exactly.
    prior=json.loads(Path('archive/study-reviews/2026-10-11-fresh-zygote-stopped/independent_review.json').read_text())
    prior_trace_hashes={}
    for case in ('h1103_off_development','h1103_off_development_half'):
        path=PARENT/case/'network.jsonl';expected=prior['original_files_sha256'][str(path.resolve())]
        assert sha(path)==expected
        prior_trace_hashes[str(path)]=expected
    status_path=STUDY/'status.json';status=json.loads(status_path.read_text())
    finished=datetime.fromtimestamp(status_path.stat().st_mtime,timezone.utc)
    local=finished.astimezone(ZoneInfo('America/New_York'))
    pair=review['pair'];dense=pair['dense'];old=review['original_dense_failure'];peak=dense['peak']
    assessment=dict(assessed_utc=datetime.now(timezone.utc).isoformat(),review_completed=True,
        scoped_refinement_accepted=True,original_cohort_accepted=False,
        qualification_scope=review['qualification_scope'],protocol_sha256=review['protocol_sha256'],
        execution_completed=True,finished_utc=finished.isoformat(),finished_local=local.isoformat(),
        seconds=status['seconds'],new_scientific_histories=0,new_numerical_paths=1,reused_numerical_paths=1,
        pair=pair,paths=review['paths'],endpoints=review['endpoints'],
        original_dense_failure=old,prior_failure_time_comparison=review['prior_failure_time_comparison'],
        new_peak_cutoff_edges=review['new_peak_cutoff_edges'],original_failure_retained=True,
        prior_dense_trace_hashes_verified=prior_trace_hashes,
        source_hashes_verified=review['source_count'],input_hashes_verified=review['input_count'],
        independent_review=dict(path=str(review_path),sha256=sha(review_path)),review_scope=review['scope'],
        biological_identity_qualified=False,general_developmental_convergence_accepted=False,
        general_continuum_qualified=False,post_cleavage_started=False,
        next_test=dict(history=1103,branches=['tension_adhesion','polarity','full'],
            dts=[.000234375,.0001171875],grid=120,duration=18.,new_numerical_paths=6,
            new_independent_histories=0,not_launched=True,
            policy='Register a separate matched-control completion study. Retain the original unexecuted design and failed status; reuse accepted feedback-off paths read-only. No original pilot resumption or retrospective cohort acceptance. Long continuation remains separately gated.'))
    out.write_text(json.dumps(assessment,indent=2)+'\n')
    series=review['spectral_time_series'];x=np.array([v['time'] for v in series])
    old_y=np.array([100*v['old_relative_spectrum'] if v['old_relative_spectrum'] is not None else np.nan for v in series])
    new_y=np.array([100*v['new_relative_spectrum'] if v['new_relative_spectrum'] is not None else np.nan for v in series])
    fig,axes=plt.subplots(1,2,figsize=(11,4.2),layout='constrained')
    for ax in axes:
        ax.plot(x,old_y,color='#b23b32',linewidth=1.3,label='Original pair (retained failure)')
        ax.plot(x,new_y,color='#267b72',linewidth=1.3,label='Finer pair (scoped pass)')
        ax.axhline(1.,color='#555555',linestyle=':',linewidth=1.,label='Unchanged 1% limit')
        ax.set_xlabel('Developmental time');ax.set_ylabel('Relative spectral discrepancy (%)')
        ax.set_ylim(0,1.3);ax.grid(alpha=.2);ax.spines[['top','right']].set_visible(False)
    axes[0].set_xlim(0,18);axes[0].set_title('Whole developmental interval')
    axes[0].legend(loc='upper left',fontsize=8)
    axes[1].set_xlim(11.25,11.60);axes[1].set_title('Previous failure neighborhood')
    axes[1].scatter([11.43],[100*old['peak']['relative_spectrum']],s=24,color='#b23b32',zorder=3)
    axes[1].annotate('Original maximum: 1.184%',xy=(11.43,100*old['peak']['relative_spectrum']),
        xytext=(11.27,1.20),fontsize=9,arrowprops=dict(arrowstyle='-',color='#555555'))
    axes[1].scatter([peak['time']],[100*peak['relative_spectrum']],s=24,color='#267b72',zorder=3)
    axes[1].annotate('Finer maximum: 0.494%',xy=(peak['time'],100*peak['relative_spectrum']),
        xytext=(11.46,.70),fontsize=9,arrowprops=dict(arrowstyle='-',color='#555555'))
    fig.suptitle('History 1103, feedback off: another timestep halving passes the registered screen',fontsize=12)
    for extension in ('png','svg'):fig.savefig(HERE/f'timestep-refinement.{extension}',dpi=170)
    plt.close(fig)
    limits=json.loads((STUDY/'protocol.json').read_text())['time_limits']
    metrics=[('Reconstructed density fields',pair['field_max'],limits['field_max'],True),
        ('Regular transport spectra',pair['errors']['relative_spectrum_max'],.01,True),
        ('Dense equal-count transport spectra',peak['relative_spectrum'],.01,True),
        ('Log chemical distribution distance',pair['errors']['log_wasserstein_max'],.01,False),
        ('Chemical contrast difference',pair['errors']['contrast_abs_max'],.005,False),
        ('Shape-axis discrepancy',pair['errors']['relative_axis_max'],.01,True),
        ('Total-volume discrepancy',pair['errors']['relative_total_volume_max'],.005,True),
        ('Total-species-amount discrepancy',pair['errors']['whole_amount_relative_max'],.005,True),
        ('Mean-polarity discrepancy',pair['errors']['mean_polarity_abs_max'],.01,False),
        ('Maturity-time discrepancy',pair['maturity_time_error'],.05,False)]
    rows=[]
    for name,value,limit,percent in metrics:
        display=f'{100*value:.4f}%' if percent else f'{value:.7g}'
        bound=f'{100*limit:g}%' if percent else f'{limit:g}'
        rows.append(f'| {name} | {display} | {bound} | Pass |')
    endpoints=[]
    for e in review['endpoints']:endpoints.append(f"| {e['dt']} | {e['log_activator_sd']:.8f} | {e['axis_ratio']:.8f} | {e['frozen_uniform_unstable_modes']} |")
    text=f'''# History-1103 timestep refinement: completed assessment

The independent saved-array review **accepts the registered finer timestep pair** for history 1103, feedback off, grid 120³, t=0–18. Maximum dense spectral discrepancy is **{100*peak['relative_spectrum']:.4f}%**, below the unchanged **1%** limit. Reconstructed density-field discrepancy is **{100*pair['field_max']:.4f}%**. Every inherited gate passes. This is scoped numerical acceptance of an existing history, not completion or retrospective acceptance of the original three-history pilot.

[Frozen diagnostic design](fresh_zygote_time_refinement.md), [machine-readable assessment](fresh_zygote_time_refinement_assessment.json), [independent review](../{review_path}), [review log](../{HERE}/review.log). The new protocol is `{review['protocol_sha256']}`. Original production records remain unmodified, including their provisional/pending-review fields; this separate report records the independent decision.

## Design and completion

One new trajectory repeated history 1103 from its original analytic zygote at **dt=0.0001171875** (153,600 steps). The saved **dt=0.000234375** trajectory is reused read-only. Only timestep and integer clocks differ; geometry, chemistry, polarity, numerical corrections, lineage, cell-cycle schedule and random streams match at t=0. The conservative transport, contact cutoff 0.02, solver/kernel binary, model and limits are unchanged.

The coordinator finished **{local:%d October %Y at %H:%M America/New_York}**, after **{status['seconds']/60:.1f} minutes**. Both paths reach t=18 with sixteen cells, fifteen abscissions, no active furrows and at least three mature units. Zero new independent histories are added. Branches, cells and timestep repeats remain nested within history 1103.

## Numerical evidence

| Comparison | Maximum | Limit | Decision |
| --- | ---: | ---: | --- |
{chr(10).join(rows)}

The dense maximum occurs at **t={peak['time']}**, with eight cells and **{peak['edges'][0]} versus {peak['edges'][1]} edges**. Regular histories have one unequal-count observation at t=8.7; dense records have four at t=5.955, 8.7, 8.745 and 11.925. These remain recorded and excluded only from equal-dimension spectral comparisons. Edge-support mismatches at t=9.39, 9.495, 11.445 and 15.915 are retained.

The new path's maximum volume error is **1.7112%**, minimum effective radius **8.4995 voxels**, dilution amount error at most **4.44e-16**, and clipping zero. All event, boundary/connectivity and conservative graph checks pass. The exact eight-step live-cleavage restarts at 2→4, 4→8 and 8→16 and the independent actual-start prefix pass their recorded gates.

The independent NumPy/SciPy review imports no embryo solver or production comparison. It verifies **43 source and 405 input hashes**, both paths' events/lineage, regular histories, all dense graph and pressure records, all fourteen regular checkpoint arrays and all six saved live-cleavage states. It independently recomputes distribution/shape/volume/amount/polarity comparisons, graph spectra and all seven paired density-field reconstructions. Actual-start and eight-step restart results are verified from saved evidence; no new GPU stepping, full fine CPU ladder or full remaining-horizon replay is performed.

## What changed at the earlier failure

At **t=11.43**, the original dt=0.00046875/0.000234375 pair remains a **1.1836% failure**, with fifteen versus seventeen contact edges. The new dt=0.000234375/0.0001171875 pair has **seventeen edges in both paths** and **0.2715%** spectral discrepancy at that same time.

![Original failure retained and finer comparison accepted under the same limit](../{HERE}/timestep-refinement.png)

The new maximum at t=11.445 still includes a contact-cutoff switch: cells 10 and 13 have overlap-minus-cutoff margins −7.55e-8 and +2.56e-7, giving zero versus 0.0082056 conductance. Refinement reduces the measured discrepancy sufficiently to pass the screen; it does **not** remove the hard-cutoff discontinuity or demonstrate a convergence order, uniform accuracy across histories, general geometric closure, or full coupled-system stability. Peak times differ between the two comparisons, so their ratio is not an error-order estimate.

Both original dense traces match their preceding independent-review hashes. The original failed pilot/status/comparison and its six unexecuted controls remain unchanged; no failed point, seed, threshold or mechanism has been replaced.

## Scientific scope of the short endpoint

| Timestep | SD of log activator at t=18 | Aggregate axis ratio | Frozen uniform unstable modes |
| --- | ---: | ---: | ---: |
{chr(10).join(endpoints)}

Chemistry is still near uniform. The closely agreeing endpoint measurements are consistent with numerical reproducibility of this short feedback-off development. They establish neither sustained chemical organization nor emergent biological identities. Shape anisotropy occurs without regulatory mechanical feedback; it is not evidence for chemically driven shape symmetry breaking. The nine frozen unstable modes concern local linear chemistry on a fixed endpoint, not the evolving coupled system or eventual nonlinear formation.

## Next test

Register a **separate history-1103 matched-control completion study** with tension/adhesion, polarity and full coupling, each at **dt=0.000234375 and 0.0001171875** through t=18 on grid 120³: **six new numerical paths**, zero new independent histories. Reuse the accepted feedback-off paths only for their declared comparison/control role. Preserve the original pilot's unexecuted design rather than resuming or overwriting it, and retain exact initial-state matching, correction-aware restarts and every physical/temporal limit.

After the matched developmental controls pass their independent review, qualify efficient mature timesteps on actual endpoints before the declared long formation/maintenance windows. No new simulation is launched by this review. Ledger and manuscript scientific claims are unchanged.
'''
    Path('docs/fresh_zygote_time_refinement_assessment.md').write_text(text)
    (HERE/'assessment-provenance.json').write_text(json.dumps(dict(utc=datetime.now(timezone.utc).isoformat(),
        protocol_sha256=sha(STUDY/'protocol.json'),final_status_sha256=sha(status_path),
        original_failed_status_sha256=sha(PARENT/'status.json'),original_failed_pair_sha256=sha(PARENT/'time-h1103-off.json'),
        independent_review_sha256=sha(review_path),assessment_json_sha256=sha(out),
        assessment_markdown_sha256=sha('docs/fresh_zygote_time_refinement_assessment.md'),
        figure_sha256={ext:sha(HERE/f'timestep-refinement.{ext}') for ext in ('png','svg')},
        new_simulation_launched=False),indent=2)+'\n')
    print('Wrote scoped accepted assessment and comparison figure. Original pilot remains failed.')


if __name__=='__main__':main()
