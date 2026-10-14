"""Write a scoped control assessment from the independent completed audit."""
from datetime import datetime, timezone
from pathlib import Path
from zoneinfo import ZoneInfo
import hashlib
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent.relative_to(Path.cwd().resolve())
STUDY=Path('output/studies/fresh-zygote-control-completion')
ORIGINAL=Path('output/studies/fresh-zygote-pilot')
BRANCHES=('off','tension_adhesion','polarity','full')
LABELS=dict(off='Feedback off',tension_adhesion='Tension + adhesion',polarity='Polarity',full='Full combination')
COLORS=dict(off='#465661',tension_adhesion='#c57422',polarity='#267b72',full='#775099')


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def main():
    out=Path('docs/fresh_zygote_control_completion_assessment.json')
    assert not out.exists(), 'Never overwrite a completed assessment'
    review_path=HERE/'independent_review.json';review=json.loads(review_path.read_text())
    assert review['review_completed'] and review['scoped_control_completion_accepted']
    assert review['scoped_combined_developmental_preparations_accepted'] and not review['original_pilot_accepted']
    status=json.loads((STUDY/'status.json').read_text())
    assert status['state']=='completed_provisional' and status['completed']==6
    finished=datetime.fromtimestamp((STUDY/'status.json').stat().st_mtime,timezone.utc)
    local=finished.astimezone(ZoneInfo('America/New_York'))
    endpoints=review['combined_primary_endpoints'];pairs={p['branch']:p for p in review['pairs']}
    total_sd=[e['log_activator_sd'] for e in endpoints];modes=[e['frozen_uniform_unstable_modes'] for e in endpoints]
    assessment=dict(assessed_utc=datetime.now(timezone.utc).isoformat(),review_completed=True,
        scoped_control_completion_accepted=True,scoped_combined_developmental_preparations_accepted=True,
        original_pilot_accepted=False,original_failure_retained=True,
        protocol_sha256=review['protocol_sha256'],finished_utc=finished.isoformat(),finished_local=local.isoformat(),
        runtime_hours=status['seconds']/3600,new_scientific_histories=0,new_numerical_paths=6,reused_numerical_paths=2,
        registered_developmental_histories=[1101,1102,1103],history_branch_preparations=12,
        primary_endpoint_dt=.000234375,refinement_dts_by_history={'1101':[.00046875,.000234375],
            '1102':[.00046875,.000234375],'1103':[.000234375,.0001171875]},
        pairs=review['pairs'],paths=review['paths'],endpoints=review['endpoints'],
        combined_preparations=review['combined_preparations'],combined_primary_endpoints=endpoints,
        source_hashes_verified=review['source_count'],input_hashes_verified=review['input_count'],
        independently_reconstructed_density_checkpoints=len(review['checkpoints']),
        live_checkpoint_states_checked=len(review['live_checkpoints']),
        independent_review=dict(path=str(review_path),sha256=sha(review_path)),review_scope=review['scope'],
        biological_identity_qualified=False,general_continuum_qualified=False,
        general_developmental_convergence_accepted=False,post_cleavage_started=False,
        next_test=dict(question='Qualify efficient mature timesteps on the actual saved developmental endpoints.',
            histories=[1101,1102,1103],mechanical_contexts=list(BRANCHES),
            candidate_dts=[.00375,.001875,.0009375],requires_separate_registration=True,
            preserve='No chemistry reset or supplied pattern; retain geometry, polarity, phase carry, chemical correction, lineage/events and RNG.',
            planned_subsequent_windows=[60.,240.],not_launched=True))
    out.write_text(json.dumps(assessment,indent=2)+'\n')
    fig,axes=plt.subplots(2,3,figsize=(12,6),sharex=True,sharey='row',layout='constrained')
    histories={}
    for col,history in enumerate((1101,1102,1103)):
        for branch in BRANCHES:
            folder=(ORIGINAL if history!=1103 else STUDY)/f'h{history}_{branch}_development_half'
            data=json.loads((folder/'gpu-history.json').read_text());histories[str(folder/'gpu-history.json')]=sha(folder/'gpu-history.json')
            times=[v['time'] for v in data]
            axes[0,col].plot(times,[v['log_activator_sd'] for v in data],color=COLORS[branch],linewidth=1.25,label=LABELS[branch])
            axes[1,col].plot(times,[v['axis_ratio'] for v in data],color=COLORS[branch],linewidth=1.25)
            if history==1103:
                quarter=STUDY/f'h1103_{branch}_development_quarter/gpu-history.json'
                fine=json.loads(quarter.read_text());histories[str(quarter)]=sha(quarter)
                for row,key in ((0,'log_activator_sd'),(1,'axis_ratio')):
                    axes[row,col].plot([v['time'] for v in fine],[v[key] for v in fine],color=COLORS[branch],linewidth=.9,linestyle='--',alpha=.8)
        axes[0,col].set_title(f'History {history}')
        axes[1,col].set_xlabel('Developmental time')
        for ax in axes[:,col]:
            ax.set_xlim(0,18);ax.grid(alpha=.2);ax.spines[['top','right']].set_visible(False)
    axes[0,0].set_ylabel('SD of log activator')
    axes[1,0].set_ylabel('Aggregate shape axis ratio')
    axes[0,0].legend(fontsize=8,loc='upper left')
    fig.suptitle('Fresh zygotes through cleavage: small chemical differences, packing anisotropy',fontsize=13)
    fig.supxlabel('All solid curves: dt=0.000234375. Dashed history-1103 curves: dt=0.0001171875.',fontsize=9)
    for ext in ('png','svg'):fig.savefig(HERE/f'developmental-controls.{ext}',dpi=170)
    plt.close(fig)
    rows=[]
    for branch in BRANCHES:
        p=pairs[branch];rows.append(f"| {LABELS[branch]} | {100*p['field_max']:.4f}% | {100*p['errors']['relative_spectrum_max']:.4f}% | {100*p['dense']['peak']['relative_spectrum']:.4f}% | Pass |")
    endpoint_rows=[]
    for history in (1101,1102,1103):
        for branch in BRANCHES:
            e=next(e for e in endpoints if e['history']==history and e['branch']==branch)
            endpoint_rows.append(f"| {history} | {LABELS[branch]} | {e['log_activator_sd']:.8f} | {e['axis_ratio']:.8f} | {e['frozen_uniform_unstable_modes']} |")
    text=f'''# Finer matched developmental controls: completed assessment

The independent saved-array audit **accepts all three new history-1103 control pairs** under the unchanged gates. Together with its previously accepted feedback-off refinement and the retained passing pairs for histories 1101/1102, **all twelve history–branch preparations now have scoped numerical support through t=18**. This combines separately declared timestep ladders; the original failed pilot remains failed and is not retrospectively accepted.

[Frozen design](fresh_zygote_control_completion.md), [machine-readable assessment](fresh_zygote_control_completion_assessment.json), [independent saved-array audit](../{review_path}), [audit log](../{HERE}/review.log). Current protocol SHA-256: `{review['protocol_sha256']}`. Original coordinator records, including their provisional/pending-review fields, are unchanged; this separate assessment records the independent decision.

## Completion and unchanged limits

The six new paths completed **{local:%d October %Y at %H:%M America/New_York}**, after **{status['seconds']/3600:.2f} hours** (about 5 hours 40 minutes). Every path reaches t=18 with sixteen cells, fifteen completed abscissions, no active furrows and at least three mature units. The new study adds **zero independent histories**: its six paths and two reused feedback-off references are nested within existing history 1103.

All history-1103 pairs compare **dt=0.000234375 versus 0.0001171875** at grid 120³. Density fields and both spectral checks retain their **1%** limits. The accepted feedback-off pair below is reused; only the other three pairs are new.

| Branch | Maximum density-field error | Maximum regular-spectrum error | Maximum dense-spectrum error | Decision |
| --- | ---: | ---: | ---: | --- |
{chr(10).join(rows)}

All chemical distribution, contrast, shape, volume/amount, polarity, maturity-time, physical, boundary/connectivity, conservative graph and cleavage gates pass. Every actual-start prefix and recorded eight-step live-cleavage restart passes. Contact-support switches and unequal-count observations remain visible: a pass does not eliminate the hard contact cutoff or prove a convergence order. The tension/adhesion dense maximum is at t=11.43 with sixteen versus seventeen edges; it remains below the same 1% limit.

## Independent audit and combined preparation scope

The NumPy/SciPy audit calls no production comparison or embryo solver and performs no fresh GPU stepping. It verifies **{review['source_count']} source hashes and {review['input_count']} input hashes**, all eight history-1103 paths' saved history/graph/pressure/event evidence, **56 regular checkpoint states and density reconstructions**, and **24 live-cleavage checkpoint states**. It independently recomputes every history-1103 distribution/shape/volume/amount/polarity comparison and regular/dense graph spectrum.

The earlier review of histories 1101/1102 had recomputed their history and graph comparisons but retained recorded passing density-field errors. This audit additionally rebuilds **all 112 checkpoint states/density fields** for those eight pairs and recomputes their regular/dense comparisons; they pass without changing any gate. Thus the combined statement does not assume independently checked density fields where they were previously only recorded. Earlier pressure/event checks are inherited through their immutable independently reviewed evidence, not replayed mechanically.

There are **three registered fresh histories (1101, 1102, 1103)**, each with four mechanical branches. They share a primary saved endpoint timestep **0.000234375**. Histories 1101/1102 compare 0.00046875/0.000234375; history 1103 compares 0.000234375/0.0001171875 under the separately registered refinement. Report these different ladders explicitly. All twelve preparations pass their scoped timestep screens; three histories remain a pilot, not a population-frequency estimate or universal developmental convergence result.

The original history-1103 coarse comparison still fails at 1.1836%. The original failed status, protocol and six unexecuted control cases remain unchanged. The new runs belong to this separate finer design. Recorded actual-start and short restart tests are verified from saved evidence; no new full fine CPU ladder or full remaining-horizon mechanics replay is claimed.

## What the development shows

The following table uses the **common primary endpoint timestep 0.000234375** rather than mixing finest endpoints across histories.

| History | Branch | SD of log activator at t=18 | Aggregate axis ratio | Frozen uniform unstable modes |
| --- | --- | ---: | ---: | ---: |
{chr(10).join(endpoint_rows)}

![Chemical contrast and shape trajectories across the three fresh histories](../{HERE}/developmental-controls.png)

**Small chemical differences arise during cleavage, but sustained differentiation has not yet been demonstrated.** Endpoint log-activator SD spans **{min(total_sd):.5f}–{max(total_sd):.5f}** across the twelve preparations. Chemistry remains close to uniform. These differences and branch effects are continuous descriptive measurements, not two prescribed identities, autonomous memory, lineage commitment or biological cell types. The short mature window here is insufficient to establish persistent organization.

**Packing and cleavage produce shape anisotropy in every branch, including feedback off.** Polarity changes the endpoint shape only slightly in these runs. This is not evidence for a large chemically selected axis or activator–inhibitor-specific shape organization. The figure shows recorded trajectories rather than independent biological replicas or a parameter sweep.

**The frozen endpoint chemistry still has supported growing modes:** {min(modes)}–{max(modes)} nonuniform modes across these preparations. This calculation holds contact conductances and volumes fixed and linearizes about uniform chemistry. It does not establish that fluctuations will grow during continued moving mechanics, or predict whole-system stability. The longer moving experiment is still needed to test whether contact remodeling changes formation versus maintenance in fresh histories.

## Next test

Register the **mature-timestep qualification on the actual twelve saved primary endpoints**, retaining each history's matched mechanical context. Test efficient mature candidates 0.00375 and 0.001875, with 0.0009375 only if required under a separately declared failure/diagnostic rule. Keep chemistry, geometry, polarity, lineage/event and random state, phase carry and chemical correction; do not reset chemistry or supply a pattern. Require actual-start CPU/GPU prefixes, correction-aware restarts, physical checks and matched timestep comparisons in the mature context.

After independent acceptance of that mature numerical method, begin the planned **60-elapsed-unit formation window**, with the **240-unit continuation** separately gated. Measure first sustained formation, retention/loss and contact-network evolution by history. These developmental qualification runs alone cannot decide initiation-versus-maintenance or biological identity.

No new simulation is launched by this review. Ledger and manuscript scientific snapshots remain unchanged. The README's numerical scope and lab notebook are updated separately, while all original failures remain accessible.
'''
    Path('docs/fresh_zygote_control_completion_assessment.md').write_text(text)
    (HERE/'assessment-provenance.json').write_text(json.dumps(dict(utc=datetime.now(timezone.utc).isoformat(),
        protocol_sha256=sha(STUDY/'protocol.json'),final_status_sha256=sha(STUDY/'status.json'),
        original_failed_status_sha256=sha(ORIGINAL/'status.json'),original_failed_pair_sha256=sha(ORIGINAL/'time-h1103-off.json'),
        independent_review_sha256=sha(review_path),assessment_json_sha256=sha(out),
        assessment_markdown_sha256=sha('docs/fresh_zygote_control_completion_assessment.md'),
        figure_sha256={ext:sha(HERE/f'developmental-controls.{ext}') for ext in ('png','svg')},
        figure_history_sha256=histories,new_simulation_launched=False),indent=2)+'\n')
    print('Wrote accepted scoped control assessment, combined preparation inventory and developmental figure.')


if __name__=='__main__':main()
