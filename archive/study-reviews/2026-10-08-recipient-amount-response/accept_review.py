"""Archive accepted raw-array review; no simulations or changed thresholds."""
import argparse,hashlib,json
from pathlib import Path
from datetime import datetime,timezone
import numpy as np

def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
def read(f):return json.loads(Path(f).read_text())
def write(f,v):Path(f).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def main(phase):
 stem='phase_carry_recipient_amount_'+('response' if phase=='immediate' else 'delayed_response')
 name='recipient-amount-'+('response' if phase=='immediate' else 'delayed-response')
 root=Path('output/studies/phase-carry-'+name);review=Path('archive/study-reviews/2026-10-08-'+name);launch=Path('archive/study-launches/2026-10-08-'+name)
 p=read(root/'protocol.json');r=read(review/'independent_review.json');summary=read(root/'summary.json');status=read(root/'status.json')
 assert r['passed'] and r['inputs_unchanged'] and summary['passed'] and status['state']=='completed' and r['protocol_sha256']==sha(root/'protocol.json')
 for collection in (p['source_sha256'],p['input_sha256'],r['inputs_sha256']):
  for f,h in collection.items():assert sha(f)==h,f
 snap=Path('archive/study-document-snapshots/2026-10-08-before-recipient-amount-response-launch/manifest.json');s=read(snap)
 guards={**s['unchanged_guards'],**{f:h for f,h in s['files_sha256'].items() if f not in ('README.md','docs/plan.md')}}
 if phase=='delayed':
  prior=read('docs/phase_carry_recipient_amount_response_assessment.json');prov=read(Path(prior['review_archive'])/'review-provenance.json')
  guards.update({f:h for f,h in prov['files_sha256'].items() if f not in ('README.md','docs/plan.md')})
 for f,h in guards.items():assert sha(f)==h,f
 rows=r['scientific_comparisons'];histories=r['history_outcomes'];ranges=[]
 for seed,cell in [(7,22),(7,23),(8,19),(8,20),(9,29),(9,30)]:
  selected=[x for x in rows if (x['seed'],x['cell'])==(seed,cell)];vals=[x['full_effects']['fine'] for x in selected]
  ranges.append(dict(history=seed,cell=cell,min=min(vals),max=max(vals),all_above=all(x['above_screen'] for x in selected),all_below=all(not x['above_screen'] for x in selected)))
 metrics=[m for x in rows for m in x['fine_pulse_metrics'].values()]
 recovery={k:dict(recovered=sum(m[k] is not None for m in metrics),logical_pulse_labels=len(metrics),min=min((m[k] for m in metrics if m[k] is not None),default=None),max=max((m[k] for m in metrics if m[k] is not None),default=None)) for k in ('target_recovery_time','network_recovery_time')}
 amounts=r['amount_evolution'];amount_diagnostic=dict(pairs=len(amounts),maximum_initial_relative_amount_error=max(x['initial_paired_amount_relative_error'] for x in amounts),final_global_ratio_range=[np.min([x['final_global_species_amount_ratios'] for x in amounts],axis=0).tolist(),np.max([x['final_global_species_amount_ratios'] for x in amounts],axis=0).tolist()],maximum_global_log_separation=np.max([x['maximum_global_species_log_amount_separation'] for x in amounts],axis=0).tolist(),full_curves=str((review/'independent_review.json').resolve()))
 elapsed=read(root/'timing.json');elapsed['gpu_production_seconds']=sum(x['wall_seconds'] for x in r['stage_checks'])
 compact=dict(full_study_accepted=True,phase=phase,protocol_sha256=sha(root/'protocol.json'),original_summary_sha256=sha(root/'summary.json'),original_status_sha256=sha(root/'status.json'),counts=r['counts'],gate_counts=r['gate_counts'],maxima=r['maxima'],source_pin_count=len(p['source_sha256']),input_pin_count=len(p['input_sha256']),immutable_consumed_artifacts=len(r['inputs_sha256']),history_outcomes=histories,scientific_comparisons=rows,recipient_effect_ranges=ranges,baseline_descriptions=r['baseline_descriptions'],recovery=recovery,amount_diagnostic=amount_diagnostic,recorded_metrics_compared=r['recorded_metrics_compared'],maximum_recorded_metric_difference=r['maximum_recorded_metric_difference'],timing=elapsed,independent_review=str((review/'independent_review.json').resolve()),independent_review_sha256=sha(review/'independent_review.json'),review_archive=str(review),prelaunch_snapshot=str(snap),scope=p['scope'])
 if phase=='delayed':compact.update(immediate_vs_delayed=r['immediate_vs_delayed'],matching_reference_response_contrasts=r['matching_reference_response_contrasts'],control_drift=r['control_drift'],waiting_source_checks=r['input_audit']['waiting_sources'])
 write('docs/'+stem+'_assessment.json',compact)
 text='# '+('Immediate' if phase=='immediate' else 'Delayed')+' response after initial surrounding-amount matching\n\n'
 text+='All prespecified numerical gates pass, including three timesteps, independent DOP853/Radau solves, CPU/GPU agreement and a separate reconstruction from saved arrays. '+str(r['counts']['independent_histories'])+' existing developmental histories supply the inputs; interventions and numerical repeats add zero histories.\n\n'
 text+='The intervention at t=510 matches the recipient’s activator/inhibitor pair and each species’ total surrounding amount, with physical context held common. Each species’ allocation fractions among the surrounding cells are preserved. '+('The pulse follows immediately.' if phase=='immediate' else 'The pulse follows a fixed60-unit wait, from qualified unpulsed endpoints. Neither target nor amounts are rematched at t=570.')+' Reactions then change amounts freely.\n\n'
 text+='| History | Recipient | Fine response RMS range | Primary screen across all combinations |\n|---|---|---|---|\n'
 for x in ranges:text+=f"| {x['history']} | {x['cell']} | {x['min']:.8f}–{x['max']:.8f} | {'Above' if x['all_above'] else 'Below' if x['all_below'] else 'Mixed'} |\n"
 text+='\nRanges span both original matching references and both shared physical contexts'+(' and both doses' if phase=='delayed' else '')+'. The strict primary screen is E>0.01; below-screen is not equivalence. Decisions agree across all five methods. Physical contrasts and interactions remain descriptive.\n\n'
 text+='| Numerical evidence | Value |\n|---|---|\n'
 for label,value in [('GPU paths',r['counts']['production_paths']),('CPU verification paths',r['counts']['cpu_paths']),('Dual-solver trajectories',r['counts']['reference_trajectories']),('Logical starts',r['counts']['logical_starts']),('Within-label distinct problems',r['counts']['unique_initial_value_problems']),('Physically distinct problems ignoring labels',r['counts']['unique_physical_initial_value_problems']),('Raw CPU/GPU maximum',f"{r['maxima']['cpu_torch']:.3g}"),('Production/reference maximum raw log error',f"{r['maxima']['production_reference']:.3g}"),('Resolved order cases',r['gate_counts']['resolved_orders']),('Full coordinator time (minutes)',f"{elapsed['coordinator_wall_seconds']/60:.2f}")]:text+=f'| {label} | {value} |\n'
 text+='\nUnresolved orders remain unresolved; no acceptance limit is changed. Recovery is measured relative to each background’s own evolving control, with censoring retained.\n\n'
 for k,x in recovery.items():text+=f"{k.replace('_',' ')}: {x['recovered']}/{x['logical_pulse_labels']} labels recover; times {x['min']}–{x['max']}.\n\n"
 if phase=='delayed':
  checks=r['input_audit']['waiting_sources'];text+=f"All {len(checks)} delayed sources pass the1e-8 endpoint gate. {sum(x['settled'] for x in checks)}/{len(checks)} pass the predefined descriptive settling screens. Largest preceding24-unit raw log drift: {max(x['late_max_log_drift'] for x in checks):.3g}; largest instantaneous absolute log rate: {max(x['maximum_relative_rhs'] for x in checks):.3g}.\n\n"
 text+='Later amount differences are outcomes of the freely evolving reactions, not evidence of a continuing matched budget. Full paired whole-network and surrounding amount curves are reconstructed in the raw-array review. This tests initial allocation at a specified common initial budget; it does not isolate nearest-neighbor signaling, dose independence at every later time, autonomous identity or new development.\n\n'
 text+='[Protocol](phase_carry_recipient_amount_response.md), [machine-readable assessment]('+stem+'_assessment.json), [separate array review](../'+str(review/'independent_review.json')+'), [archived waveforms](../'+str(review/'waveforms.npz')+').\n'
 Path('docs/'+stem+'_assessment.md').write_text(text)
 assets=[f for folder in (review,launch) for f in folder.rglob('*') if f.is_file() and '__pycache__' not in f.parts and f.name!='review-provenance.json']
 assets += [Path('docs/'+stem+'_assessment.md'),Path('docs/'+stem+'_assessment.json'),Path('docs/phase_carry_recipient_amount_response.md'),Path('docs/phase_carry_recipient_amount_response_launch.md'),snap]
 provenance=dict(verified_utc=datetime.now(timezone.utc).isoformat(),passed=True,protocol_sha256=compact['protocol_sha256'],original_summary_sha256=compact['original_summary_sha256'],original_status_sha256=compact['original_status_sha256'],independent_review_sha256=compact['independent_review_sha256'],independent_review_script_sha256=sha(review/'independent_review.py'),source_and_input_pins_unchanged=True,source_pin_count=len(p['source_sha256']),input_pin_count=len(p['input_sha256']),immutable_consumed_artifacts=len(r['inputs_sha256']),ledger_manuscript_and_parent_unchanged=True,unchanged_guards=guards,files_sha256={str(f):sha(f) for f in sorted(set(assets))})
 write(review/'review-provenance.json',provenance)
 print(json.dumps(dict(accepted=phase,recipient_effect_ranges=ranges,recovery=recovery,amount_diagnostic=amount_diagnostic,timing=elapsed),indent=2))
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('phase',choices=('immediate','delayed'));main(parser.parse_args().phase)
