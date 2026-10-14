"""Write human-readable reporting after a complete passing independent audit."""
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).parent
ROOT = Path('output/studies/development-gpu-cleavage-qualification-v2')


def read(path):
    return json.loads(Path(path).read_text())


def main():
    review = read(HERE / 'independent_review.json')
    assert review['passed'] and not review['independent_review_pending']
    protocol = read(ROOT / 'protocol.json')
    assert hashlib.sha256((ROOT / 'protocol.json').read_bytes()).hexdigest() == review['protocol_sha256']
    backend = review['backend_rows']; timesteps = review['timestep_rows']
    maximum = lambda rows, key: max(row['errors'][key] for row in rows)
    metrics = ('chemical_log_max', 'polarity_abs_max', 'relative_axis_max',
        'relative_volume_max', 'relative_transport_max', 'final_phi_abs_max',
        'ancestor_chemical_log_max', 'whole_amount_relative_max')
    names = ('Cellwise log concentration', 'Polarity component', 'Relative tissue axis ratio',
        'Relative cell volume', 'Relative transport-operator norm', 'Final phase field',
        'Ancestor-averaged log concentration', 'Relative whole species amount')
    reports = [read(ROOT / job['key'] / (arm + '-result.json')) for job in protocol['jobs'] for arm in ('gpu', 'cpu')]
    archive_path = HERE.resolve().relative_to(Path.cwd().resolve()).as_posix()
    events = [event for report in reports for event in report['events']]
    audits = [report['audit'] for report in reports]
    summary = dict(passed=True,scope=review['scope'],new_scientific_histories=0,
        backend_cases=16,timestep_pairs=8,protocol_sha256=review['protocol_sha256'],
        independent_review=str(HERE / 'independent_review.json'),
        independent_review_sha256=hashlib.sha256((HERE / 'independent_review.json').read_bytes()).hexdigest(),
        backend_maxima={key:maximum(backend,key) for key in metrics},
        timestep_maxima={key:maximum(timesteps,key) for key in metrics},
        native_prefix_phi_max=max(row['native_prefix_phi_error'] for row in backend),
        native_prefix_chemical_max=max(row['native_prefix_chemical_error'] for row in backend),
        volume_error_max=max(a['max_volume_error'] for a in audits),
        radius_min=min(a['min_radius'] for a in audits),
        event_amount_jump_max=max(e['amount_relative_jump'] for e in events),
        event_volume_jump_max=max(e['volume_relative_jump'] for e in events),
        represented_event_amount_jump_max=max(r['represented_gpu_event_amount_jump'] for r in backend),
        event_target_jump_max=max(e['target_relative_jump'] for e in events),
        represented_occupancy_jump_max=max(e['effective_occupancy_absolute_jump'] for e in events),
        dilution_amount_error_max=max(a['dilution_amount_error'] for a in audits),
        total_seconds=read(ROOT / 'timing.json')['seconds'],source_pins=review['source_pins'],input_pins=review['input_pins'],
        previous_evidence_pins=review['previous_evidence_pins'],consumed_pins=len(review['consumed_sha256']))
    table = '\n'.join(f'| {name} | {maximum(backend,key):.3g} / {protocol["limits"][key]:.3g} | {maximum(timesteps,key):.3g} / {protocol["time_limits"][key]:.3g} |'
        for name,key in zip(names,metrics))
    text = f'''# Development GPU single-cleavage qualification: completed assessment

The independent saved-array review accepts all **16 CPU/GPU cases and eight timestep pairs** under the unchanged prospective limits. The isolated developmental extension is qualified for the prescribed single-cleavage design at grid 72 through t=3. This is **zero new scientific developmental histories**. Naturally scheduled repeated cleavage and developmental spatial/time refinement remain unqualified.

[Prospective protocol](development_gpu_qualification.md), [launch record](development_gpu_qualification_launch.md), [literal independent review](../{archive_path}/independent_review.json), [review code](../{archive_path}/review_saved_arrays.py), [machine-readable assessment](development_gpu_qualification_assessment.json).

## Method and controls

Resident PyTorch arrays and chemical/transport matrix operations are retained, with a separate custom CUDA contractile-ring mechanics kernel. The accepted mature sources are unchanged. The new backend preserves phase carry through pressure correction and daughter partition, conserves chemical amounts and volume targets at abscission, resizes cell buffers, and retains polarity, lineage, cycle/partition RNG states and numerical corrections on restart. Missing developmental schema/correction fields are rejected.

An analytic radius-0.4 diffuse cell starts with uniform activator 1.03 and inhibitor 0.98. Two prescribed cleavage directions, four mechanical branches and dt=0.00375/0.001875 produce sixteen synthetic cases. All start one progressive division at t=0, are capped at two cells and are followed through t=3. The branches are coupling off, tension/adhesion only, directional polarity tension only and the full combination; chemistry remains active throughout. The oblique cases have a specified initial polarity vector. These axes and polarity vectors exercise numerical mechanics and are not spontaneous-axis evidence.

The independent NumPy reference implements the developmental equations without calling the GPU solver. A separate 0.06-unit historical native prefix checks the numerical extension against prior arithmetic. The mother stays a well-mixed chemical compartment until abscission; gradients between connected lobes are omitted. Beginning-geometry splitting remains first order. Pressure correction and occupancy splitting preserve the represented field, extending historical float32 cleavage rather than claiming identical old-method trajectories.

## Numerical results

Each entry gives maximum observed discrepancy / declared limit. Cellwise comparisons use identical IDs; ancestor-averaged chemistry and whole amounts cover all observations. **No unmatched observations occurred** in the sixteen backend or eight timestep comparisons.

| Metric | NumPy/GPU | Coarse/fine GPU |
|---|---|---|
{table}

The largest native-prefix phase discrepancy is {summary['native_prefix_phi_max']:.3g} against 2e-5, and its largest chemical log discrepancy is {summary['native_prefix_chemical_max']:.3g} against 1e-5. All sixteen independent eight-step chemical/dilution checks pass the 1e-11 limit. The existing mature branch and the new no-cleavage branch agree bit-for-bit for eight steps, including phase carry, chemical correction and geometry.

Every case completes one resolved division: coarse abscission at t=0.90375 and fine at t=0.901875. CPU/GPU event times agree; the coarse/fine shift 0.001875 is below 0.0075. Every daughter has one connected component at the event boundary and saved nondividing observations. The pinching mother is excluded from that connectivity requirement.

Maximum visible-volume jump at abscission is {summary['event_volume_jump_max']:.3g} against 1e-6; recorded chemical-amount jump is {summary['event_amount_jump_max']:.3g} against 2e-14. Independent long-double reconstruction of the GPU represented chemical pairs gives maximum amount jump {summary['represented_event_amount_jump_max']:.3g}. Target-volume jump is {summary['event_target_jump_max']:.3g}; represented occupancy discrepancy is {summary['represented_occupancy_jump_max']:.3g} against 1e-14. Maximum every-step dilution closure error is {summary['dilution_amount_error_max']:.3g}. Maximum volume-target error is {100*summary['volume_error_max']:.3f}% against 5%, minimum equivalent radius is {summary['radius_min']:.3f} grid intervals against four, and no mechanical clipping occurs. Every saved boundary screen passes. All saved transport operators preserve constants and volume-weighted amounts, with symmetric volume-weighted conductances and nonnegative exchange coefficients.

## Restart and preservation evidence

All sixteen paths pass eight exact restored/uninterrupted steps during the live furrow and eight immediately after abscission. A separate full replay from t=0.45 through t=3 reproduces the final tensors, corrections, geometry, cell count, IDs, lineage, events, clocks and RNG states bit-for-bit in **all sixteen cases**. [Replay record](../{archive_path}/restart-replay.json). The initial replay script compared the decimal clock to exact 0.45; correcting that audit-only comparison to the existing 1e-12 clock tolerance resolved the script assertion. The saved-array review initially required every due time to be finite, overlooking the model's deliberate +infinity marker for an already dividing mother. Its corrected check permits that marker only for active mothers. Both first attempts and scripts are archived; solver code, data and acceptance thresholds were unchanged.

Two additional tests on the GPU pass: an ellipsoid divides along its long axis without changing fields at division initiation; and naturally scheduled development reaches four cells with an exact restart from the three-cell stage while another mother is dividing. These are software checks, not repeated-development accuracy or new histories. Twenty-eight unit tests and five explicitly enabled existing CUDA backend tests also pass.

All {summary['source_pins']} source pins, {summary['input_pins']} input pins and {summary['previous_evidence_pins']} previous scientific/ledger/manuscript pins verify. The review records {summary['consumed_pins']} consumed hashes and the runtime binary hash. Original production reports retain their independent-review-pending flag; this separate review supplies the acceptance decision. The preliminary unsuffixed prepared directory was never run, and its earlier source version is preserved. The qualification took {summary['total_seconds']/60:.2f} minutes with four parallel CPU references; this is a qualification runtime, not a full-development speed benchmark.

## Scope and next gate

Immediate GPU event amounts, targets and occupancy volumes are reconstructed from saved arrays. The pointwise represented-occupancy jump, intermediate component counts, every-step extrema and CPU event jumps are recorded runtime evidence rather than raw-field reconstruction of every step. Final CPU/GPU and immediate GPU daughter fields are independently relabeled. Restart replay uses the candidate solver; the independent physical comparison is the separate NumPy trajectory.

Next qualify **naturally scheduled repeated cleavage from one cell toward sixteen**, retaining live-event restarts, matched controls, temporal refinement, cell/interface resolution and boundary screens. Then declare unused seeds and a separate protocol for three new independent developmental histories. This single-cleavage pass does not establish spatial convergence, a tissue continuum, chemical pattern initiation, inherited identity or biological lineage commitment. Existing scientific claims and manuscript remain unchanged.
'''
    Path('docs/development_gpu_qualification_assessment.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    Path('docs/development_gpu_qualification_assessment.md').write_text(text)
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    main()
