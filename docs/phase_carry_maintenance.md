# Carry-corrected moving maintenance

The [completed polarity/conductance controls](polarity_conductance_controls_assessment.md) qualify initiation on three existing mature geometries and identify a tested contribution of evolving conductances to its suppression. The next question is whether **already developed chemical differences survive the same moving mechanics after retaining sub-float32 phase updates**. Earlier maintenance results used the original no-carry method; they remain original-method evidence until this matched test passes its own gates.

## Experimental question and inputs

Use the original mature 16-cell states at physical time 150 from histories 7, 8 and 9. Within each history, the developed and near-uniform chemical preparations are the exact supplied inputs from the earlier polarity studies; do not generate new noise, patterns or preconditioning. Keep cell geometry, polarity, volumes, IDs, lineage and random states matched. Coarse/fine sources differ only in numerical clock/timestep bookkeeping.

Developed initial contrasts, measured as the across-cell standard deviation of natural-log activator, are approximately 1.18186, 1.20362 and 1.28962 for histories 7, 8 and 9. They are chemical states, rather than assigned identities or autonomous fates. Their prior formation and the inherited mature geometry used the old method; the new experiment qualifies their continuation, not their developmental origin.

| Chemical preparation | chi=0 | chi=0.35 | Numerical levels |
|---|---|---|---|
| Developed pattern | New moving maintenance runs | New moving maintenance runs | dt=0.00375 and 0.001875 |
| Near-uniform | Qualified carry initiation reuse | Qualified carry initiation reuse | Same two timesteps |

There are **twelve new maintenance paths and twelve reused initiation paths nested within three histories**. Histories are the replication units. Contrasts, timesteps, cells and endpoint preparations do not add replicas. History 9 is scheduled first because it is the one history where the qualified initiating control forms and positive directional tension suppresses formation.

## Model and compute path

Keep beta=2, D_a=0.02, D_b=0.55 (ratio 27.5), activity-tension/adhesion coefficients 0.25/0.35, grid 72 cubed, extent 2.24 and interface parameter 0.085. chi changes only the directional polarity tension action; polarity dynamics continue even at chi=0. Activator-dependent mechanics remain active. Cell shapes and contacts evolve in every new path.

Unlike the preceding intervention, **chemical conductances evolve with geometry** in this test. The conservative chemical operator is

$$
\Delta(t)=-M(t)^{-1}K(t),
\qquad M(t)=\operatorname{diag}(V_i(t)),
\qquad K(t)=\operatorname{diag}(G(t)\mathbf 1)-G(t).
$$

Mechanical concentration conversion by V_old/V_new remains active, with no artificial volume amount source. Chemistry and mechanics co-evolve; there is no fixed reservoir or fixed transport override.

Use resident PyTorch arrays/matrix operations and the unchanged custom CUDA carry mechanics and geometry/polarity kernels. Visible fields and contact accumulations remain float32, chemical/geometry arrays use their existing precisions, and phase-update residuals are float64. Residuals begin at zero when creating a new mature-state continuation, as in the qualified initiation study; checkpoint resumption must preserve them exactly. No previous scientific source or protocol is edited.

## Qualification before launch

Verify the completed initiation/contact study, its final reviewed evidence and all inherited source/input hashes. The reused near-uniform paths must match history, chemical preparation, point, timestep, physical configuration and actual endpoint; both their pilot and full raw comparisons must pass.

Every new developed context receives:

- Sixteen steps reproducing the original GPU arithmetic exactly in the numerical baseline, including an exact first carry step from zero residual.
- Eight independent NumPy evaluations of chemistry and mechanical amount accounting, using the actual GPU transport matrix and measured old/new volumes. Chemical log error must be at most 1e-11 and accounting error at most 2e-14.
- A carry checkpoint and four exactly matching resumed steps, including chemistry, geometry, polarity, residuals, rounding counters and transport.
- A 0.6-unit native CPU versus actual carry GPU coupled comparison at the exact developed preparation, geometry, polarity, chi and timestep. Retain the original strict native/GPU discrepancy criteria, including chemical/polarity errors 1e-5, relative transport 1e-5 and final field error 2e-5.

Also compare mechanics against the independent NumPy/SciPy float64 reference at the prespecified **developed history-9 chi=0.35** start, for 0.15 units and dt=0.00375/0.001875/0.0009375. Concentrations and polarity are held fixed in this component check, while geometry and mechanics evolve. Keep the original field, volume, clipping and reference convergence criteria. These three nested component checks supplement the twelve actual coupled context gates; they do not validate all developed states by themselves.

Freeze protocol, sources, inputs, preparation evidence and decision rules after these checks pass. A failed preparation does not authorize scientific continuation.

## Moving experiment and acceptance

Run each new coarse/fine pair to a 60-unit pilot. Continue that pair to 240 elapsed units only if it passes the raw pilot refinement gate. A failed pair remains unresolved at the pilot; independent pairs may still run. Save observations every 0.15 units and atomic checkpoints every three units. End at physical time 390. Resume checkpointed chemistry, polarity and residuals; never replace an evolving state with its original preparation.

Use the original full-window limits: maximum raw chemical log discrepancy 0.01, polarity absolute error 0.01, relative axis error 0.01, relative volume error 0.005, relative transport error 0.01, frozen growth error 0.001 and growth/onset crossing-time error 0.3. Compare physical times directly, with no temporal alignment or relaxed thresholds. Both levels must agree on the late contrast decision. For developed starts, also require the same continuous-retention decision and first contrast-loss time agreement within 0.3 units, including the same presence/absence of loss.

Require maximum target-volume error below 5%, equivalent radius at least four grid spacings, zero clipping, finite positive chemistry, sampled boundary occupancy below 0.01 and dilution amount-accounting error at most 2e-14. A scientific loss of contrast can pass numerical quality; it is not itself a solver failure.

## Outcomes and reporting rules

Define chemical contrast

$$
S(t)=\operatorname{SD}_{i}\bigl(\log a_i(t)\bigr).
$$

The primary maintenance outcome is S(t)>0.1 at every saved observation throughout elapsed 216–240. Report the late-window minimum and final contrast. Separately report the minimum over the entire window, whether all saved observations remain above 0.1, and the first downward crossing of 0.1. A late recovery after transient loss is not continuous retention. Finite observation spacing limits detection of unsampled transients.

For each qualified history, report both maintenance branches and the reused initiation branches separately:

| Developed-start outcomes | Interpretation |
|---|---|
| Both contrasts maintain | Developed chemistry persists with or without directional polarity tension in this moving context. Initiation may still differ. |
| chi=0 maintains; chi=0.35 loses | Positive directional tension reduces maintenance under the tested continuation. |
| Neither maintains | The tested moving system does not retain either supplied developed pattern over the declared late window. Do not label this specifically a polarity effect. |
| Only chi=0.35 maintains | Report the nonmonotonic maintenance result explicitly. |
| Numerical/physical/endpoint gate fails | That history is unresolved; keep outcomes and failed checks without accepting the scientific conclusion. |

Run the original dual-solver frozen chemical endpoint assay on each new path's actual final operator, with the same stationarity, Jacobian stability and local-return criteria. Reuse qualified initiation endpoints read-only. Twenty-four endpoint graphs and 120 chemical preparations remain nested within the three histories. Distinguish stable patterned chemistry, stable uniform chemistry and local uniform/patterned coexistence. Frozen chemical stability is not stability of the complete moving system.

The test does not establish biological identity, autonomy, inheritance, fresh carry development, large emergent shape symmetry breaking, general transport closure, spatial convergence or a broad backend qualification. Broader exchange, neighbor/reservoir and parameter claims still require their own carry revalidation. Preserve original ledger decisions and old numerical failures.

## Execution

Runner: [phase_carry_maintenance.py](../embryo/phase_carry_maintenance.py). Targeted verification: [tests](../tests/test_phase_carry_maintenance.py). Output directory: `outputs/phase-carry-maintenance/`.

Twelve new GPU paths run sequentially on the qualified GTX 1080 Ti. Two CPU workers overlap independent frozen endpoint solves. This follows the measured single-GPU scheduling policy; concurrent GPU trajectories would compete for memory bandwidth. The protocol, status, run log, launch metadata and context/reference evidence record execution. Preparation and launch status are reported separately from scientific acceptance.
