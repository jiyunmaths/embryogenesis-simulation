# Implementation plan

## Objective

Determine how far a small, interpretable set of mathematical equations can explain the emergence of organized development from a single zygote. Explicit activator–inhibitor feedback is a core requirement: activation, induced inhibition, and spatial communication must be represented dynamically and coupled to cell identity and changing 3D geometry. Start with two identities and a generic model before choosing organism-specific hypotheses.

The [scientific objective and criteria in the README](../README.md#scientific-objective) govern the remaining work. Keep identity differentiation, signaling symmetry breaking, and geometry/shape symmetry breaking as separate measured outcomes. Treat complex structures as hypotheses to explain rather than shapes to prescribe.

## Next scientific priority and decision rules

Gierer–Meinhardt signaling with conservative volume-weighted contact transport and apical–basal mechanical feedback are now implemented. Finite-graph spectra and cleavage mode transfer are checked before interpreting patterns. The next priority is to establish developmental timescale, parameter, and resolution regimes where these equations yield persistent organization. Activator/inhibitor variables must be distinguished from A/B identity labels.

For each model extension, record the proposed mechanism, quantities already imposed by the model, an ablation or perturbation that could challenge the explanation, and quantitative acceptance criteria. Use unbiased initial fluctuations for spontaneous-symmetry-breaking experiments. Keep externally imposed signals as labeled controls. Analyze stability and spatial modes before claiming a Turing mechanism; report regimes without organization as well as successful patterns.

The first [signaling-versus-cleavage timescale screen](timescales.md) is complete for seed 7: the time-15 baseline develops large signal differences by time 60, and doubling the cycle interval mainly delays reaching the larger graph. This resolves the initial short-run observation for one parameter set; contact-threshold sensitivity, mechanical/grid refinement, and independent-seed replication remain the next acceptance work.

### Domain-size verification for the conservative model

The [fixed-spacing domain study](domain_size.md) compares half-widths 1.60, 1.92, and 2.24 on grids 40, 48, and 56. The larger dashboard preset is available in `configs/large_domain.json`. The matched zygote screen through $t=30$ is complete: the largest-domain boundary maximum is 1.90e-11, and the two largest domains differ by at most 1.41e-7 in relative axis ratio. Both declared screens pass, while the minimum-cell-radius resolution screen remains below four spacings. The time-180 continuation is also complete: the original box fails the boundary screen from 131.4 onward, while the larger pair passes and differs by at most 0.000673% in axis ratio. This is a single-seed domain result; the small-cell-resolution shortfall remains. Adding empty space does not replace cell-resolution refinement, and the historical normalized-model boundary failure is not treated as a result for the new model.

### Coupled resolution screen and projection audit (complete; independent confirmation passed)

The [manufactured 16-cell protocol](coupled_resolution.md) independently varies spatial grid (56, 72, 88 at fixed domain and interface width) and time step (0.015, 0.0075, 0.00375). Frozen-geometry quadrature and spectra agree closely under refinement, and finer grids clear the four-spacing radius criterion. Five short full-coupling cases are complete. Eleven of twelve checks pass; the 1.43% finest-pair spatial field discrepancy fails the 1% criterion. The projection audit is complete: the original linear diagnostic already differs by 1.49% at the identical initial state, whereas higher-order final comparisons range from 0.24% to 0.32%. The [independent 112³ confirmation](refinement_confirmation.md) passed all 11 criteria recorded before evolution, with finest-pair field differences of 0.162–0.207%. This supports the revised diagnostic for this short manufactured-state screen. Keep the archived failure unchanged; do not substitute an uncalibrated diagnostic or loosen the threshold. Cleavage and long-time pattern convergence remain separate tests.

### Contact-cutoff sensitivity (frozen and dynamic screens complete)

The [seven-state frozen screen](contact_sensitivity.md) retains graph connectivity and unstable-mode counts throughout the tested cutoff range, but fails the declared all-spectral-growth-rate tolerance at developmental times 30, 60, 90, and 120. The fastest spatial growth changes by at most 0.0016 in the local cutoff neighborhood. The [matched dynamic branches](cutoff_dynamics.md) completed t=30 to t=90 with all five sensitivity checks passing while polarity cutoff was fixed at 0.02. All sampled identity labels match; maximum relative axis-ratio difference is 1.19e-6. The starting 12 A / 4 B labels make this an identity-retention experiment; the known coarse-radius failure remains a separate quality decision. Preserve the failed quantitative screen.

### Controlled cleavage refinement (complete; revised measurement validation passed)

The [small-cell cleavage screen](cleavage_resolution.md) isolates one division with common analytic initialization, fixed interface width, grids 56/72/88, independent time refinement, and axial/oblique cuts. It passed 24/25 checks, including conservation immediately around abscission, daughter connectivity, event timing, and shape. Initial reconstruction of the 56³ field failed its calibration tolerance (0.477% versus 0.25%); a separate [native-phase reconstruction validation](cleavage_measurement.md) now passes all 10 checks under unchanged tolerances. The original failed report remains preserved. This addresses division absent from the earlier manufactured-state refinement; it does not replace repeated-cleavage or long-time developmental convergence.

### Full developmental refinement (running)

At the user’s request, the [full trajectory study](development_refinement.md) now advances five fresh analytic zygotes through t=90. Spatial grids 56/72/88 use dt=0.0075; two additional 72³ runs use dt=0.015 and 0.00375. All physical parameters and continuous target volumes remain fixed. Compare shape, signal distributions/contrast, fate fractions, cell counts, and exact last-division times; retain lineage and spatial records without assuming numeric IDs are homologous. Candidate 72³/88³ quality and coarse-reference limitations are reported separately. The earlier cleavage calibration issue is now resolved in a separate native-phase reconstruction validation; the original failed result and this study’s frozen protocol remain unchanged.

### Geometric transport closure (benchmark complete; general-geometry checks fail)

The [known-geometry validation](geometric_transport.md) passes flat calibration, numerical integration, conservation, and orthogonal flux checks, but fails curved-area, positive-gap, and nonorthogonal linear-flux tolerances. The [face-normal-aware prototype](skew_flux.md) is complete: seven of eight checks pass, including corrected linear fluxes and second-order refinement, but positivity fails even under exact discrete evolution. The [positive wider-stencil correction](positive_skew_flux.md) now passes all 11 checks on periodic affine meshes. The [reflecting-wall benchmark](skew_boundary.md) passes 8/9 checks but fails the convergence-order criterion at intermediate shears. The [derived wall correction](skew_boundary_correction.md) now passes all nine checks, including independent wall families and 64³ confirmations near second order. Next validate nonuniform capacities and graded skew geometry before live interface extraction and coupling. Preserve the unchanged current solver while assessing its refinement; do not infer general continuum consistency from conservation or cutoff insensitivity.

### Causal feedback screen (first fixed-geometry stage complete)

At the user’s request, the [paired causal screen](causal_signaling.md) tests seven interventions on the resolved 72³ post-cleavage graph using twenty matched chemical perturbation seeds. Full-loop persistent signals appear in 20/20, while no self-activation, no transport, and equal diffusion suppress contrast. Both fate labels nevertheless appear in all those controls; the separate bistable switch can preserve/amplify transient activity differences. Signal-to-fate ablation prevents commitment. Inhibition knockouts become locally unstable and hit the predeclared ceiling, not an organized-pattern success. Historical fine-time checks failed continuous-fate tolerances in two controls. The [joint signal/fate validation](joint_fate.md) now passes against tightened DOP853 references for all five tested arms and twenty seeds, without relaxing tolerances. Four [moving-geometry controls](moving_causal.md) are running through t=78: full, no mechanical feedback, no self-activation, and no signal-to-fate coupling. Their shape comparison remains pending. This priority shift does not clear the remaining geometric-transport limitations.

### Live conservative coupling (implemented; general geometry convergence pending)

[Conservative live transport](live_transport.md) now includes contact-area calibration, measured capacities, dilution, conservative cleavage, volume-weighted spectra, dashboard preflight, and exact prescribed-geometry replay. Flat-contact refinement recovers the fixed-domain continuum unstable mode with approximately second-order spectral convergence. Next: vary voxel spacing at fixed cell geometry/interface width; independently reduce interface width; test contact cutoff, gaps, curvature, and nonorthogonal faces; refine coupled time steps; rerun timescale/shape controls and seed ensembles with physical diffusivities fixed. Historical shape outcomes do not validate this new model.

### Discrete-to-continuum bridge (fixed-domain signaling and prescribed motion verified)

The next scientific direction is to connect the small-cell model to a continuum description while keeping signaling, identity, and shape as separate explanatory targets. A [fixed-domain conservative transport benchmark](continuum_bridge.md) now verifies volume-weighted conservation, approximately second-order diffusion convergence, and early GM growth against discrete and continuum mode predictions on successively refined 3D meshes. It does not establish a continuum limit of the actual embryo contact graph.

A second [irregular-domain benchmark](irregular_bridge.md) now verifies smooth diffusion convergence on a nonuniform orthogonal mesh in a 3D L-shaped prism. It includes true face conductances, unequal capacities, correct no-flux walls, and an independent pulse that crosses between the arms and equilibrates. All ten default criteria pass; the finest graded refinement has observed order 1.992. This is an irregular-boundary test, not validation of arbitrary unstructured cells or an emergent embryo shape.

The [weighted-spectrum and early-signaling benchmark](irregular_signaling.md) is also complete. Full discrete spectra give 9, 7, 4, and 4 unstable spatial modes as the mesh is refined; selected growing and decaying nonlinear probes agree with their linear predictions. All nine default criteria pass, including a 1.11% last-refinement change among the first twelve positive eigenvalues. This does not certify an exact continuum spectrum or nonlinear pattern persistence.

The [nonlinear persistence benchmark](nonlinear_bridge.md) is now complete for one prescribed continuous perturbation. Conservative spatial comparison reduces the final field discrepancy from 5.43% to 1.26%; independent time refinement gives a 0.095% maximum sampled difference for the smaller step pair. All eleven criteria pass, including late contrast, less than 1% late field change, and decay under equal diffusivities. This is finite-window evidence, not an ensemble or attractor-stability test.

The [moving-domain benchmark](moving_domain.md) now passes eleven checks, including inverse-volume dilution, amount drift below $2.2\times10^{-14}$, and approximately second-order space/time convergence for three prescribed affine motions. This holds connectivity fixed and has no chemical reactions.

Next: conservative remapping during compartment splits/refinement and coarsening, followed by a hybrid field/cell coupling. Nonlinear robustness across physical perturbations, parameters, and disturbances remains pending. Choose bulk, membrane-limited, or extracellular communication explicitly before interpreting contact weights as molecular conductances. Preserve fate memory and mechanics during this transition. Contact-cutoff sensitivity and the other validation requirements below remain necessary.

## Milestone 1 — mechanics and cleavage (implemented; coarse numerical checks)

- Diffuse deformable cells with volume constraints, exclusion, and interface attraction.
- One-cell, pair-contact, and cleavage checks.
- Shape-oriented spindles, progressive contractile-ring surrogate, conservative abscission, and explicit lineage tracking.
- Reproducible output and volume/boundary diagnostics.

Next acceptance work: establish a resolution and time-step range with acceptable individual volume error, cell connectivity, contact geometry, and relaxation speed. Compare pair-contact angles as attraction changes. Test rotated configurations for grid bias. Use grids fine enough to resolve the smallest cells with several interface-independent interior voxels.

## Milestone 2 — activator–inhibitor signaling and two identities (graph model implemented; robustness pending)

- Gierer–Meinhardt signaling, constant-preserving normalized graph exchange, a downstream fate switch, discrete-mode stability, and cleavage spectral transfer.
- Legacy independent fate noise and geometry bias are off by default.
- Paired controls with independent random streams.
- Thresholded fate counts and continuous state recording.

Required next validation (the graph kinetics, transport, and linear analysis below now have an implementation):

- Explicit Gierer–Meinhardt activities, production, turnover, and inhibition are implemented. Validate their parameter ranges and timescales relative to cleavage and mechanics.
- Conservative volume-weighted exchange, geometric dilution, and amount-preserving cleavage are implemented. Test contact-threshold sensitivity, geometric closure accuracy, and refinement-driven mode changes; normalized exchange remains a historical control.
- The homogeneous equilibrium, Jacobian, discrete mode growth rates, and cleavage mode transfer are implemented and checked. Extend analysis to evolving-graph transient amplification and nonlinear regimes; do not assume a fixed physical wavelength or a single pole.
- Activator-biased fate dynamics are implemented. Measure activity distributions and dwell times, withdraw partition noise after differentiation, and perturb or reposition cells to assess persistence and reversibility.
- Test removal of self-activation, inhibitor production/action, and signal transport separately. Map outcomes over parameters and seeds before expanding the identity count or asserting robustness.

Acceptance: a reproducible account of when the loop amplifies or suppresses spatial fluctuations, whether the resulting signals generate persistent fate differences, and which mechanisms are necessary. Two identities need not imply two signaling domains or a changed embryo shape.

## Milestone 3 — two-way geometry feedback (prototype implemented)

- Contact/exposure influences regulatory state.
- Activity changes surface tension and interface attraction.
- Compare mechanics-only, geometry-to-fate only, and full-feedback controls.

Next acceptance work: couple the validated activator–inhibitor module to mechanical properties and update its transport/sensing as geometry evolves. Add contact enrichment, fate-versus-exposure plots, and signal spatial correlations; isolate adhesion, tension, and each feedback direction. Run at least 20 seeds per screened parameter set before drawing robustness conclusions. Test fixed preassigned activities as a distinct sorting experiment, clearly separated from emergent differentiation.

## Milestone 4 — 32–64 cells and reliable geometry (pending)

- Benchmark memory and time; avoid prematurely increasing cell count on a fixed coarse grid.
- Profile local field storage, neighbor culling, and accelerated kernels.
- Refine the grid at fixed physical interface width, then study the diffuse-interface approximation separately.
- Add explicit surface reconstruction, area, cell shape axes, and connected-component checks.
- Maintain free embryo boundaries; model a deformable external envelope only as a separate experiment.

## Milestone 5 — persistent axis formation (controlled pilot implemented; validation pending)

The [coupled shape-persistence pilot](shape_persistence.md) measures identifiable long-axis memory after cleavage and compares matched mature-state interventions with a separate development-without-feedback control. Persistence, feedback-specific excess, and numerical quality are separate decisions. It uses the existing contact-graph model; those historical results used normalized exchange; new runs use conservative transport. The completed seed-7 screen finds persistent ratios 1.323–1.331 in all four controls, no positive feedback-specific excess, and axes within 1.37° of the first-cleavage direction. All runs narrowly fail the four-grid-spacing radius criterion. Prioritize mechanical refinement and cleavage-memory controls before interpreting or strengthening the coupling; robustness remains unestablished.

The next diagnostic pair is implemented: [unchanged-parameter continuation to $t=180$](shape_persistence.md#long-time-continuation), and [external chemical-patch forcing followed by release](signal_patch.md), including uniform controls and a bounded fate–tension-strength comparison. The mechanics kernel optimization is checked against the previous trajectory, separately from the scientific response. Both protocols are now complete. The patch/uniform response grows with the fate–tension coefficient, but the patched aggregate remains near ratio 1.329 at the end of forcing; the feedback-enabled branches also fail the boundary screen. These tests do not change the spontaneous-axis acceptance criteria. The completed $t=180$ continuation still lacks feedback-specific elongation and exposes boundary-screen failure in the full run from $t=124.2$. Increase domain extent at fixed grid spacing and establish spatial/time convergence before interpreting later shape changes; adding cells at the same grid resolution is not the next validation step.

- Apical–basal polarity magnitude/direction, exposed-cortex cues, neighbor alignment, and directional cortical tension are implemented; validate their timescales and influence on tissue-scale organization.
- Shape-aligned division and an isotropic control are implemented; next compare them with polarity-aligned division and a calibrated tensile-stress-based rule.
- Polarity-dependent cortical tension is implemented with its spatial-gradient contribution; test this feedback against no-polarity controls and refine its mechanical interpretation.
- Test whether the validated activator–inhibitor module generates a persistent pole or axis on the changing tissue, separately from externally imposed gradients.
- Derive/check instability conditions and domain-size effects before expecting a single pole.

Acceptance: persistent shape anisotropy beyond division-only controls; identify onset, duration, and axis memory; isotropic axis directions across seeds in an unbiased environment; reproduce after grid rotation and refinement. Check whether the first cleavage axis seeds later asymmetry.

## Milestone 6 — biological and geometric extensions (pending)

Choose only after the prior mechanisms are understood: a resolved extracellular medium beyond a contact-graph signaling model, a lumen with pressure/transport, growth after cleavage, asymmetric cleavage, or localized contraction and folding. Select each extension to test a specific limit of the equations' explanatory power. Two identities alone do not imply any one of these processes.

## Experiment record

Every run records complete configuration, seed, lineage, metrics, visualization trajectory, and final checkpoint. A comparison must name what is held fixed, what changes, the numerical resolution, and the number of independent seeds. With geometry-dependent cytokinesis, using the same seed does not guarantee identical division directions or completion times across feedback controls. Attractive morphology alone is not a success criterion.
