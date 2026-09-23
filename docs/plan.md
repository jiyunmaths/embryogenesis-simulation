# Implementation plan

## Objective

Determine how far a small, interpretable set of mathematical equations can explain the emergence of organized development from a single zygote. Explicit activator–inhibitor feedback is a core requirement: activation, induced inhibition, and spatial communication must be represented dynamically and coupled to cell identity and changing 3D geometry. Start with two identities and a generic model before choosing organism-specific hypotheses.

The [scientific objective and criteria in the README](../README.md#scientific-objective) govern the remaining work. Keep identity differentiation, signaling symmetry breaking, and geometry/shape symmetry breaking as separate measured outcomes. Treat complex structures as hypotheses to explain rather than shapes to prescribe.

## Next scientific priority and decision rules

Gierer–Meinhardt signaling on the normalized contact graph and apical–basal mechanical feedback are now implemented. Finite-graph spectra and cleavage mode transfer are checked before interpreting patterns. The next priority is to establish developmental timescale, parameter, and resolution regimes where these equations yield persistent organization. Activator/inhibitor variables must be distinguished from A/B identity labels.

For each model extension, record the proposed mechanism, quantities already imposed by the model, an ablation or perturbation that could challenge the explanation, and quantitative acceptance criteria. Use unbiased initial fluctuations for spontaneous-symmetry-breaking experiments. Keep externally imposed signals as labeled controls. Analyze stability and spatial modes before claiming a Turing mechanism; report regimes without organization as well as successful patterns.

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
- Normalized contact-graph exchange and activity inheritance are implemented. Test contact-threshold sensitivity and refinement-driven mode changes; assess when volume-capacity-based molecular transport is needed instead of normalized activity exchange.
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

## Milestone 5 — persistent axis formation (pending)

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
