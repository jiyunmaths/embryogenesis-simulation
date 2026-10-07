# Embryogenesis in 3D

## Scientific objective

**To what extent can mathematical equations explain the emergence of organized living structure from an initially simple zygote?** This project investigates how small fluctuations in an approximately uniform initial state can develop into reproducible cell identities, spatial signals, and complex, changing three-dimensional shapes.

The central mechanism to investigate is an **explicit activator–inhibitor feedback loop**: an activator promotes its own production and the production of an inhibitor; the inhibitor suppresses activation. We study how kinetics, interaction ranges, and coupling to cell mechanics determine whether this feedback amplifies fluctuations into order, suppresses them, or produces unstable patterns. Activator and inhibitor are continuous chemical activities, not cell identities.

We will evaluate three linked outcomes separately:

| Outcome | Scientific question | Evidence to measure |
|---|---|---|
| Cell identity differentiation | Can initially similar cells acquire distinct, persistent behavior through local regulatory dynamics and signaling? | Continuous attributes, persistence, responses to perturbations/exchange, context dependence, and eventual inheritance tests |
| Signaling symmetry breaking | Can nearly uniform signaling develop persistent spatial domains, poles, or axes without a prescribed directional cue? | Pattern onset, spatial correlations, wavelength, domain number, and orientation across independent runs |
| Geometry and shape symmetry breaking | Can signaling and cell behavior produce sustained changes in tissue shape while the embryo divides and rearranges? | Shape anisotropy, axis persistence, cell organization, and later folding or cavity formation where the model supports them |

The intended feedback is reciprocal: **signaling influences identity and cell mechanics; cell movement, division, contacts, and shape alter signaling in return.** The embryo boundary must evolve with these interactions.

## How we will assess explanatory power

The aim is to identify the smallest interpretable set of equations that accounts for progressively more developmental organization, and to document where that explanation fails. Two chemical regulators and a small, resolved cell population are the starting point. Whether discrete identities emerge, and how many, are questions for the experiments rather than supplied classifications.

- **Separate emergence from assumptions.** Record which features arise from the dynamics and which are supplied through initial conditions, boundary conditions, cell-cycle rules, or prescribed forces. Begin spontaneous-symmetry-breaking experiments with unbiased fluctuations; label imposed gradients and asymmetries as separate controls.
- **Establish the mechanism.** Analyze steady states and their stability, identify parameter regimes where spatial perturbations grow, and test the necessity of activation, inhibition, signal transport, and mechanical feedback by disabling them individually. An activator–inhibitor loop does not automatically imply a Turing instability or a single developmental axis.
- **Measure robustness and limits.** Compare independent seeds, parameter ranges, perturbations, and numerical resolutions. Report uniform, mixed, fragmented, and failed outcomes as well as organized ones; distinguish reproducible structure from transient noise and grid artifacts.
- **Increase biological specificity only with evidence.** Start with a generic, dimensionless model. Quantitative claims about a particular embryo will require experimental calibration and independent validation. Active signaling, proliferation, and force generation remain explicit biological inputs to the model.

This objective guides subsequent model changes and experiments. Each extension should state the phenomenon it aims to explain, the feedback it introduces, a control that could challenge the proposed explanation, and a measurable success criterion. Visual resemblance alone is insufficient.

## Consolidated evidence and compute paths

**Phase-update precision audit:** the [ledger now separates outcome-level and trajectory-level claims](docs/results_ledger.md#phase-update-precision-audit). All curated biological results used no-carry moving paths or inherited their geometries/states. Previous timestep passes qualify that method; they do not establish insensitivity to the carry correction. Persistence, coexistence and donor-nearer signs remain original-method observations with carry revalidation pending. Onsets, spectral crossing times, conductance declines, response amplitudes and the history-9 **262.1-fold** frozen/moving contrast ratio remain explicitly **old-method measurements**, pending rederivation. Matching final contrast in one history does not validate the other claims. The [completed carry polarity assessment](docs/phase_carry_polarity_assessment.md) qualifies the near-uniform mature-state continuations at chi=0/0.35 in three existing histories; it does not broadly revalidate the original ledger.

**Latest qualified initiation result:** the [completed polarity/conductance controls](docs/polarity_conductance_controls_assessment.md) restore sustained chemical patterns at chi=0.35 in **all three existing histories** when initial chemical conductances are preserved. Matched fixed-conductance chi=0 controls also form, with onset differences below 0.02 model units. The original polarity-specific suppression remains **1 of 3 histories (history 9)**: histories 7 and 8 had neither original branch forming. The intervention removes the demonstrated loss in history 9 and rescues the separate context failures in histories 7/8; it does not turn those latter failures into original polarity-specific effects.

Mechanics, polarity, chemical feedback, dilution and measured volumes continue to evolve. The control preserves symmetric exchange capacities, using $\Delta(t)=-M(t)^{-1}K(0)$ with current volumes. Evolving chemical graphs lose their unstable spatial modes, while fixed-conductance graphs retain two throughout the observed window. This supports **conductance remodeling as a tested route connecting polarity-dependent mechanics to restricted initiation**, without establishing a unique universal cause. The [preceding dilution/contact controls](docs/moving_initiation_controls_assessment.md) also show that removing dilution alone does not restore formation in histories 7/8. These are chemical-state results on inherited mature geometries, not new zygote-to-identity development.

All twelve pilot and twelve full timestep comparisons and all 24 frozen endpoint assays pass. Maximum raw full-window log discrepancy is 0.000410 against the unchanged 0.01 gate; the four new contexts have maximum 5.94e-6. Eight new GPU paths and sixteen reused paths remain nested within three histories. Fixed-conductance endpoints support stable patterns with unstable uniform chemistry; the evolving chi=0 history-9 endpoint supports local uniform/patterned coexistence. The [protocol](docs/polarity_conductance_controls.md) and [launch verification](docs/polarity_conductance_controls_verification.json) remain preserved. Carry maintenance is now qualified at the tested point; fresh-exchange and network-reset behavior also pass at chi=0.35. Delayed responses and dose-controlled context tests remain priorities before expansion.

**Completed carry network-context reset:** the [independent assessment](docs/phase_carry_network_context_assessment.md) accepts all 42 paths (33 new, nine exact reused) across three existing histories through t=450–510. Preserving one recipient while resetting all other cells shifts both selected recipients’ late chemical states in every history. Immediate -10% pulse responses change above the prespecified screen in 1/2 recipients in history 7 and 2/2 in histories 8 and 9; both recipients remain donor-nearer after reset in all three histories, and all pulses recover toward their own moving controls. The largest recipient activator reductions are about 13.6-fold and 17.5-fold relative to sham. This supports maintained, context-sensitive chemical phenotypes with a resilient response tendency associated with prior preparation. Over 99% of the immediate-response difference occurs in the first six units, so delayed-pulse behavior remains untested. The reset changes chemical amounts as well as their distribution and subsequent mechanics; it does not isolate spatial neighborhood, establish network necessity or qualify autonomous/inherited biological identity. [Protocol](docs/phase_carry_network_context.md) and [result figure](outputs/phase-carry-network-context-review/2026-10-06T210547Z-final/network-context.png).

**Active delayed moving-response test:** [New pulses from the actual t=510 sham/reset endpoints](docs/phase_carry_delayed_response.md) test whether response differences persist after chemical-context reorganization. Each background has its own moving control; compare -10% fractional pulses with the same feasible removed amount across backgrounds, normalizing each by its actual log amplitude. The 66 paths are nested within three existing histories and continue through t=570. Implementation tests, dose/state accounting and four full-size native/carry prefixes pass; individual context gates, per-history pilots and full-window paired refinement still govern scientific acceptance. This does not assume stationary states or control the dosage of the original reset. [Launch checks and progress](docs/phase_carry_delayed_response_launch.md).

**Completed carry exchange/response:** the [independent assessment](docs/phase_carry_exchange_response_assessment.md) accepts all 48 moving paths in histories 7–9 at chi=0.35. Both recipients are donor-nearer after fresh, initially amount-preserving exchange and 60 units of coupled evolution; all six comparisons are informative and agree at both timesteps. Concentrations reorganize rather than exactly preserving the transferred reference, while -10% activator response tendencies follow the transferred high/low chemical states. Every pulse recovers toward its own moving control. All 48 context checks and all six pilot/six full-stage reports pass; maximum raw log-concentration discrepancy is 7.16e-5 and pulse-normalized response discrepancy is 0.000673, both below 0.01. The 48 paths and six recipient comparisons are nested within three existing histories. This qualifies the new t=390–510 window, not autonomous or inherited biological identity or older t=150–270 amplitudes. Next test carry-corrected moving neighbor-context/reset dependence. [Original protocol](docs/phase_carry_exchange_response.md).

**Completed carry maintenance:** the [final moving-maintenance assessment](docs/phase_carry_maintenance_assessment.md) qualifies all twelve new developed-start paths across histories 7–9 at chi=0/0.35 and two timesteps. Every path retains chemical contrast at all saved observations through elapsed 240, with no observed loss-and-recovery episode. All six new pilot and six new full comparisons pass; the maximum raw full-window log-concentration discrepancy is 3.51e-5 against the unchanged 0.01 limit. All twelve new frozen endpoint graphs support local uniform/patterned coexistence. Including reused initiation, all twelve pilot and twelve full comparisons and all 24 endpoint assays pass. Preparations, sources and the [launch recovery](docs/phase_carry_maintenance_launch_recovery.md) remain preserved.

Developed chemistry therefore survives moving mechanics that prevents near-uniform initiation at chi=0.35 in all three tested starting contexts. History 9 remains the only original polarity-specific initiation suppression case; histories 7/8 fail to initiate at both original contrasts. These are finite-horizon mature chemical-state results on inherited developmental preparations, rather than fresh carry zygote development, biological identity or full moving-system stability. **Carry-corrected chemical exchange/response and the moving network-context/reset study both pass at this point.** The broader original ledger and other parameter/response claims remain pending their own revalidation.

The [results ledger](docs/results_ledger.md) records 22 curated studies with conclusions, numerical decisions, limitations, source/report hashes, and explicit developmental-history IDs. **The central evidence represents three distinct histories (7, 8, 9), reused across assays.** Cells, paired branches, endpoint graphs, chemical preparations, pulses, and numerical repeats are nested measurements. They do not add developmental replicas, and intervention fractions are not population success probabilities. The [machine-readable ledger](docs/results_ledger.json) also inventories historical result directories without promoting completion status to acceptance.

The current narrative separates **initiation**, **maintenance**, and **network-supported context sensitivity**. Polar mechanics restricts initiation on one tested component geometry; developed patterns survive the feedback switch in three histories; chemical exchange and neighbor interventions show that persistent cell-associated chemistry can depend on the collective network. Both selected immediate-pulse nonrecoveries in history 8 recover when challenged after settling under percentage and amount controls. These are chemical-state results, not autonomous or inherited biological cell identities. General geometric transport closure and full developmental convergence remain unresolved, including their recorded failed checks.

[Backend policy](docs/backends.md) names the active NumPy/SciPy developmental reference, C++/OpenMP mature reference, resident PyTorch/custom-CUDA mature backend, and independent frozen-ODE solvers. [Retired prototype snapshots](archive/backends/README.md) preserve the earlier C polarity, host-transfer CUDA, and pure-PyTorch snapshot experiments. Pinned originals and shared kernels remain at their recorded paths so historical evidence and the active references stay verifiable. Consolidation launches no new scientific simulation.

## Current starting point

The current scientific model, `AttributeSimulation`, develops one cell into a deformable multicellular aggregate. Gierer–Meinhardt activator and inhibitor activities evolve through **conservative volume-weighted contact transport**, reactions, and mechanical dilution. Chemical activity changes tension and adhesion continuously. Apical–basal polarity develops from exposed cortex and neighboring orientations and changes cortical tension directionally. There is no independently integrated fate variable.

**Finite-graph stability precedes pattern interpretation.** Each supported eigenvalue of the conservative Laplacian is checked against the reaction–diffusion Jacobian, with spectral mode transfer recorded at cleavage. A continuous unstable band can contain no supported modes on a small graph. Stability about uniform chemistry describes initiation, not the persistence of a finite-amplitude pattern or stability of the complete moving system. See [live conservative signaling](docs/live_transport.md) and [polarity mechanics](docs/graph_signaling.md#apicalbasal-polarity-and-mechanics).

This is an exploratory model in dimensionless units, not a reconstruction of a particular organism. The simulation evolves each cell's shape on a 3D grid. It does not prescribe an embryo outline or assign daughter identities.

## Definitions and interpretation

**Cell identity means a recognizable, relatively persistent profile of attributes and responses under stated conditions.** It can depend on support from surrounding cells. Autonomy, developmental commitment, and transmission through division are additional properties to test, rather than assumptions hidden in the word “identity.” The current model investigates candidates for such identity through continuous measurements; it has not established biological cell types.

| Term | Definition and use in this project |
|---|---|
| Cell and cell ID | A deformable, well-mixed chemical compartment described by a phase field. Its integer ID tracks that individual through time; a parent ID records ancestry. These identifiers do not label a cell type. |
| Cell state | The instantaneous values of a cell's chemical activities, polarity, shape, and division bookkeeping. State can change without the cell becoming a different biological type. |
| Chemical state | The pair $(a_i,b_i)$ of positive activator and inhibitor activities in cell $i$. Each is spatially uniform within that cell. Concentration/activity is distinguished from amount, $V_i a_i$ or $V_i b_i$. Neither species is an identity label or an identified gene. |
| Cell attributes and phenotype | Attributes are measured properties; phenotype here is their observable profile and behavior under specified conditions. We record log activities, polarity magnitude, cell axis ratio, and asphericity; separate pulse assays measure response behavior. Tension and adhesion computed directly from activity supply no independent evidence of identity. |
| Cell identity | A distinguishable profile that persists over a declared interval and has reproducible responses to specified challenges. Current persistence, recovery, exchange, and context assays test parts of this definition. High activator, a distinct color, a cluster, or contrast above a threshold alone does not establish identity. |
| Differentiation | Development of lasting differences between initially similar cells. **Chemical differentiation** means persistent differences in the reduced chemical variables; **identity differentiation** requires the broader attribute/response evidence above. A continuous range of states does not automatically constitute two cell types. |
| Geometric and chemical context | Geometric context includes position, volume, exposed cortex, and contacts; chemical context includes neighboring cells' activities and exchange conditions. Lineage and developmental history are recorded covariates, not extracellular signals. Context can support or rewrite a cell-associated chemical state. |
| Response behavior | The time course following a defined perturbation, compared with that same background's unperturbed continuation. Donor-nearer pulse responses are a measured similarity, not proof of equivalent identities or unchanged chemical transfer. |
| Fate and commitment | Fate concerns a cell's future developmental outcome or potential. Commitment would mean retention of that outcome after a specified change in instructive conditions. The research model supplies no lineage outcomes or fate switch, and does not establish commitment. |
| Memory and autonomy | Memory means a lasting effect of prior conditions after the initiating intervention ends; it can reside in a collective state. Autonomy would require retention under interventions that remove or standardize environmental support. Persistence on the original contact graph is not an autonomy test. |
| Inheritance | Transmission of an attribute or response profile through subsequent division. Tracking a parent ID, copying a variable at cleavage, or observing mature nondividing cells is not evidence that an emergent identity is inherited. |

### Polarity is orientation within a cell

**Polarity is a directional asymmetry in a cell's organization.** Our reduced apical–basal polarity variable is a vector $\mathbf p_i$: its direction points toward the modeled apical side, and $\lVert\mathbf p_i\rVert\in[0,1]$ measures the strength of that orientation. At zero magnitude the cell has no modeled polarity direction. Its magnitude is a normalized model measure, not a molecular concentration, probability, or apical surface fraction. “Apical” denotes the side selected by exposed-cortex cues; “basal” denotes the opposing direction. These are local cell directions, not the embryo's head–tail axis. The vector does not resolve molecular apical domains, basal machinery, or planar polarity.

Exposed cortex is the portion of a cell boundary not occupied by other cells. Its surface normals provide the geometric cue driving polarity. Unequal exposure can therefore create polarity at uniform chemistry. Polarity then modulates cortical tension directionally. This is a simplified representation motivated by [apical-domain experiments in mouse embryos](https://pmc.ncbi.nlm.nih.gov/articles/PMC5300053/), not a simulation of their molecular lineage program.

Polarity and elongation are different measurements: the polarity vector distinguishes an oriented side, while a shape axis describes the longest geometric extent. A polar cell need not be elongated, and an elongated cell need not have nonzero modeled polarity. The coefficient $\chi$ sets how strongly polarity changes tension; setting $\chi=0$ removes that mechanical action while retaining polarity dynamics. A polar cell also need not have a distinct chemical state or identity.

The **cell axis ratio** compares the largest and smallest principal lengths of its occupancy-weighted spatial covariance. **Asphericity** measures inequality of those covariance eigenvalues; both are insensitive to rotating the cell. Axis ratio one and asphericity zero mean isotropic second moments, which alone do not prove a spherical boundary. The **aggregate axis ratio** applies the corresponding measurement to the whole tissue and is a separate observable. The manuscript gives the exact cell-shape formulas.

### Organization, stability, and symmetry breaking

| Term | Operational meaning and limit |
|---|---|
| Initiation / formation | Development of chemical contrast from a specified near-uniform start with small perturbations. A mature-state restart tests formation on that geometry; it is not a new zygote trajectory. |
| Maintenance / persistence | Retention of a developed difference over a declared observation window. Recovery after perturbation, independence from context, and inheritance require separate tests. |
| Attractor / basin | An attractor is a state or set approached by the dynamics; its basin is the set of initial conditions that approach it. Frozen chemical assays test local attraction on a fixed measured graph. |
| Local bistability | Locally stable uniform and patterned chemical equilibria coexist on the same frozen graph. This is a whole-network result; it does not imply two autonomous fates inside every cell or full moving-system stability. |
| Signaling symmetry breaking | Development of differences from approximately uniform chemical activities. Our main contrast measure is the across-cell SD of $\log a$. It does not assign identities or by itself demonstrate spatial domains, spontaneous axis selection, or a Turing mechanism. |
| Geometry / shape symmetry breaking | Loss of a specified geometric symmetry, such as approximate rotational symmetry of the aggregate. Cell shapes, contacts, and aggregate axis ratios measure different scales of geometry; unequal chemistry does not establish a new global shape axis. |
| Spontaneous and emergent organization | A feature is emergent when the dynamics produce it rather than assigning it as a label or imposing its final form. Calling an axis spontaneous additionally requires ruling out a prescribed cue and checking bias from starting geometry, boundaries, and the numerical grid. The regulatory laws and mechanical couplings remain supplied assumptions. |

For example, a cell can retain high activity and recover after a pulse because its neighbors sustain a collective chemical pattern. That supports a persistent, context-dependent phenotype. It does not demonstrate isolated-cell memory, a specified developmental fate, or inherited identity. These distinctions govern the claims below and the [manuscript definitions](manuscript/introduction_methods.md#operational-definitions-state-phenotype-identity-and-context).

## Cell identity as emergent attributes

The [attribute-based developmental experiment](docs/attribute_development.md) completed two fresh zygote-to-t=90 branches: direct regulatory mechanical feedback and a matched no-feedback branch. We record log activator, log inhibitor, polarity magnitude, cell axis ratio, and asphericity as continuous attributes. Position, lineage, exposure, and contact geometry are context. Tension and adhesion are derived from activity and therefore do not count as independent evidence of differentiation. No clustering or number of identities is imposed.

The method now combines developmental controls, frozen-context stability/recovery, chemical-state exchange, and matched pulse responses under moving geometry. These tests ask separately whether differences form, persist, recover, or follow transferred chemistry. Environment-supported collective organization is distinguished from an autonomous cell identity. The [completed moving-history replication](docs/cell_response_moving.md#replication-across-developmental-histories) repeats exchange and response tests on seeds 8 and 9, using seed 7 as a historical reference. [Timestep refinement](docs/exchange_response_refinement.md) now tests the new histories' formation/retention and selected responses from identical physical starts.

The historical `Simulation` and dashboard still integrate a supplied bistable fate switch and can report A/B. They are available for reproducing earlier experiments. Their fate labels are not used by the current attribute-based identity studies, and GPU acceleration is not automatically enabled in the dashboard.

### Experimental stages and what they test

| Stage | Starting state and controlled change | Question answered |
|---|---|---|
| Development | One zygote, initially unit chemistry and zero polarity; matched feedback-on/off branches | Do continuous attributes become different during cleavage and changing shape? |
| Component controls | One common mature geometry; vary tension, adhesion, or polarity's mechanical action with either small chemical perturbations or a prepared pattern | Which mechanical term changes the opportunity to form a pattern, and does that differ from maintaining one? |
| Moving survival | Actual developed chemistry at t=90; switch feedback on or keep it off through t=150 | Can a developed pattern survive evolving mechanics, transport, and dilution? |
| Frozen recovery and bistability | Hold measured volumes and the contact operator fixed; perturb equilibria or change chemical initial conditions | Which chemical states are locally stable on that geometry? |
| Chemical exchange | Conservatively exchange both species between selected cells; compare fresh and frozen pre-relaxed exchanges | Do states return to their destination values, follow transferred chemistry, or reorganize collectively? |
| Matched response | At each moving background's t=210 endpoint, pulse one cell's activator by ±10%; subtract that background's unpulsed continuation | Does response behavior resemble the chemical donor or the original destination? |
| Environment assays | Remove transport, clamp a shared reservoir, or let a finite shared reservoir evolve | Are differences autonomous or supported by reciprocal environmental exchange? |

The mature exchange/response sequence uses prepared chemical states on developed geometry; it is not a second zygote-to-identity experiment. Chemical exchange leaves cell shape, polarity, position, and lineage in place. Reservoir assays are separate compartment models, not an extracellular field in the live embryo. Formation, maintenance, response similarity, and inherited identity therefore remain distinct claims. [Detailed assay methods](docs/attribute_development.md#how-the-identity-assays-build-on-development).

## Method and mathematical model

Read the equations as coupled rules: **geometry determines contacts and polarity cues; contacts transport chemicals; chemical activity modulates polarity, tension, and adhesion; mechanics changes geometry and volume; volume changes dilute or concentrate chemicals**. The model does not assume that any one link is sufficient to produce an organized embryo. The equations below describe the current attribute model; the supplied fate switch is identified separately as a historical alternative.

In the notation below, $i,j$ label cells, $\mathbf{x}$ is a location in the 3D box $\Omega$, a dot means a time derivative, and an integral adds a quantity over the box. A star on $V_i^\star$ denotes a prescribed target, not a measured volume. The occupancy function $h(\phi)$ is unrelated to the inhibitor, which is written $b_i$ here (`inhibitor` and sometimes `h` in the code). Adhesion strength is $A_{ij}$; geometric contact area is $\mathcal A_{ij}$. See [the annotated model](docs/model.md#reading-the-equations) for derivative notation and a term-by-term mechanical derivation.

### State variables and initial conditions

The simulation combines a continuous description of cell shape with a discrete description of signaling between cells. Each cell occupies a deformable region of a common three-dimensional computational domain; the contact network is reconstructed from those regions as they move and divide.

| Variable | Meaning | Initial zygote |
|---|---|---|
| $\phi_i(\mathbf{x},t)$ | Diffuse cell indicator: approximately one inside cell $i$, zero outside | Smooth sphere of radius 0.8 |
| $a_i(t), b_i(t)$ | Nonnegative activator and positive inhibitor activities | Both one |
| $\mathbf{p}_i(t)$ | Apical–basal polarity vector, pointing toward the apical side | Zero |
| $V_i^\star$ | Target cell volume | Measured initial zygote volume |
| Cell ID, parent ID, cycle state | Lineage and division bookkeeping | One founder cell |

All quantities are **dimensionless**. Activities are reduced regulatory variables, not identified genes or measured concentrations. Cell-cycle timing, constitutive laws, and noise amplitudes are supplied assumptions. The model tests their consequences; it does not derive living matter, metabolism, or the cell cycle from chemistry.

**Resource availability is an assumption.** Nutrients, biosynthetic precursors, and energy availability are assumed permissive throughout the simulated interval. Metabolism and resource limitation are not modeled. Cleavage partitions existing cell volume rather than increasing tissue biomass. The total preferred cell volume is conserved across divisions, while actual volumes can deviate under mechanics. Division timing, signaling reactions, and active mechanics do not depend on nutrient, oxygen, or ATP availability. Results therefore describe organization under permissive conditions, not nutrient-dependent growth, division arrest, or metabolic self-sufficiency.

Current developmental studies use $[-2.24,2.24]^3$ on a $72^3$ grid, with zero-normal-gradient mechanical boundaries and a sixteen-cell cap. The embryo outline is not prescribed. Initially uniform chemistry receives small amount-balanced partition perturbations at cleavage; there is no imposed chemical gradient in the baseline. Random division orientation is used within geometrically degenerate long-axis subspaces, unless the isotropic control is explicitly selected. The smaller $40^3$ core/dashboard default is a demonstration preset, not the current research resolution.

### Deformable cells and mechanical interactions

A phase field lets each cell change shape without treating it as a rigid sphere or prescribing its surface mesh. This representation is motivated by multicellular phase-field work such as [MorphoSim (2023)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9938209/); the energy and solver here are an independent simplified implementation.

Define a smooth occupancy function and an interface potential:

$$
h(\phi)=\phi^2(3-2\phi),
\qquad q(\phi)=\phi^2(1-\phi)^2,
\qquad V_i=\int_\Omega h(\phi_i)\,\mathrm{d}\mathbf{x}.
$$

**Occupancy is a smooth volume-counting rule.** It gives $h(0)=0$, $h(1)=1$, and $h(1/2)=1/2$: a grid location in the diffuse boundary contributes a fractional amount to cell volume. The cubic is the lowest-degree polynomial that also has zero slope at both endpoints. Consequently, volume-restoring forces proportional to $h'(\phi)=6\phi(1-\phi)$ act mainly at the boundary, rather than changing the uniform cell interior or empty space. Occupancy is not a probability of cell identity.

**The potential makes inside and outside preferable to an intermediate phase.** The nonnegative quartic $q$ has equal minima at 0 and 1 and a maximum at $1/2$. It penalizes a broad region of intermediate phase; the gradient term below opposes an infinitely sharp boundary. Their competition produces a finite-width interface. These polynomials are simple modeling choices with useful smoothness and symmetry, not unique biological laws.

![Occupancy rises smoothly from zero to one; the potential has minima at zero and one and a maximum at one half.](docs/images/occupancy-potential.svg)

The panels have different vertical scales: occupancy counts volume, whereas the potential contributes to interface energy. They are different functions with different jobs.

With regulatory states and the spatial tension coefficients frozen during each mechanical update, the baseline energy is

$$
\begin{aligned}
E={}&\sum_i\int_\Omega\gamma_i(\mathbf{x})
\left[\frac{\epsilon^2}{2}\lVert\nabla\phi_i\rVert^2+q(\phi_i)\right]\,\mathrm{d}\mathbf{x}\\
&+\sum_i\frac{K_V}{2V_i^\star}(V_i-V_i^\star)^2\\
&+\frac{R}{2}\sum_{i<j}\int_\Omega\phi_i^2\phi_j^2\,\mathrm{d}\mathbf{x}\\
&-\sum_{i<j}A_{ij}\int_\Omega q(\phi_i)q(\phi_j)\,\mathrm{d}\mathbf{x}.
\end{aligned}
$$

“Energy” here means a scalar mechanical cost used to generate shape-restoring forces. Lower values are preferred by the passive shape update. It is not ATP use, a developmental objective, or a claim that an embryo minimizes one fixed energy throughout development.

| Term, in equation order | Why this form? | What competes with it? |
|---|---|---|
| Interface cost | $\lVert\nabla\phi\rVert^2$ penalizes abrupt spatial changes; $q(\phi)$ favors inside/outside values. Together they assign a cost to cell boundary area. | Without a volume constraint, a cell can lower this cost by shrinking. |
| Volume penalty | The square penalizes both swelling and shrinkage and vanishes at the target. Dividing by $V_i^\star$ makes its derivative depend on relative volume error. | Surface and contact forces can sustain a small volume error because this constraint is soft. |
| Overlap repulsion | $\phi_i^2\phi_j^2$ is large where two cell interiors occupy the same region. The positive sign makes overlap costly. | It limits the overlap favored indirectly by cell attraction. |
| Interface attraction | $q(\phi_i)q(\phi_j)$ is appreciable only where both diffuse boundaries meet. The negative sign makes such contact favorable. | It competes with repulsion and interface cost; it does not impose a specific aggregate shape. |

The sum $i<j$ counts each cell pair once. The factors $1/2$ are coefficient conventions that simplify derivatives; they do not represent half a physical interaction. In particular, the attraction coefficient is not a contact area or a molecular binding constant.

| Symbol | Configuration name (default) | Interpretation when increased, with other parameters fixed |
|---|---|---|
| $\epsilon$ | `interface_width` (0.085) | Broadens diffuse interfaces and also changes their energy per area; it is not simply grid resolution. |
| $\gamma_0$ | `surface_tension` (1) | Raises the baseline cost of the cell boundary and the associated restoring force. |
| $K_V$ | `volume_stiffness` (12) | Resists fractional volume changes more strongly; can require smaller explicit time steps. |
| $R$ | `repulsion` (3) | Penalizes interpenetration more strongly. |
| $A_0$ | `adhesion` (4) | Strengthens the preference for overlapping diffuse interfaces. |
| $V_i^\star$ | Stored cell target, not a global configuration constant | Sets the preferred cell volume; daughters share the mother's target and do not grow between divisions. |

These tendencies do not guarantee monotonic changes in whole-embryo shape. With this energy normalization, an isolated flat equilibrium interface with constant $\gamma$ has energy per area $\gamma\epsilon\sqrt{2}/6$. Thus `surface_tension` is a model energy coefficient, not directly a measured cortical tension. The interface attraction is a contact surrogate, not a resolved cadherin or tight-junction model. [Derivation and force signs](docs/model.md#how-the-energy-generates-shape-change).

Mechanical evolution is overdamped, with mobility set to one:

$$
\frac{\partial\phi_i}{\partial t}
=-\frac{\delta E}{\delta\phi_i}+F_{\mathrm{furrow},i}.
$$

The minus sign means the passive field moves in the direction that reduces the instantaneous energy. The functional derivative $\delta E/\delta\phi_i$ asks how the energy changes when the shape field changes locally. “Overdamped” means there is no separate acceleration or momentum equation; mobility converts restoring force into the rate of shape change and is absorbed into the time scale here.

The furrow force is present only during division. Volume constraints are soft for nondividing cells; dividing mothers additionally receive the volume-preserving correction described below. Because signaling changes material properties and cytokinesis supplies active forcing, the complete simulation is not passive relaxation of a single fixed energy.

### Activator–inhibitor signaling on the changing contact graph

For interface shell $s_i=\phi_i(1-\phi_i)$, the raw contact weight is

$$
W_{ij}^{\mathrm{raw}}=\int_\Omega s_i s_j\,\mathrm{d}\mathbf{x},
\qquad W_{ii}=0.
$$

The live simulation defaults to **conservative concentration transport** (`signal_transport="conservative"`). Contacts below 2% of the largest current overlap are removed symmetrically. For complementary flat equilibrium interfaces,

$$
\frac{W_{ij}}{\mathcal A_{ij}}=\int_{-\infty}^{\infty}\phi^2(1-\phi)^2\,\mathrm{d}s
=\frac{\epsilon}{6\sqrt{2}}.
$$

We therefore estimate interface area and conductance by

$$
\widehat{\mathcal A}_{ij}=\frac{6\sqrt{2}}{\epsilon}W_{ij},\qquad
\ell_{ij}=\lVert\mathbf c_j-\mathbf c_i\rVert,\qquad
g_{ij}=\frac{\widehat{\mathcal A}_{ij}}{\ell_{ij}}.
$$

The shell weight $s_i$ is zero in the interior and exterior and largest at the interface, so $W_{ij}$ measures diffuse boundary overlap. The factor $6\sqrt{2}/\epsilon$ converts that overlap volume into an estimated area under the flat-interface assumption. Conductance increases with area and decreases with transport distance: a wider connection permits more exchange, while a longer path permits less. Here $\mathbf c_i$ is the occupancy-weighted cell center. The area estimate is a calibrated closure, not an exact face reconstruction. Gaps, overlapping or curved interfaces, and nonorthogonal center-to-face directions can bias transport. Coincident centers with positive contact are rejected.

With measured cell volumes $M=\operatorname{diag}(V_i)$ and $K=\operatorname{diag}(G\mathbf1)-G$,

$$
\Delta_V=-M^{-1}K,\qquad
(\Delta_V c)_i=\frac{1}{V_i}\sum_jg_{ij}(c_j-c_i).
$$

The difference $c_j-c_i$ sends regulator from higher to lower concentration. Dividing the net incoming amount flux by $V_i$ converts it to a concentration rate; the same incoming amount changes a small cell's concentration more than a large cell's. $M$ stores these volume capacities, and $K$ assembles the equal-and-opposite pairwise fluxes.

Exchange preserves constants and total amount $\sum_i V_i c_i$ on frozen geometry. Isolated cells have zero exchange. The Gierer–Meinhardt equations on moving compartments are

$$
\begin{aligned}
\frac{\mathrm d(V_i a_i)}{\mathrm dt}
 &=V_i\left(\frac{a_i^2}{b_i}-a_i\right)+D_a\sum_jg_{ij}(a_j-a_i),\\
\frac{\mathrm d(V_i b_i)}{\mathrm dt}
 &=V_i\beta(a_i^2-b_i)+D_b\sum_jg_{ij}(b_j-b_i).
\end{aligned}
$$

| Reaction term | Meaning and modeling choice |
|---|---|
| $a_i^2/b_i$ | Activator promotes its own production through the quadratic numerator; inhibitor suppresses that production through the denominator. The reciprocal law is a simplified inhibition rule requiring $b_i>0$. |
| $-a_i$ | Linear activator turnover; its coefficient is set to one by the chosen time scale. |
| $\beta a_i^2$ | Activator induces inhibitor, closing the negative feedback loop. |
| $-\beta b_i$ | Linear inhibitor turnover. The same $\beta$ scales inhibitor production and loss, changing its reaction time scale without shifting the positive uniform equilibrium $(1,1)$. |
| $D_a,D_b$ | Transport strengths (`signal_da`, `signal_dh`). Faster inhibitor transport can oppose broad activation while allowing localized activation; it does not guarantee a pattern on the available graph. |

The powers and coefficients specify a minimal reduced feedback system, not identified biochemical reaction steps. A quadratic production law is one way to supply nonlinear self-amplification; it is not derived from the mere existence of an activator–inhibitor loop. The model must test whether this assumed loop explains the observed outcome.

Thus concentration equations include dilution $-c_i\dot V_i/V_i$. Each step advances reaction/exchange on frozen pre-step geometry with positivity-preserving SSP-RK2 substeps, advances mechanics, then rescales concentrations by $V_i^{old}/V_i^{new}$. This is first-order splitting of the coupled moving problem; the internal RK2 solver does not make the full simulation second order. Abscission partitions measured amounts conservatively. Externally clamped signal experiments instead supply/remove regulators to maintain the imposed concentration.

Defaults are $\beta=2$, $D_a=0.02$, $D_b=0.4$ in model length-squared/time units. These are exploratory diffusivities, not calibrated molecular measurements or a conversion of the old exchange rates. The explicit `random_walk` option preserves the historical model; old checkpoints without the selector restore that option and their recorded coefficients. Polarity neighbor alignment still uses normalized orientation averaging, separately from molecular transport.

### Discrete linear stability before 3D simulation

The positive homogeneous equilibrium is $(a_*,b_*)=(1,1)$, with reaction Jacobian

$$
J=\begin{pmatrix}1&-1\\2\beta&-\beta\end{pmatrix},
\qquad \operatorname{tr}J=1-\beta,
\qquad \det J=\beta.
$$

Each Jacobian entry is the response of one reaction rate to a small change in one activity. The positive upper-left entry represents self-amplification, the negative upper-right entry inhibition, and the positive lower-left entry inhibitor induction. Local stability means a small perturbation in a well-mixed cell decays through their combined action, despite the activator's individual positive feedback. Local kinetics are stable for $\beta>1$. For frozen-geometry spectral analysis, the code diagonalizes

$$
S=M^{-1/2}KM^{-1/2}.
$$

Its eigenvalues are nonnegative, have inverse-length-squared units, and equal the eigenvalues of $-\Delta_V$. Physical right modes are $M^{-1/2}\mathbf q_k$, where $\mathbf q_k$ are orthonormal eigenvectors of $S$. There is no normalized upper bound of two. Each mode has a two-variable linear system

$$
M_k=J-\lambda_k\operatorname{diag}(D_a,D_b),
\qquad r_k=\max\operatorname{Re}\operatorname{eig}(M_k).
$$

A graph mode is a pattern of cell-to-cell activity differences, and $\lambda_k$ measures how strongly transport damps that mode. The zero mode is constant on each connected component. The rate $r_k$ tells whether a small perturbation grows ($r_k>0$) or decays ($r_k<0$); early linear growth is proportional to $\exp(r_k t)$. A mode with small positive growth may need much longer than one cell cycle to become visible.

A diffusion-driven instability requires stable local kinetics and at least one **supported nonzero graph eigenvalue** with $r_k>0$. The determinant is

$$
\det M_k=\beta+(\beta D_a-D_b)\lambda_k+D_aD_b\lambda_k^2.
$$

For the defaults $\beta=2$, $D_a=0.02$, and $D_b=0.4$, the unstable interval is $6.492189<\lambda_k<38.507811$. Its physical units and coefficients stay fixed during refinement; whether discrete eigenvalues enter the interval must converge and is not guaranteed on a coarse mesh. The dashboard and batch preflight use fixed unit-box meshes, followed by actual live-geometry spectra. Every abscission records spectra and volume-weighted mode transfer.

The [calibrated-contact refinement experiment](docs/live_transport.md) uses the same contact adapter on manufactured 3D slabs with unit transverse area. First-mode error decreases from 1.27% at 8 compartments to 0.020% at 64 (observed final order 1.998). The continuum has one unstable mode; 8 compartments incorrectly support two, while 16, 32, and 64 recover one. Measured small-perturbation growth agrees with the discrete prediction within $7.6\times10^{-11}$. This establishes convergence for flat complementary contacts, not arbitrary embryo geometry.

Frozen-geometry spectra describe the signaling subsystem near $(1,1)$. Volume changes introduce dilution; graph changes, mode mixing, and finite growth time require time-dependent analysis. These spectra are not a stability proof for the coupled chemistry–polarity–mechanics system. Historical normalized-graph experiments remain documented in [graph_signaling.md](docs/graph_signaling.md).

### From chemical activity to mechanics

The current attribute model uses an instantaneous, bounded material response:

$$
\begin{aligned}
r_i&=\tanh(a_i-1),\\
\gamma_i^0&=\gamma_0(1+c_\gamma r_i),\\
A_{ij}&=A_0(1+c_A r_i r_j),\qquad A_{ii}=0.
\end{aligned}
$$

Subtracting the homogeneous activator reference, one, makes baseline chemistry mechanically unbiased. The tanh bounds the response and prevents unbounded material coefficients. Higher activity increases baseline tension; cells with responses of the same sign have stronger attraction than opposite-sign pairs. The contrasts $c_\gamma=0.25$ and $c_A=0.35$ retain configuration names `fate_tension` and `fate_adhesion` for compatibility, but **no fate variable enters this calculation**. At these contrasts, tension and attraction remain positive.

The response $r_i$ has no differential equation, stored memory, or threshold classification. These constitutive laws assume prompt activity-dependent mechanics; any persistent phenotype must arise from the coupled chemical/geometric dynamics. In the no-feedback control, tension and attraction retain their baseline values, and polarity cannot change mechanical tension; chemical and polarity dynamics remain active. [Detailed definitions and developmental controls](docs/attribute_development.md).

### Historical fate-switch alternative

The core/dashboard model instead integrates a supplied signed switch after four cells are present:

$$
\dot f_i=r_f[f_i-f_i^3+g_a(a_i-1)].
$$

Its two stable states at $f=\pm1$ are an assumption. It derives material responses from $\tanh(f_i)$ and reports A/B using $f_i>0.55$ or $f_i<-0.55$. Earlier experiments showed that transient chemical differences can select these labels without sustained Turing patterning. This motivated removing the switch from identity-focused research. The attribute model keeps only zero-valued compatibility arrays for checkpoint bookkeeping; no fate drift or A/B classification is performed. [Historical equations and interpretation](docs/model.md#regulatory-activity).

### Apical–basal polarity and directional mechanics

A cell's exposed cortex supplies a local orientation cue. Define the unoccupied-cortex weight, outward normal, and average cue by

$$
\begin{aligned}
e_i(\mathbf{x})&=1-\operatorname{clip}\left(2\sum_{j\ne i}h(\phi_j),0,1\right),\\
\mathbf{n}_i&=-\frac{\nabla\phi_i}{\max(\lVert\nabla\phi_i\rVert,10^{-10})},\\
\mathbf{q}_i&=\frac{\int_\Omega s_i e_i\mathbf{n}_i\,\mathrm{d}\mathbf{x}}
{\int_\Omega s_i\,\mathrm{d}\mathbf{x}}.
\end{aligned}
$$

An isolated sphere has no net cue. Contacts can create a cue toward free cortex without specifying a global embryo axis. The polarity vector evolves as

$$
\dot{\mathbf{p}}_i=
\alpha\frac{2a_i}{1+a_i}\mathbf{q}_i
+\eta(\Delta_{\mathrm{rw}}\mathbf{p})_i
-(\mu+\lVert\mathbf{p}_i\rVert^2)\mathbf{p}_i.
$$

Here $\alpha$ (`polarity_rate`, 1) controls response to free cortex, $\eta$ (`polarity_alignment`, 0.25) controls averaging toward neighboring orientations, and $\mu$ (`polarity_decay`, 0.5) removes polarity when cues are absent. The activity factor $2a/(1+a)$ equals one at $a=1$ and saturates at two, so signaling modulates rather than indefinitely amplifies the geometric cue. The cubic vector term $-\lVert\mathbf p\rVert^2\mathbf p$ increasingly opposes large polarity. These are response and saturation assumptions, not a resolved polarity-protein network. Each numerical update caps the vector magnitude at one. Daughters inherit the mother's vector and subsequently adapt to their new geometry.

Polarity modifies cortical tension spatially:

$$
\begin{aligned}
\gamma_i(\mathbf{x})&=\gamma_i^0
[1-\chi\mathbf{p}_i\cdot\widehat{\mathbf{r}}_i(\mathbf{x})],\\
\widehat{\mathbf{r}}_i&=\frac{\mathbf{x}-\mathbf{c}_i}
{\sqrt{\lVert\mathbf{x}-\mathbf{c}_i\rVert^2+\epsilon^2}},\\
F_i^{\mathrm{surface}}&=\epsilon^2\nabla\cdot(\gamma_i\nabla\phi_i)-\gamma_i q'(\phi_i).
\end{aligned}
$$

Here $\mathbf{c}_i$ is the cell centroid and $\chi$ (`polarity_tension`, 0.35) is the strength of the directional tension contrast. The dot product selects position relative to the apical direction; the $\epsilon^2$ in the denominator prevents division by zero at the centroid. Positive $\chi$ lowers effective tension on the apical side and raises it on the basal side; $\chi<1$ keeps tension positive. Conservative face fluxes include the spatial gradient of tension. Centroids and polarity are held fixed within each mechanical update.

[Nissen et al. (2018)](https://elifesciences.org/articles/38407) motivates investigating polarity-dependent morphology, but this particular cue and tension law are project-specific assumptions. There is no planar cell polarity or resolved apical protein network. Geometry can polarize cells even with uniform activator, so local polarity is not evidence of a chemically selected global axis.

### Shape-aligned division and progressive cytokinesis

Division uses the largest-eigenvalue direction of the occupancy-weighted cell covariance:

$$
\mathbf{C}_i=\frac{1}{V_i}\int_\Omega
(\mathbf{x}-\mathbf{c}_i)(\mathbf{x}-\mathbf{c}_i)^{\mathsf T}
h(\phi_i)\,\mathrm{d}\mathbf{x}.
$$

The spindle follows this long axis; the cleavage plane is perpendicular to it and initially bisects occupancy. Nearly equal maximal eigenvalues define a subspace in which orientation is sampled without privileging grid eigenvectors. This is a Hertwig-style shape rule, not a tensile-stress calculation or a simulation of spindle alignment. Experiments demonstrate a role for ring mechanics in achieving long-axis alignment in particular embryos, but their cell-rotation mechanism is not implemented here. [Middelkoop et al. (2024)](https://pubmed.ncbi.nlm.nih.gov/38870057/)

A cycle timer initiates division while retaining one mother field. With onset $t_0$ and nominal duration $T$, constriction follows

$$
u=\operatorname{clip}((t-t_0)/T,0,1),
\qquad S(u)=u^2(3-2u),
\qquad R_{\mathrm{ring}}(t)=R_0[1-S(u)].
$$

For signed plane distance $z_i$ and distance $r_\perp$ from the spindle axis, the prescribed equatorial force is

$$
F_{\mathrm{furrow},i}=-\kappa S(u)
\exp\left(-\frac{z_i^2}{2\epsilon^2}\right)
\frac{1+\tanh[(r_\perp-R_{\mathrm{ring}})/\epsilon]}{2}
h'(\phi_i).
$$

The smooth ramp starts and ends with zero slope, avoiding an abrupt switch in forcing. The Gaussian confines the force near the cleavage plane; the tanh factor selects tissue outside the shrinking ring; and $h'(\phi_i)$ localizes its action to the cell boundary. The negative sign removes occupancy from that equatorial region, while volume preservation redistributes it into the lobes. $T$ (`cytokinesis_duration`, 0.9) sets the nominal constriction time and $\kappa$ (`ring_strength`, 8) sets force amplitude. Neither alone determines actual abscission time, which also depends on the resolved neck.

This is a contracting-ring surrogate, not a resolved actomyosin network. After each mechanical update, an interface-local scalar correction preserves the mother's measured onset volume. That correction is a numerical constraint, not a hydrostatic pressure solve.

Abscission requires elapsed constriction time, a sufficiently thin resolved neck, low prospective daughter overlap, and two nontrivial lobes. Only then does an occupancy partition replace the mother:

$$
h(\phi_1)=w\,h(\phi_{\mathrm{mother}}),
\qquad h(\phi_2)=(1-w)h(\phi_{\mathrm{mother}}),
\qquad w=\frac{1+\tanh(z_i/\epsilon)}{2}.
$$

Occupancy is conserved pointwise by this construction, up to numerical precision. Daughter target volumes follow their measured lobe fractions; there is no growth between divisions. Signal partition perturbations preserve measured volume-weighted regulator amounts at cleavage; reactions can change these amounts, whereas exchange and mechanical dilution conserve them. Unresolved divisions remain active rather than being forcibly cut. See [cytokinesis details and validation](docs/cytokinesis.md).

### Numerical workflow and research scales

Each attribute-model step reconstructs contacts and volumes, advances reactions/exchange on that geometry, updates polarity, computes material coefficients from the updated activator, and advances mechanics. Concentrations are then multiplied by $V_i^{\mathrm{old}}/V_i^{\mathrm{new}}$, preserving each species' amount during the volume change. CPU development additionally handles completed and newly scheduled divisions. Chemistry uses positivity-restricted SSP-RK2 substeps; polarity and mechanics use explicit updates. **The complete coupled scheme is first-order operator splitting**, even though frozen-geometry chemical integration is second order.

| Parameter group | Current research values |
|---|---|
| Domain and population | $[-2.24,2.24]^3$; $72^3$ grid; sixteen-cell cap |
| Time | Development/component pilots: $\Delta t=0.0075$; accepted mature response/GPU studies: $\Delta t=0.00375$; response refinement: $0.001875$ |
| Mechanics | $\epsilon=0.085$, $\gamma_0=1$, $K_V=12$, $R=3$, $A_0=4$ |
| Signaling | $\beta=2$, $D_a=0.02$, $D_b=0.4$; partition-noise scale 0.001 |
| Material response | $r_i=\tanh(a_i-1)$, $c_\gamma=0.25$, $c_A=0.35$; no fate switch |
| Polarity | $\alpha=1$, $\eta=0.25$, $\mu=0.5$, $\chi=0.35$ |
| Division | Mean cycle interval 2; nominal constriction duration 0.9; maximum 16 cells |
| Mature response observations | Every 0.15 units; restart checkpoint every three units; 60-unit moving response horizon |

These are dimensionless experimental values, not calibrated biology or universal accuracy guarantees. Signaling growth time, cell-cycle time, and mechanical relaxation time must be compared. [Config](embryo/model.py) and [default.json](examples/default.json) describe the smaller historical demonstration defaults; frozen protocols record each research configuration and its controls.

Every-step production screens require per-cell volume error below 5%, equivalent radius at least four grid spacings, zero clipping, and positive finite chemistry. Sampled boundary occupancy must remain below 0.01. Clipping is a numerical safeguard and its activation fails these studies; staying finite is insufficient. Numerical failures are not classified as biological loss of organization.

Shape measurements use the capped aggregate occupancy $\rho=\min(\sum_i h(\phi_i),1)$. If its spatial covariance eigenvalues are $\ell_1\le\ell_2\le\ell_3$, the principal axis ratio is $\sqrt{\ell_3/\ell_1}$, equal to one for a sphere. Continuous chemical attributes, lineage, polarity, graph spectra, and volume errors are recorded alongside shape. Elongation alone cannot identify whether its cause is cleavage, regulatory feedback, or numerical anisotropy.

### What this model can explain, and what remains missing

The implemented system tests whether local activation and inhibition, a changing contact network, and polarity-dependent mechanics are sufficient for specific forms of organization. Formation and maintenance must be tested separately: polarity-mediated mechanics can suppress initiation on a given geometry while established nonlinear chemical patterns persist. Continuous heterogeneity, donor-nearer response behavior, and local polarity each require further evidence before being called cell identity or a selected developmental axis.

**Blastocoel cavitation is not implemented.** Repulsion and cell-volume penalties can leave geometric gaps, but they do not explain accumulation of pressurized extracellular fluid. Na⁺/K⁺-ATPases actively transport ions; aquaporins conduct water passively; a low-leak epithelial barrier permits fluid accumulation. Experimental work also links pump signaling to tight-junction function. [Giannatselis et al. (2011)](https://pubmed.ncbi.nlm.nih.gov/21901128/)

A future cavity model needs solute and water balances, barrier permeability, lumen pressure coupled to cell mechanics, and communication between microlumens. Hydraulic opening of contacts and microlumen coarsening can contribute to cavity positioning, as shown in mouse embryos. [Dumortier et al. (2019)](https://pubmed.ncbi.nlm.nih.gov/31371608/) A centered cavity can preserve rotational symmetry; selecting its position is a separate explanatory target. Signaling and polarity could regulate these processes, but cannot substitute for fluid transport equations.

The remaining requirements include convergence of the full attribute-based developmental trajectory and geometric transport closure, exchange-specific refinement, additional histories, context/neighbor-dependence and inheritance tests, and biological calibration. Completed component, conservation, and backend checks do not replace these. References motivate individual mechanisms; none establishes that the combined implementation reproduces an actual embryo.

## Run

From this directory, using Python 3.10 or newer:

```bash
python -m pip install -e '.[test]'
```

### Current attribute-based research workflow

Development is run on CPU, using the attribute model rather than the historical core CLI:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.attribute_development prepare \
  --output outputs/attribute-development-repeat
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.attribute_development run \
  --output outputs/attribute-development-repeat
```

The moving-history replication requires the completed survival, frozen-endpoint/exchange, and GPU-validation artifacts already present locally. It uses a separate output directory and verifies their hashes:

```bash
CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python -m embryo.exchange_response_histories prepare \
  --output outputs/exchange-response-histories-repeat

CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python -m embryo.exchange_response_histories run \
  --output outputs/exchange-response-histories-repeat
```

The GPU path additionally requires the compatible PyTorch/CUDA environment and CUDA compiler described in [the backend method and validation](docs/model.md#resident-gpu-backend-validation). Device selection above is specific to this machine's dedicated GTX 1080 Ti. It is not a portable GPU-index convention. Preparation requires a fresh directory; restart uses `run` on the existing prepared directory without changing its frozen sources or inputs. Full scientific outputs are local artifacts, not included by installing the Python package.

### Historical core dashboard

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.dashboard
```

Open [http://127.0.0.1:8765](http://127.0.0.1:8765). The dashboard advances the actual 3D simulation while displaying shaded cell surfaces, signaling variation, cell identities, shape anisotropy, and numerical diagnostics. **Show → Cell surfaces** renders the actual closed cell boundaries; **Point samples** retains the previous display. Restart an existing server and refresh the browser after updating. The [surface-rendering documentation](docs/dashboard.md#closed-cell-surfaces) explains mesh extraction, cutaway behavior, and how to inspect an existing checkpoint.

- **Run** starts the solver; **Pause** stops at a numerical step boundary; **Resume** continues the same state and random streams.
- Edit parameters before running, or pause, edit, then **Reset** to apply them to a new zygote. **Defaults** restores the form values; Reset applies them.
- Orbit, zoom, change cell colors, display polarity/contact links, and cut through the embryo. Scrub retained frames to inspect earlier states; **Follow live** returns to the newest snapshot.
- Download retained visualization frames and their parameters as JSON. This download is not a restart checkpoint. The session is held in memory and ends when the server stops.

The server binds only to this computer and needs no browser dependencies or internet connection. Use `--config examples/default.json` for an initial configuration or `--port 8766` for another port. Reference graph stability is computed before stepping; the dashboard also shows the current contact graph's growing-mode count. See [dashboard controls and limits](docs/dashboard.md).

### Historical core batch simulation and offline playback

```bash
# Analyze finite reference graphs before a coupled run (no 3D simulation):
OPENBLAS_NUM_THREADS=1 python -m embryo.graph_analysis --output outputs/my-graph-analysis
# Every coupled run also writes a finite-graph preflight before its first step:
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo --output outputs/my-run
```

NumPy, SciPy, and Matplotlib are the runtime dependencies. The thread settings avoid excessive BLAS threading overhead for the small contact matrices. Existing output directories are never overwritten.

Open `outputs/my-run/viewer.html` in a browser. It is self-contained and works offline:

- Drag to rotate; scroll to zoom.
- Play or scrub the timeline from the zygote through cleavage and differentiation.
- Move the cutaway plane to inspect internal cells.
- Color by fate, activator, inhibitor, or cell ID; show apical polarity arrows.
- Follow the principal axis ratio and volume diagnostics.

Signal colors are relative to the homogeneous activity of one, not automatically rescaled to exaggerate small variations. The viewer displays samples of the simulated cell surfaces. Its apparent surface texture is a rendering approximation; it does not enter the dynamics.

## Experiments

The commands immediately below exercise the historical core/fate model. Current identity-focused protocols use `AttributeSimulation` and the research workflow above. Older experiments remain useful controls and numerical evidence, but their A/B outcomes are not reinterpreted as emergent identities in the new model.

```bash
# Same seed and parameters, without fate/polarity-dependent mechanical properties:
OPENBLAS_NUM_THREADS=1 python -m embryo --no-feedback --output outputs/no-feedback

# Cleavage and mechanics, with signaling, polarity, and differentiation disabled:
OPENBLAS_NUM_THREADS=1 python -m embryo --mechanics-only --output outputs/mechanics

# Separate signaling/polarity ablations:
OPENBLAS_NUM_THREADS=1 python -m embryo --no-signaling --output outputs/no-signaling
OPENBLAS_NUM_THREADS=1 python -m embryo --no-polarity --output outputs/no-polarity

# More cells and finer geometry (substantially more expensive):
OPENBLAS_NUM_THREADS=1 python -m embryo --grid 64 --max-cells 32 --steps 1200 --output outputs/larger

# Customize coefficients with a JSON file:
OPENBLAS_NUM_THREADS=1 python -m embryo --config examples/default.json --output outputs/custom

# Small paired-seed ablation screen; this is not statistical validation:
OPENBLAS_NUM_THREADS=1 python -m embryo.experiments --seeds 7 11 23 --output outputs/screen

# Deterministic four-cell spatial/time refinement checks:
OPENBLAS_NUM_THREADS=1 python -m embryo.convergence --output outputs/convergence.json
```

`--no-feedback` retains signaling and polarity dynamics but disables their fate/polarity-dependent mechanical effects. `--mechanics-only` also disables signaling, polarity, and differentiation. Mechanical and regulatory RNG streams are separate. Shape-dependent orientation and mechanically gated abscission can change division directions and timing across feedback controls, even at a fixed seed.

### Signaling growth versus cleavage timescale

The timescale study compares longer signaling on frozen actual 4-, 8-, and 16-cell graphs, coupled runs with normal and doubled cell-cycle intervals, and replay of an identical contact history with stretched timing. Every-step graph recording allows the replay to be checked against the original coupled signals before interpreting timing effects. See [the protocol, commands, and results](docs/timescales.md).

In the first seed-7 screen, the normal coupled run reaches activator standard deviation 0.1 at time 36.14; doubling the cycle interval delays this to 44.21. Measured from first reaching 16 cells, the delays are nearly equal: 23.76 and 23.34. Thus the original weak time-15 signals primarily reflect insufficient growth time in this case. The frozen 4-cell graph suppresses fluctuations, while the 8- and 16-cell graphs amplify them on different timescales. Multi-seed robustness, mechanical/grid convergence, and fate persistence remain untested by this screen.

### Larger domain at unchanged cell resolution

The new conservative model has a [fixed-spacing domain-size protocol](docs/domain_size.md). The larger dashboard preset uses `grid=56`, `extent=2.24`, preserving the original voxel spacing 0.08 while adding 0.64 length units on each side:

```bash
python -m embryo.dashboard --config configs/large_domain.json
```

This starts a new zygote; the preset runs to $t=30$. The completed conservative-model study compares grids 40, 48, and 56 from the same initial state. Through $t=30$, peak boundary occupancy falls from $7.39\times10^{-4}$ to $1.90\times10^{-11}$; the two larger domains agree within 0.0000141% in axis ratio and 0.000948% in final common-box cell fields. All three boxes pass the boundary screen in this window. Increasing the domain does not improve cell resolution or establish clearance for longer runs. See the protocol for boundary and domain-sensitivity checks.

The $t=180$ domain continuation is complete. The original box crosses the boundary threshold at $t=131.4$; the two larger boxes pass their screens and agree within 0.000673% in axis ratio. The 56³ run retains 12 A / 4 B labels while activator contrast decays almost to zero, so persistent labels do not establish persistent signaling. Full results are in `outputs/domain-conservative-extended/RESULTS.md`. See [continuation and monitoring](docs/domain_size.md#long-time-continuation-protocol).

### Coupled spatial and temporal resolution

The [completed resolution screen](docs/coupled_resolution.md) tests grids 56³, 72³, and 88³ at fixed half-width 2.24 and interface width 0.085, with a separate time-step sweep. Every grid samples the same manufactured 16-cell geometry directly. Frozen-geometry spectra agree within 0.000193% between the finer grids. The five coupled cases are complete: 11/12 checks pass, but the finest spatial occupancy-field discrepancy is 1.43%, above the 1% criterion. The completed [projection-error audit](docs/projection_audit.md) finds 1.49% discrepancy even before evolution under the original linear reconstruction; higher-order final comparisons range from 0.24% to 0.32%. This identifies substantial measurement contamination. The original failed result remains unchanged. An [independent 112³ confirmation](docs/refinement_confirmation.md) passed all 11 prospectively specified checks; finest-pair field differences are 0.162–0.207%. This supports short-time consistency for the manufactured state, not developmental convergence.

### Discrete-to-continuum transport bridge

A new fixed-domain benchmark checks conservative compartment exchange against the same continuum diffusion and Gierer–Meinhardt equations at increasing spatial resolution. It uses explicit compartment volumes and face conductances, with a sparse operator that preserves molecular amount under pure diffusion. This is a transport verification experiment; its compartments are numerical elements, not a simulation of tens of thousands of biological cells.

```bash
OPENBLAS_NUM_THREADS=1 python -m embryo.continuum --output outputs/continuum-bridge
```

Across 64 to 32,768 compartments, diffusion error falls from 0.002871 to 0.00004994 with approximately second-order convergence. Coarse grids predict ten unstable modes; finer grids recover the continuum prediction of seven. Fixed normalized exchange rates do not approximate the same bulk diffusivity under this refinement. The live simulator now uses conservative transport through a calibrated diffuse-contact adapter; its geometry closure has separate validation requirements. See [equations, protocol, results, and the remaining bridge to changing geometry](docs/continuum_bridge.md).

The completed irregular-domain bridge verifies diffusion on a **3D L-shaped domain with unequal compartment volumes**. It uses exact continuum cell averages, a uniform-mesh control, and a separate pulse traveling between the two arms:

```bash
OPENBLAS_NUM_THREADS=1 python -m embryo.irregular --output outputs/irregular-bridge
```

Across 48 to 24,576 active compartments, graded-mesh error falls from 0.003924 to 0.0001186; the final observed order is 1.992. Pure-diffusion mass drift stays below $1.5\times10^{-14}$ in the refinement runs. The pulse reaches the other arm and approaches the correct volume-weighted equilibrium. All ten default acceptance checks pass. These results establish smooth-solution transport convergence on this fixed nonconvex domain; nonlinear pattern formation, moving geometry, and cell/field coupling remain separate tests. See [the irregular-domain protocol and results](docs/irregular_bridge.md).

The [irregular-domain signaling study](docs/irregular_signaling.md) now computes the complete discrete weighted spectrum and verifies early Gierer-Meinhardt growth:

```bash
OPENBLAS_NUM_THREADS=1 python -m embryo.irregular_signaling \
  --output outputs/irregular-signaling-repeat
```

Refinement changes the predicted unstable-mode count from **9 to 7 to 4 to 4**. On the finest mesh, the fastest predicted growth rate is 0.215193833 and the nonlinear solver measures 0.215193831. Both growing and decaying probes pass verification, and all nine default checks pass. Agreement of the finest discrete counts is refinement evidence; an exact continuum instability count remains a separate question.

The [nonlinear persistence and convergence study](docs/nonlinear_bridge.md) now follows one physical initial perturbation through $t=100$ on the same graded L-domain:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.nonlinear \
  --output outputs/nonlinear-bridge-repeat
```

All **eleven default criteria pass**. The final two-species spatial difference decreases from 5.43% ($n=8$ versus $16$) to **1.26%** ($n=16$ versus $32$). Independent time refinement gives a **0.095% maximum sampled difference** between the two smaller steps. Late field change stays below 0.41%, while an equal-diffusivity control returns to uniformity. The benchmark uses conservative fine-to-coarse averaging and an independently tested implicit-diffusion solver; it saves concentration snapshots, histories, and comparison figures. This verifies finite-window nonlinear signaling persistence for one prescribed perturbation on fixed geometry. Moving-domain dilution and amount balance are now tested below; initial-condition robustness, fate persistence, and emergent shape remain separate tests.

The dashboard’s **Continuum signaling** workspace now plays the saved concentration fields with 3D cutaways, X/Y/Z cross-sections, activator/inhibitor colors, run selection, and convergence charts. Start with `python -m embryo.dashboard --benchmark outputs/nonlinear-bridge-verified`; an existing server needs restarting. Playback is separate from the live embryo solver. See the [dashboard guide](docs/dashboard.md#continuum-signaling-nonlinear-benchmark-playback).


The [moving-domain transport test](docs/moving_domain.md) now verifies conservation and dilution during prescribed 3D expansion and deformation:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.moving \
  --output outputs/moving-domain-repeat
```

All **eleven checks pass**. Isotropic expansion, unequal axis expansion, and volume-preserving stretching conserve total amount to within $2.2\times10^{-14}$ relative drift. Uniform concentration follows inverse volume; spatial and temporal errors approach second-order convergence. The largest finest-mesh transport error is **0.01069%** relative to the exact full concentration field. Omitting dilution instead creates an artificial **82.21% amount gain** in the expanding control. This tests prescribed material motion with fixed connectivity, not emergent shape. Conservative remapping at live-cell abscission is now implemented; general remeshing remains separate.


The [persistent shape pilot](docs/shape_persistence.md) now tests the coupled embryo model through $t=90$, including matched mature-state feedback removal and a separate zygote developed without feedback. All four trajectories retain an identifiable long axis, with final shape ratios **1.323–1.331**. Full feedback does **not** exceed the controls; all final axes remain within **1.37°** of the first-cleavage axis. This is evidence of persistent elongation, not a demonstrated activator–inhibitor-specific shape axis. Every run also narrowly fails the four-grid-spacing smallest-cell-radius screen, so mechanical refinement remains necessary. The experiment exports offline 3D playbacks, axis histories, and a matched-camera comparison; [the protocol](docs/shape_persistence.md#reproduce) reproduces all controls.


The [long-time continuation](docs/shape_persistence.md#long-time-continuation) extends the three mature branches to $t=180$ without increasing the timestep. A separate [controlled chemical-patch experiment](docs/signal_patch.md) tests the response of the existing mechanics to imposed signaling, compares default and stronger fate-dependent tension against matched uniform controls, and then releases the clamp. It is explicitly a forced-response diagnostic, not evidence of spontaneous symmetry breaking.

The completed extension gives final axis ratios **1.328** (full feedback) and **1.341** (both ablations), with no positive feedback-specific excess. The full run crosses the boundary-occupancy screen at $t=124.2$, and all branches retain the smallest-cell resolution failure. Longer time has therefore exposed a domain-size limitation as well as modest shape evolution; [the continuation results](docs/shape_persistence.md#completed-continuation-results) report both.

The completed [patch-and-release pilot](docs/signal_patch.md#completed-pilot-results) also shows modest deformation: patched aggregates remain near axis ratio **1.329** at the end of forcing, at both tested fate–tension strengths. The patch/uniform difference increases with coupling, partly because the stronger uniform control becomes less elongated. All patch runs retain the resolution failure; the four feedback-enabled branches also cross the boundary screen. These are historical fate-model results. Later larger-domain and resolution studies address selected numerical issues but do not validate activator–inhibitor-specific shape formation.

The mechanics kernel now avoids discarded Laplacians and repeated occupancy calculations. It reproduces the previous kernel's tested trajectories exactly; the mature 16-cell benchmark measured **1.34× throughput (25% less runtime)**. This optimization also applies to the live dashboard after restarting its server. The [benchmark and protocol](docs/signal_patch.md#equation-preserving-kernel-optimization) document the scope and reproduction commands.


### Current computing methodology

| Path | Role | Validation scope |
|---|---|---|
| Python/NumPy/SciPy attribute model | CPU development and reference equations | Original scientific protocols and developmental diagnostics |
| Native C++/OpenMP | Faster mature CPU mechanics and geometry/polarity reference | Complete tested control/pulse replays; retains the standard attribute checkpoint format |
| Resident PyTorch + custom CUDA | Current mature direct-feedback continuations | Four full 60-unit CPU/GPU comparisons at the tested parameters, plus checks at new starting contexts |

**PyTorch owns GPU arrays and matrix operations; custom CUDA computes mechanics and spatial geometry/polarity.** Spatial fields stay resident between steps. Fields, occupancy, and shell/contact products use float32; chemistry, transport, polarity, mechanical force arithmetic, and geometry accumulation use float64. Shell-squared adhesion buffers are float64. TF32 and fast-math are not enabled. Small diagnostics cross to the host routinely; complete spatial fields are copied for standard checkpoints. CPU/GPU agreement is measured with explicit tolerances, rather than assumed bitwise equality.

All four [full-horizon validation replays](docs/model.md#resident-gpu-backend-validation) passed against accepted CPU trajectories. Maximum chemical log difference was 4.92e-7 and normalized response difference was below 1.57e-6; target/network recovery times were unchanged. The gated adapter independently rechecks the saved evidence, physical parameters, device/software environment, source hashes, and starting state before scientific execution. GPU support currently covers mature, nondividing, polar, direct-feedback attribute states with conservative chemistry. Developmental cleavage, cue forcing, and no-feedback/nonpolar regimes require separate GPU work and validation.

A matched 120-step coupled benchmark measured **22.3× throughput** against the four-thread native CPU backend: **9.97 versus 222.37 ms/step**, including routine quality checks and sampled diagnostics. Checkpoint writes were measured separately and take about 1.4 seconds on both backends. This is a short mature-step benchmark, not a full developmental or multi-worker speedup. Earlier C/CUDA/PyTorch component measurements and execution transitions remain in [the model documentation](docs/model.md#opt-in-mechanics-optimization-and-validation).

The [seed-7 moving exchange-response study](docs/cell_response_moving.md#completed-moving-exchange-response-study) completed all fifteen runs with passing numerical quality. All eight exchanged-cell comparisons were closer to the chemical donor response than the destination reference. Six [targeted timestep-refinement continuations](docs/cell_response_moving.md#targeted-exchange-response-timestep-refinement) passed, preserving recovery times and nearest-reference classification. The [completed replication](docs/cell_response_moving.md#replication-across-developmental-histories) adds seeds 8 and 9: all six formation/retention controls and thirty response continuations passed numerical quality, and all sixteen new exchanged-cell comparisons favored the donor. This gives 24 comparisons within three developmental histories. Fresh and pre-relaxed exchanges reorganize concentrations; donor-nearer responses do not imply an unchanged transferred state or donor equivalence.

The [completed timestep study](docs/exchange_response_refinement.md) halves dt from 0.00375 to 0.001875 for six formation/retention and twelve same-state response continuations. All eighteen pass, after four full native CPU/GPU comparisons and ten starting-context checks pass. Donor-nearer classifications and recovery times are unchanged; maximum normalized waveform discrepancy is 0.0221%, below the 1% limit. Responses start from the original t=210 physical states, separately from formation refinement. These are interventions nested within existing histories, not new developmental replicas or evidence of spatial convergence.

The [completed neighbor-context assay](docs/neighbor_context.md) holds each target's initial chemistry and frozen geometry fixed while changing only its neighbors' initial chemistry. All sixty arm jobs pass independent DOP853/Radau checks. Eleven of 36 challenges shift late chemistry and 46 of 72 pulse comparisons shift response. Conservative averaging changes the low-state recipient's late concentrations in every fresh-exchange history. Two seed-8 negative pulses select different stable chemical endpoints from their own controls, whereas untouched equivalents recover.

The [completed delayed-pulse follow-up](docs/neighbor_context_delayed.md) finds that **all 24 pulses after settling recover**, including both percentage and original-amount dose controls, untouched references, and the alternative endpoints. The earlier nonrecoveries therefore reflect sensitivity during network reorganization; they do not show that the resulting equilibria fail to resist the tested pulses. This is a case-selected frozen-geometry follow-up within seed 8. Moving confirmation remains separate. Four CPU workers handle the small frozen graph; moving follow-ups retain the validated PyTorch/custom-CUDA architecture.

## Outputs

Current research experiments retain the following evidence, with exact locations/schema specified by each protocol:

| Artifact | Meaning |
|---|---|
| `protocol.json` | Frozen parameters, interventions, thresholds, source/input hashes, and interpretation limits |
| `status.json` | Root/child progress and failures; completion alone does not establish a biological claim |
| `initial_states.npz`, source checkpoint | Prepared chemistry, IDs, measured capacities, and the full physical start |
| Per-job `history.json` | Matched-time chemistry, volumes, centers, polarity, transport, and shape diagnostics; developmental histories additionally record attributes/context |
| Per-job `latest_state.npz` | Full standard attribute checkpoint with restart metadata, clock, lineage, and random streams |
| Per-job `result.json` | Protocol association and numerical-quality audit |
| `formation_comparison.json`, `comparison.json`, refinement reports | Stage-specific outcomes or aggregate assessment after required jobs finish |
| GPU-validation and prefix reports | Independent CPU/GPU discrepancies, endpoint comparisons, and evidence hashes |

Restore current checkpoints with the attribute class, preserving its constitutive law:

```python
from embryo.attribute_development import AttributeSimulation

sim = AttributeSimulation.restore(
    "outputs/exchange-response-histories/seed-8/formation/"
    "unexchanged_control/latest_state.npz"
)
for _ in range(100):
    sim.step()
sim.checkpoint("outputs/continued-attribute.npz")
```

This illustrates a manual CPU continuation, not a replacement for a hashed scientific protocol. Recorded GPU studies use the gated runner. Restoring an attribute checkpoint into the core `Simulation` would select the wrong material-response implementation.

The historical core batch CLI writes the following visualization and fate-model outputs:

| File | Contents |
|---|---|
| `viewer.html` | Offline, interactive 3D playback |
| `summary.png` | Shape snapshots, fate counts, and axis ratio |
| `config.json` | Full experiment parameters |
| `metrics.csv` | Volume errors, shape measures, and fate counts over time |
| `trajectory.json` | Sampled surfaces and activities for playback |
| `graph_preflight.json` | Discrete-mode stability on 4-, 8-, and 16-node reference graphs, computed before 3D stepping |
| `graph_history.json` | Actual weighted contacts, activities, polarity, and frozen-graph stability over time |
| `cleavage_spectra.json` | Before/after graphs, inherited-mode transfer, and signal partition jumps at abscission |
| `lineage.json` | Every cell's birth, parent, and division time |
| `final_state.npz` | Full final fields and RNG states for continuation |
| `diagnostics.json` | Run timing and numerical diagnostics |

The sampled trajectory is for visualization, not full-field reconstruction at earlier times. Core schema-3 checkpoints can be continued with their original class; older unsupported schemas are rejected. For a historical core checkpoint only:

```python
from embryo import Simulation
sim = Simulation.restore("outputs/my-run/final_state.npz")
for _ in range(100):
    sim.step()
sim.checkpoint("outputs/continued.npz")
```

## What is implemented

- Deformable 3D cell phase fields, soft volume constraints, exclusion, and interface attraction.
- Shape-aligned division, progressive equatorial constriction, mechanically gated abscission, conserved occupancy partitioning, and lineage tracking on CPU.
- Gierer–Meinhardt reactions, conservative geometry-dependent transport, mechanical dilution, and amount-balanced chemical partitioning.
- Finite-graph stability and cleavage spectral-transfer diagnostics; orientation averaging for polarity remains separate from chemical transport.
- Continuous activity-dependent material response and apical–basal polarity, without fate drift or A/B classification in the attribute model.
- Developmental ablations, frozen-context recovery/bistability, state exchange, matched pulse-response assays, and moving continuations.
- Opt-in native CPU acceleration and a validation-gated resident GPU backend for the supported mature regime.
- Per-step numerical audits, hashed protocols/evidence, observation histories, and resumable checkpoints.

The historical core/dashboard additionally retains its supplied fate switch for reproducing earlier studies. Current research remains capped at sixteen chemically well-mixed cells. A finer shape grid does not add chemical compartments or establish a many-cell continuum limit.

## What our current model and results can explain, and what remains missing

**The strongest current explanation is conditional collective chemical organization, with distinct mechanisms for formation and maintenance.** A chemical pattern, a persistent cell-associated state, a transferred response, and a developmental shape axis are different outcomes.

| Question | Current evidence | Limit of the claim |
|---|---|---|
| How does moving mechanics affect formation? | The matched component study finds suppression with polarity-only or full mechanics; tension plus adhesion without polar mechanics restores formation. Corresponding contact spectra move toward or away from uniform-state instability. | Necessity/sufficiency for suppression is demonstrated on one starting geometry, perturbation, coefficient set, and horizon. Polarity's chemical modulation remains active; geometry-only causation has not been isolated. |
| Can formed chemical differences survive movement? | Actual no-feedback developmental patterns survive feedback-on and feedback-off continuations from t=90 to 150 in all three planned histories, with the seed-7 timestep-halving check passing. | Finite-horizon maintenance after the sixteen-cell cap, not inherited identity or indefinite stability of the full moving system. |
| Are uniform and patterned states compatible with the same geometry? | Frozen-endpoint assays support locally stable uniform and patterned chemical equilibria on all six endpoint graphs from the three histories. | Chemical bistability conditional on a frozen graph, not a stability proof for mechanics and chemistry together. |
| Does chemical history affect subsequent behavior? | All 48 frozen exchanged-cell response comparisons and all 24 moving comparisons within three developmental histories are nearer the untouched donor response. Selected timestep checks retain classifications and recovery times. The neighbor-context assay shifts late chemistry in 11/36 challenges and responses in 46/72 comparisons, including state changes under conservative averaging. All 24 delayed pulses in the two selected nonrecovery cases return to their own settled endpoints under percentage/amount controls. | These are nested interventions. Nearest-reference behavior need not mean equivalence, unchanged transfer, or autonomy. Neighbor chemistry can change behavior even while the target's initial state and geometry are fixed. Immediate-pulse nonrecovery reflects reorganization-sensitive basin selection in these cases; moving-context confirmation remains necessary. |
| Can the system support differences without a clamped environment? | The separate finite, evolving-reservoir assay supports stable differences in selected exchange regimes with conservative cell–reservoir transfer. | That assay is not the live spatial embryo; its reservoir is not part of the current 3D/GPU state. Reactions still produce and remove chemicals. |
| Does the chemical loop explain a new shape axis? | Cleavage, packing, and polarity produce changing asymmetric geometry. Tested controls have not demonstrated additional persistent activator–inhibitor-specific elongation. | A spontaneously selected chemical shape axis, folding, or cavitation has not been established. |

See [completed mechanism and survival assessment](docs/feedback_completed_assessment.md), [frozen-endpoint replication](docs/feedback_endpoint_bistability.md), [response assays](docs/cell_response.md), and [moving responses](docs/cell_response_moving.md) for denominators, protocols, and numerical evidence.

The live operator conserves volume-weighted amounts under exchange and mechanical dilution. Manufactured regular/irregular transport benchmarks show approximately second-order convergence in their specified geometries. However, the overlap-to-area/path-length approximation fails general curved, separated, and nonorthogonal-contact closure screens. Conservative bookkeeping does not establish a continuum diffusion PDE on arbitrary embryo contacts. Historical full-developmental refinement also fails combined acceptance criteria; later response/backend checks do not retrospectively validate those trajectories. [Transport limitations](docs/geometric_transport.md), [developmental refinement](docs/development_refinement.md).

Cell identity is evaluated operationally through continuous attributes, persistence, perturbation responses, relocation, and dependence on surroundings. There is no imposed population classifier. Cell-autonomous memory, inheritance through later divisions, distinct biological function, and calibration to measured lineages remain missing. These require further experiments and ultimately biological evidence.

Order here means amplification and organization of differences in a structured mathematical system. Cells, regulatory laws, material couplings, and division rules are supplied assumptions. The project tests what follows from those assumptions; it does not derive the origin of life, metabolism, or the feedback laws themselves.

## Tests and next steps

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q

# Opt-in GPU checks require the compatible CUDA environment:
EMBRYO_CUDA_TESTS=1 CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q \
  tests/test_gpu_backend.py tests/test_gpu_response_runner.py \
  tests/test_validate_gpu_backend.py tests/test_exchange_response_histories.py
```

Software tests, numerical convergence, backend equivalence, and biological mechanism tests answer different questions. A faster backend is admitted only after whole-trajectory agreement, response/recovery comparisons, quality screens, and source/evidence verification; a backend pass is not biological validation.

The immediate sequence is:

1. Diagnose the [zero-contrast formation trajectory's numerical sensitivity](docs/polarity_refinement_assessment.md) before expanding the parameter map. The completed additional halving to dt=0.0009375 retains initiation in history 9, giving the same outcome at three timesteps, but the maximum adjacent-pair chemical discrepancy increases from 0.02574 to 0.03750 against the unchanged 0.01 tolerance. Backend, pilot, moving-quality and frozen endpoint checks pass; the original no-carry formation trajectory still fails quantitative agreement. A descriptive timing shift explains much of the transient difference without changing the failed gate. The [frozen chemical diagnosis](docs/chemistry_accuracy.md) now passes all 18 numerical cases on three measured contexts: second-order stepping, maximum log error 1.2e-6, and native/PyTorch agreement 1.1e-12. This rules out a large isolated chemical-stepping error in those fixed contexts. All eight [short moving precision controls](docs/geometry_precision.md) now complete with passing quality: preserving small phase updates reduces early transport disagreement by about 77% and endpoint phase disagreement by about 146-fold; the peak chemical discrepancy improves only 4% because it occurs at elapsed 0.15, while endpoint chemical disagreement falls about 87%. Float64 contact accumulation alone changes little. These six-unit tests do not resolve the long formation failure. The [completed full phase-carry formation pair](docs/phase_carry_formation.md#completed-formation-assessment) now passes all inherited gates through 240 elapsed units. Its raw maximum log chemical discrepancy is 0.000205, about 183-fold below the original no-carry discrepancy of 0.03750 and below the unchanged 0.01 limit. Both paths form persistent contrast and support local chemical bistability on their actual frozen endpoints. This is one existing mature-state history, not new developmental replication or general backend acceptance. The [completed three-timestep carry qualification](docs/phase_carry_convergence.md#completed-results) passes: adjacent maximum chemical discrepancies decrease from 0.000410 to 0.000205, a factor of 1.998, and all six short float64 mechanics comparisons pass with maximum field error 5.32e-8. This is a tested temporal trend in one mature history at zero directional contrast, not full-system or spatial convergence. The [completed matched carry polarity assessment](docs/phase_carry_polarity_assessment.md) qualifies all six full-window comparisons across histories 7, 8 and 9 at contrasts 0 and 0.35 and two timesteps. Polarity-dependent suppression is supported in history 9 only; histories 7 and 8 fail to initiate at chi=0 despite actual-input frozen formation. Their physical restriction remains unresolved, and both qualify for the specified dilution/contact-conductance controls. Two paths were reused and ten new jobs completed; no new histories were added. The [separate changing-geometry coupling replay](docs/geometry_chemistry_coupling.md) completes all 24 cases and independent references: beginning sampling is first order, midpoint sampling is second order and much more accurate. The local formation-window spacing check passes, but the full-start geometry-spacing difference 0.00365 fails its declared 0.001 screen. That failure remains explicit; shared-geometry timing error alone is much smaller than the original moving discrepancy. These add no independent histories. The [broader polarity study](docs/polarity_robustness_assessment.md) retains three shared histories: developed patterns survive all tested contrasts, and all longer frozen references form patterns. The [consolidated ledger](docs/results_ledger.md) remains an earlier snapshot. Moving confirmation of conservative/reset controls and pulse timing after the [neighbor screen](docs/neighbor_context.md) and [delayed assay](docs/neighbor_context_delayed.md) remains separate.
2. Test whether recovery or response survives stronger interventions, longer observation, and subsequent divisions. GPU division needs separate implementation and validation. Do not replace these tests with threshold labels.
3. Resolve geometric transport closure and full attribute-development spatial/time convergence before extending to substantially more cells or claiming a tissue continuum.
4. Isolate the chemical modulation of polarity and test shape-specific causality with matched controls; calibrate any biological interpretation independently.

Current protocols and decision rules are in [the implementation plan](docs/plan.md). Research sources and outputs are hashed before launch; numerical failures and uninformative reference pairs remain explicit outcomes. The root/child `status.json` files give live progress. Documentation describes the protocol; a complete aggregate result is produced only after all required trajectories pass their quality checks.

## Documentation guide

| Topic | Document |
|---|---|
| Consolidated conclusions, histories, validation decisions, and failures | [Results ledger](docs/results_ledger.md), [machine-readable snapshot](docs/results_ledger.json) |
| Supported compute paths and retired prototypes | [Backend policy](docs/backends.md), [backend archive](archive/backends/README.md) |
| Equations, parameter meanings, legacy alternatives, and backend method | [Model](docs/model.md) |
| Fate-free attributes, direct constitutive coupling, and developmental controls | [Attribute development](docs/attribute_development.md) |
| Conservative amounts, dilution, discrete stability, and closure limits | [Live transport](docs/live_transport.md) |
| Formation versus maintenance and polarity controls | [Completed feedback assessment](docs/feedback_completed_assessment.md), [polarity necessity control](docs/feedback_polarity_ablation.md) |
| Recovery and environment-supported memory | [Persistence](docs/attribute_persistence.md), [common environment](docs/attribute_common_environment.md), [finite reservoir](docs/attribute_finite_reservoir.md) |
| Frozen and moving exchange/response methodology | [Cell response](docs/cell_response.md), [moving response](docs/cell_response_moving.md) |
| Exchange formation and response timestep acceptance | [Current refinement protocol](docs/exchange_response_refinement.md) |
| Target state versus neighboring chemistry | [Neighbor-context assay](docs/neighbor_context.md) |
| Reorganization versus settled-state robustness | [Delayed-pulse follow-up](docs/neighbor_context_delayed.md) |
| Diffusivity ratio, geometry-conditioned basins, and moving polarity-tension checks | [Parameter robustness](docs/parameter_robustness.md), [completed assessment](docs/parameter_robustness_assessment.md) |
| Fixed-ratio polarity causality, longer horizons, and targeted timestep gates | [Polarity robustness protocol](docs/polarity_robustness.md), [completed assessment with retained failure](docs/polarity_robustness_assessment.md) |
| Third timestep for history-9 zero-contrast initiation | [Completed refinement and numerical diagnosis](docs/polarity_refinement_assessment.md) |
| Cleavage, shape, and numerical refinement | [Cytokinesis](docs/cytokinesis.md), [shape persistence](docs/shape_persistence.md), [developmental refinement](docs/development_refinement.md) |
| Discrete-to-continuum evidence and unresolved geometry | [Continuum bridge](docs/continuum_bridge.md), [geometric transport](docs/geometric_transport.md), [skew-boundary correction](docs/skew_boundary_correction.md) |
| Historical visual interface and earliest results | [Dashboard](docs/dashboard.md), [historical results](docs/results.md) |
| Priorities and acceptance criteria | [Plan](docs/plan.md) |

## Manuscript draft

The [O1-focused manuscript](manuscript/introduction_methods.md), [compiled PDF](output/pdf/embryogenesis_o1.pdf), and [build/evidence notes](manuscript/README.md) are research drafts, not a submitted paper. Use the model and assay documents above for the latest methodology and execution status; manuscript claims must be tied to the specific completed evidence included in each draft.
