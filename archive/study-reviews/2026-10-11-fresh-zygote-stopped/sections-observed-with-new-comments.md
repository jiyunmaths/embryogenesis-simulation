# Abstract {.unnumbered}

Starting a chemical pattern, maintaining it, and preserving a cell's behavior after its surroundings change are different dynamical questions. We investigate them in a three-dimensional model of deformable, dividing cells with activator-inhibitor reactions, conservative transport, and polarity-dependent mechanics. No cell identities or number of types are supplied. In one matched component geometry, polar mechanics suppresses initiation; removing its mechanical contribution restores formation. However, enabling full coupling after a pattern develops preserves chemical contrast in all three tested developmental histories. Uniform and patterned chemical equilibria coexist on both frozen endpoint graphs from each history. Chemical-state exchange produces donor-nearer pulse responses in these same three histories, without retaining the transferred concentrations unchanged. Amount-preserving neighbor averaging shifts the fresh-exchange low-state recipient's late chemistry in every history, demonstrating support by the collective network. In a follow-up selected from two immediate-pulse nonrecoveries in one history, percentage-matched and amount-matched pulses applied after settling return to their own equilibria. These findings separate initiation from maintenance and show state-sensitive basin selection during reorganization alongside robustness after settling. Cells, branches, graphs, and pulses are nested interventions; the replication unit is the developmental history. The results support collective chemical differentiation and network-supported context sensitivity, not autonomous or inherited biological identity. General geometric transport accuracy and full developmental convergence remain unresolved.

# Introduction

An embryo begins as one cell, but its descendants develop different properties and responses within an organized tissue. These differences arise while cells divide, change shape, and rearrange their contacts. Understanding development therefore requires more than an explanation of how a pattern first appears: it also requires an explanation of how differences persist as their physical and chemical surroundings change. Here we investigate how chemical feedback, communication between cells, and cell mechanics contribute to these distinct stages of organization.

Activator-inhibitor feedback offers a mechanism for generating differences without assigning cell types in advance. An activator promotes its own production and stimulates an inhibitor, which in turn suppresses the activator's activity. If inhibition spreads more widely than activation, a local increase can reinforce itself while limiting activation in surrounding regions. With suitable reaction and exchange rates, small fluctuations can therefore grow into a spatial pattern [@turing1952; @gierer1972]. Experiments support these ingredients: Lefty has a higher effective diffusivity than Nodal in zebrafish [@muller2012], an engineered Nodal-Lefty circuit produces signaling domains in genetically identical mammalian cells [@sekine2018], and a related reaction-diffusion network [what does the related reaction-diffusion network mean? please make it clear what the model is] contributes to vertebrate digit patterning [@raspopovic2014]. These examples motivate asking how local self-amplification and more widespread inhibition can organize a group of initially similar cells.

In a continuous description, chemical concentrations vary smoothly throughout tissue. When each cell is treated as a well-mixed compartment, with one concentration of each chemical, communication is instead described by a network of cells and their exchange of chemicals [@othmer1971; @nakao2010]. Small chemical differences across this network can be described as combinations of basic arrangements of high and low concentrations, called *spatial modes* [is there a citation for this statement/definition? what does the basic arrangement mean? overall, this sentence does not make sense to me]. The network has a finite set of these modes, each with a different exchange strength. An *eigenvector* lists the relative concentration differences between cells in one mode; its *eigenvalue* measures how strongly exchange tends to smooth those differences. Activator-inhibitor feedback and exchange together determine whether each mode grows or fades [@wolfrum2012; @hajali2025]. 

As cells deform and rearrange, their contact areas, separations, and volumes change, altering chemical exchange. Our hypothesis is that this remodeling can leave none of the network's modes able to grow [what does 'leave none of the network's modes able to grow' mean? say it in plain language], closing the opportunity for small fluctuations to develop into a pattern. A calculation with contacts and volumes held fixed describes only that particular network. We test the hypothesis by following chemistry and cell mechanics together and intervening directly on chemical exchange.

The route of communication matters for this hypothesis. [what is the animal?] Connexin and gap-junction interventions alter periodic feather patterning, providing experimental motivation for contact-mediated communication [@tseng2024]. In human gastruloids, Nodal activity spreads through induction of Nodal production in neighboring cells, demonstrating a relay rather than unrestricted ligand diffusion [@liu2022]. We represent chemical feedback with an activator and an inhibitor in each cell. Both chemicals are exchanged with neighboring cells only across contacts. The *conductance* of a contact describes how easily chemicals pass between the two cells. The conductance is the same in either direction, and the transport conserves the total amount of each chemical. 

Mechanics can change the geometry through which communication occurs. In early mouse embryos, contractility differences help sort cells into inner and outer positions [@maitre2016], and apical-domain manipulation shows that polarity can influence the first lineage decision [@korotkevich2017]. Fluid pressure reorganizes cell contacts during cavity formation [@dumortier2019], while receptor localization and secreted inhibitors connect tissue geometry to signaling boundaries in mammalian gastruloids [@etoc2016]. Genetic and optogenetic interventions in zebrafish show that adhesion and tissue rigidity regulate Nodal dispersal and Lefty timing [@autorino2026]. That mechanism acts through extracellular porosity, which differs from the contact conductances defined above; it nevertheless demonstrates that tissue architecture can regulate chemical transport. These findings motivate our investigation of whether mechanical remodeling changes the likelihood to initiate chemical differences among cells and whether it also disrupts differences that already exist.

Existing models establish important parts of this connection [be specific about the connection]. Mercker and colleagues couple morphogen regulation to curvature in a continuous tissue sheet [@mercker2013]. Yu and colleagues connect reaction-diffusion signaling to three-dimensional vertex mechanics and cell rearrangements, including network reconnection [@yu2024]. Polarity-based models (what are the polarity based models?) generate collective tissue shapes [@nissen2018], and phase-field frameworks represent deforming cell boundaries and division [@kuang2023]. Mechanochemical coupling and changing cell neighborhoods have therefore already been established modeling strategies. The gap addressed in our study is more specific: whether remodeling of chemical exchange conductances affects initiation and maintenance differently, and whether preserving those conductances changes the outcome while the cells continue to deform.

Theory also treats pattern formation on networks that vary in time and on evolving domains [@petit2017; @vangorder2021]. In those analyses, network or domain evolution is specified rather than generated by the coupled cell mechanics studied here. Such theory cautions against treating a sequence of frozen spectra [what are the frozen spectra?] as the stability analysis of a time-dependent system. We therefore combine spectral diagnostics with moving simulations and direct interventions on transport. This allows the proposed conductance mechanism to be tested rather than inferred solely from a correlation between geometry and chemical contrast.

Initiation and maintenance need separate tests because uniform and patterned chemical states can both be stable on the same network [@nakao2010; @wolfrum2012]. When remodeling removes all growing modes, small chemical fluctuations decay on the resulting frozen network. Yet a developed pattern with substantial concentration differences may remain stable. Remodeling can therefore close the window for onset without erasing established differences. This makes maintenance after loss of the growing modes a prediction to test. Experiments in a nonliving reaction-diffusion medium demonstrate different thresholds for the appearance and disappearance of chemical patterns [@kaminaga2006]. A biological analogy comes from intestinal organoids, where lumen inflation before crypt budding impairs budding, whereas inflation after budding preserves the structure [@xue2025]. This morphological result does not establish chemical bistability in embryos, but it shows why the timing of a mechanical intervention can matter.

A further distinction concerns what differences between cells mean. We use *chemical state* for a cell's current activator and inhibitor activities, and *cell identity* for a recognizable, relatively persistent profile of attributes and responses under stated conditions. An instantaneous concentration difference does not establish identity. Signaling history can change later behavior: WNT pre-exposure changes ACTIVIN-dependent differentiation in human stem cells [@yoney2018], epiblast stem-cell populations can retain response biases after reseeding under common differentiation conditions [@menezes2025], and time-integrated BMP signaling predicts differentiation better than an instantaneous signaling level in a human stem-cell model [@teague2024]. These observations motivate measuring responses and their history alongside concentrations.

Persistence may still require support from other cells. Neighbor interactions sustain differentiation-related signaling in small human stem-cell colonies [@nemashkalo2017]. In Xenopus transplantation experiments, early muscle precursors require neighboring cells to differentiate, whereas later precursors can differentiate singly [@katogurdon1993]. Intestinal regeneration and organoid experiments connect mechanical context, lineage-biased chromatin states, and FOXA1/Delta-Notch feedback to persistent spatial fate patterns [@schwayer2026]. These biological mechanisms include molecular and epigenetic machinery absent from our reduced model. A lasting cell-associated response can therefore be informative without demonstrating autonomous, committed, or inherited identity.

We examine these questions in a three-dimensional model of deformable, dividing cells. Each cell contains an activator and an inhibitor. Conservative exchange uses conductances estimated from current contacts and accounts for changing cell volume. Interface costs, volume-restoring forces, overlap repulsion, and baseline adhesion act alongside activity-dependent tension and adhesion. Apical-basal polarity develops from surface exposure and makes mechanics directional; it can arise even when chemistry is nearly uniform. Setting the directional-tension contrast $\chi$ to zero removes that directional action, while the other mechanical terms and activity-dependent modulation remain active. No cell types, number of types, or fate switches are supplied. The model lets us ask:

1. Does mechanically remodeled contact transport restrict initiation, and does preserving the starting conductances restore it?
2. Can an established pattern persist after that remodeling removes the modes supporting initiation?
3. After chemical states are exchanged between cells in different geometric settings, how do the transferred state and surrounding chemistry influence later responses?

We use matched mechanical controls, moving continuations, frozen-network stability tests, chemical-state exchanges, and neighborhood and pulse interventions. The transport intervention holds the starting symmetric conductances fixed while retaining evolving mechanics, actual cell volumes, and dilution; it freezes neither the tissue nor the per-volume chemical operator. Polarity controls remove directional tension while retaining polarity dynamics and activity-dependent tension and adhesion. Chemical-state exchange transfers the two regulators without relocating cells, and each pulse is compared with its own unperturbed continuation.

The numerically qualified mature-tissue continuations show why the conductance mechanism deserves this focus. Onset is tested by resetting chemistry to near-uniform activities on $t=150$ geometries previously evolved with patterned chemistry and mechanical feedback. Because this reset changes tension and adhesion coefficients on the retained geometry, it may induce a relaxation transient. With evolving conductances, polarity-specific suppression occurs only in history 9: near-uniform chemistry forms a pattern without directional polarity tension but fails with it. Histories 7 and 8 fail to initiate at either setting. Preserving the starting conductances restores formation in all three histories at both settings while mechanics and dilution remain active. Frozen spectral diagnostics are consistent with loss of the supported unstable modes in evolving branches and their retention in conductance-preserving branches. The mature initiation series lacks a matched continuation with all mechanical feedback disabled, so it cannot distinguish baseline mechanical relaxation from activity-dependent tension/adhesion as the source of remodeling at $\chi=0$. The intervention supports a conductance-mediated contribution to failed initiation in these reset preparations, without establishing its force-level cause or full moving-system stability.

Developed patterns behave differently. They retain chemical contrast for 240 model-time units under evolving mechanics in all three histories at both polarity-tension settings, and their frozen endpoint networks support local uniform/patterned coexistence. Chemical-state exchange produces donor-nearer pulse responses despite subsequent concentration reorganization. Changing surrounding chemical allocation at matched initial species amounts produces immediate and selective delayed response differences. The history is the replication unit: branches, cells, endpoint graphs, and perturbations are nested observations within histories 7, 8, and 9, not additional developmental replicas.

The contribution is a controlled account of how remodeling during continued mechanics can limit the formation of chemical organization through changing communication while allowing developed differences to persist. It also shows why persistence after chemical transfer does not imply independence from the surrounding network. These mature reset and continuation assays do not establish that the same remodeling occurs during uninterrupted development from a fresh zygote. The conclusions concern collective chemical differentiation and context-sensitive behavior. Autonomous or inherited biological identity, a chemically generated tissue axis, calibration to a living embryo, general geometric transport accuracy, and full developmental convergence remain open.

# Methods

## Model overview

The simulation begins with one nearly spherical zygote and follows division to sixteen cells in three-dimensional space. Each cell has a deformable boundary, two well-mixed chemical concentrations, and a polarity vector. “Well mixed” means that each cell has one activator value and one inhibitor value; we do not resolve chemical gradients inside an individual cell. The cell boundaries and the aggregate outline can change throughout the simulation.

| Quantity | Meaning in the model |
|---|---|
| $\phi_i(\mathbf x,t)$ | Smooth description of cell $i$: about one inside and zero outside |
| $a_i(t),b_i(t)$ | Activator and inhibitor concentrations in that cell |
| $V_i(t),V_i^\star$ | Measured volume and prescribed target volume |
| $\mathbf p_i(t)$ | Cell-level apical direction and polarity strength, with $0\le\lVert\mathbf p_i\rVert\le1$ |
| $g_{ij}(t)$ | Geometric conductance governing exchange between cells $i$ and $j$ |

A dot denotes a time derivative; $i,j$ label cells, and $\mathbf x$ is a position in the box. All quantities are dimensionless and are not calibrated to a particular organism. We record log chemical concentrations, polarity magnitude, and cell-shape descriptors as separate attributes. Position, surface exposure, volume, and lineage describe context. Mechanical coefficients derived directly from activity are not counted as independent evidence of differentiation.

The founder starts with $a=b=1$ and zero polarity. Division times, orientational choices in nearly spherical cells, and small chemical partition perturbations introduce variability without a prescribed embryo axis. Total target volume remains fixed and is divided among descendants; there is no growth between divisions. Once the sixteen-cell cap is reached, cells continue to deform and exchange chemicals but no longer divide. The old fate variable is disabled and cannot affect these experiments.

## Operational definitions: state, phenotype, identity, and context

A **cell** is a deformable compartment with its own phase field, well-mixed regulators, and polarity. A **cell ID** tracks that individual; a parent ID records ancestry. Neither is a cell-type label. The full instantaneous **cell state** includes chemical activities, polarity, shape, and division bookkeeping. For chemical assays, **chemical state** means only the positive pair $(a_i,b_i)$. Concentrations are distinct from amounts $(V_i a_i,V_i b_i)$ and from any specific molecular identity program.

An **attribute** is a measured property. The developmental measurements form the continuous profile

$$
\mathbf z_i(t)=\bigl(\log a_i,\log b_i,\lVert\mathbf p_i\rVert,\rho_i-1,B_i\bigr).
$$

Here $\rho_i$ is the ratio of the largest to smallest principal cell lengths, obtained from the occupancy-weighted spatial covariance; $B_i$ is asphericity. For covariance eigenvalues $\mu_{i1}\le\mu_{i2}\le\mu_{i3}$ and their mean $\bar\mu_i$, the definitions are

$$
\rho_i=\sqrt{\mu_{i3}/\mu_{i1}},\qquad
B_i=\frac{3}{2}\frac{\sum_{m=1}^3(\mu_{im}-\bar\mu_i)^2}{(\sum_{m=1}^3\mu_{im})^2}.
$$

The axis ratio is one and asphericity zero when the second moments are isotropic; neither descriptor alone proves a spherical boundary. These cell-scale measurements are distinct from the aggregate axis ratio. Numerical evaluation protects against division by a vanishing smallest eigenvalue. **Phenotype** here means the observable attribute profile and behavior under specified conditions. **Response behavior** is the time course after a defined challenge, relative to the same background's unperturbed continuation; its exact pulse definition appears below. Tension and adhesion are computed from activity and do not supply independent evidence for identity.

We operationalize **cell identity** as a distinguishable profile that persists over a declared interval and has reproducible responses to stated challenges. This is a research definition rather than a claim that the current reduced variables specify biological cell types. The experiments test persistence, local recovery, exchange, and context dependence separately. A high-activity cell, a colored rendering, a statistical cluster, or a contrast threshold alone is insufficient. **Differentiation** means development of lasting differences between initially similar cells; we use **chemical differentiation** for differences in the reduced regulators and reserve the broader identity interpretation for evidence about attributes and responses. No discrete type count is inferred from a continuous spread of states.

**Geometric context** comprises position, volume, surface exposure, and contacts. **Chemical context** comprises neighboring activities and exchange conditions. Lineage and developmental history are covariates, not extracellular signals. This state/context distinction follows each intervention: chemical exchange changes only the two regulators while holding destination shape, polarity, position, and lineage initially fixed. Those fixed quantities are part of the full cell state, even though they provide the background for a chemical-state assay. A phenotype may remain cell-associated while depending on collective network support.

| Stronger property | Definition and evidence needed |
|---|---|
| Fate | A future developmental outcome or potential, rather than an instantaneous marker value. No organism-specific lineage outcomes are specified here. |
| Commitment | Retention of a developmental outcome after a stated change in instructive conditions. Frozen chemical recovery alone does not demonstrate it. |
| Memory | A lasting effect of prior conditions after the initiating intervention ends. Memory can reside in a collective state; an extra memory variable is not required by this definition. |
| Cell autonomy | Retention under interventions that remove or standardize environmental support. Maintenance on the original graph is insufficient to establish autonomy. |
| Inheritance | Transmission of an attribute or response profile through subsequent division. Parent bookkeeping or an assumed partition rule is insufficient evidence; the mature persistence assays do not test it. |

For example, a high-activity cell that recovers after a pulse while its neighbors remain available demonstrates conditional persistence and recovery. It may be part of a network-supported phenotype. That observation alone establishes neither an autonomous identity nor an inherited fate.

## Operational definitions: organization and symmetry breaking

**Initiation** or **formation** means development of contrast from a specified near-uniform state with small perturbations. **Maintenance** or **persistence** means retention of an already developed difference over a declared window. A restart on mature geometry and a fresh zygote trajectory answer different initiation questions. Recovery, context independence, and inheritance are additional tests, not synonyms for persistence.

An **attractor** is a state or set approached by the dynamics; its **basin** is the set of initial conditions that approach it. **Local bistability** here means locally stable uniform and patterned chemical equilibria on the same frozen graph. This concerns the collective chemical network, not two autonomous fates inside every cell or stability of the full moving system.

**Signaling symmetry breaking** denotes development of differences from approximately uniform chemical activities. Across-cell SD of log activator measures contrast; it does not alone establish spatial domains, axis selection, or a Turing mechanism. **Geometry or shape symmetry breaking** denotes loss of a specified geometric symmetry, such as approximate rotational symmetry of an aggregate. Local cell shape, contact-network structure, and aggregate shape are distinct geometric observables. **Emergent organization** is produced by the dynamics rather than imposed as a final form. Calling an axis **spontaneous** additionally requires excluding a prescribed directional cue and testing inherited geometry, boundaries, and grid bias. Chemical differences on an asymmetric graph do not by themselves establish spontaneous geometric symmetry breaking.

## Chemical feedback: local amplification and inhibition

Within a cell, the reduced Gierer-Meinhardt reactions are

$$
R_a(a,b)=\frac{a^2}{b}-a,\qquad
R_b(a,b)=\beta(a^2-b).
$$

The first term in $R_a$ makes activator production increase with activator itself and decrease with inhibitor. The term $-a$ removes activator. In the second equation, $\beta a^2$ produces inhibitor in response to activation, and $-\beta b$ removes it. The parameter $\beta=2$ sets the inhibitor reaction timescale relative to activator turnover. The positive uniform reference is $(a,b)=(1,1)$.

These are deliberately simple regulatory laws. The quadratic term is an assumed nonlinearity, not a measured molecular stoichiometry. The equations do not track nutrients or energy, and their reactions can create or consume regulator amount. Thus “self-maintained pattern” refers to the dynamics of this supplied reaction system, not metabolic self-sufficiency.

## Communication between cells: conserve amount, account for volume

Chemical exchange depends on the current cell geometry. A larger estimated contact area increases exchange; greater separation between cell centers reduces it. These effects are combined in a symmetric conductance $g_{ij}=g_{ji}$. Its construction from diffuse cell interfaces is given in the appendix.

For either chemical, denoted $c$, the amount balance is

$$
\underbrace{\frac{\mathrm d(V_i c_i)}{\mathrm dt}}_{\text{change in amount}}
=
\underbrace{V_iR_c}_{\text{local reaction}}
+
\underbrace{D_c\sum_jg_{ij}(c_j-c_i)}_{\text{exchange with neighbors}}.
$$

$D_c$ sets the transport scale: $D_a=0.02$ for activator and $D_b=0.4$ for inhibitor. The inhibitor therefore has the larger transport coefficient. On each link, the amount gained by one cell is lost by the other. If the receiving cell has twice the volume, the same incoming amount produces half the concentration change.

Expanding the amount derivative gives the concentration equation used to interpret changing cell volumes:

$$
\dot c_i=R_c+\frac{D_c}{V_i}\sum_jg_{ij}(c_j-c_i)
-c_i\frac{\dot V_i}{V_i}.
$$

The last term is dilution or concentration by deformation. For example, expansion lowers concentration if chemical amount does not change. Numerically, this effect is applied using measured volumes before and after each mechanical step. Exchange and dilution conserve amount; the reactions do not.

For spectral analysis, the exchange rule is written as a volume-weighted graph operator,

$$
\Delta_V=-M^{-1}K,\qquad
M=\operatorname{diag}(V_i),\qquad
K=\operatorname{diag}(G\mathbf1)-G,
$$

where $G$ contains the conductances. This conservative construction follows an amount-balance principle used in finite-volume methods [@eymard2000]. It replaces normalized random-walk transport for chemistry. Nevertheless, the geometric conductances remain an approximation: conserving amount does not prove that they reproduce continuum diffusion on arbitrary deformed cell contacts.

## Cell shape, tension, adhesion, and polarity

A phase field avoids treating each cell as a rigid sphere. We measure its volume using the smooth occupancy function $H$, and distinguish bulk cell interior from exterior using the double-well function $Q$:

$$
H(\phi)=\phi^2(3-2\phi),\qquad
Q(\phi)=\phi^2(1-\phi)^2,\qquad
V_i=\int_\Omega H(\phi_i)\,\mathrm d\mathbf x.
$$

$H$ equals zero outside and one inside, with zero slope at both limits. $Q$ has minima at those same limits and penalizes intermediate values. A gradient penalty prevents an infinitely sharp boundary. The mechanical cost also penalizes departures from target volume and overlap between cell interiors, while adhesion favors contact between interfaces. The complete cost and the resulting force are provided in the appendix, with each term identified.

Chemical activity changes baseline tension and adhesion through

$$
r_i=\tanh(a_i-1),\qquad
\gamma_i^0=\gamma_0(1+c_\gamma r_i),\qquad
A_{ij}=A_0(1+c_A r_i r_j).
$$

Here $\gamma_0=1$ and $A_0=4$ are the baseline tension and attraction strengths. The response $r_i$ is zero at unit activator and varies smoothly; it is not an identity threshold. The contrasts $c_\gamma=0.25$ and $c_A=0.35$ determine how strongly chemistry changes those material properties. Same-sign responses increase attraction under this chosen law. No separate memory variable is introduced.

Polarity supplies a distinct mechanical route.

**Polarity** means directional asymmetry in a cell's organization. In this reduced model, $\mathbf p_i$ points toward the modeled **apical** side; **basal** denotes the opposing direction. These are local cell directions, not an embryo-wide head-tail axis. The magnitude $\lVert\mathbf p_i\rVert\in[0,1]$ measures polarity strength, with no defined orientation at zero; the numerical update projects the vector into the unit ball. Its value is a normalized model measure, not a molecular concentration, probability, or apical surface fraction. Exposed cortex is boundary unoccupied by neighboring cells, whose outward normals supply the orientation cue. This geometric abstraction is motivated by contact-dependent apical-domain formation [@korotkevich2017]. It does not explicitly resolve apical molecular machinery, a basal domain, or planar polarity.

A polarity vector distinguishes an oriented side, whereas a long shape axis describes geometric extent. A cell may be polar without elongation, elongated without modeled polarity, or polar while chemically similar to its neighbors. Thus polarity, morphology, chemical differentiation, and identity are separate measurements. The coefficient $\chi$ controls polarity's mechanical action; $\chi=0$ leaves polarity dynamics active while removing its directional effect on tension.

The vector is driven by the exposed-cortex cue $\mathbf q_i$, aligned with neighboring polarities, and limited by decay and saturation:

$$
\dot{\mathbf p}_i=
\frac{2a_i}{1+a_i}\mathbf q_i
+0.25(\Delta_{\mathrm{rw}}\mathbf p)_i
-(0.5+|\mathbf p_i|^2)\mathbf p_i.
$$

The normalized operator $\Delta_{\mathrm{rw}}$ is used only to average neighboring orientations, not to transport chemicals. The appendix defines the exposure cue and directional tension field. Polarity reduces tension toward its oriented side with strength $\chi=0.35$ under full coupling.

At uniform chemistry, $2a/(1+a)=1$, so geometric exposure can drive polarity without a pre-existing chemical pattern. In contrast, $r_i$ is then zero, and the activity-dependent material changes vanish. This distinction motivates the component controls. Because activity still modulates polarity away from the uniform reference, a polarity-only mechanical control is not a test of complete independence from chemical regulation.

## The experimental comparisons

We distinguish a pattern's **formation** from its **maintenance**, and a fresh developmental history from a mature-state intervention. The completed comparisons are:

| Experiment | Starting state and intervention |
|---|---|
| Developmental comparison | Same seed-7 zygote; full mechanical coupling versus constant-material mechanics through time 90 |
| Mechanical components | Same mature time-18 geometry; six coupling choices, each with near-uniform or prepared patterned chemistry; continue to time 78 |
| Developed-pattern survival | Clone an actual no-feedback time-90 state; enable full coupling in one copy, leave the other unchanged; continue both to time 150 |
| Endpoint attractors | Freeze each time-150 geometry from seeds 7, 8, and 9; perturb either uniform or developed chemistry and test its return |
| Numerical and history checks | Halve the seed-7 survival timestep; repeat development and paired survival with seeds 8 and 9 |
| Exchange and response | Transfer chemistry on mature endpoint graphs; compare destination and donor responses under frozen and moving geometry |
| Neighbor challenges | Keep target initial chemistry and geometry fixed; change only neighbors once, then allow the whole chemical network to evolve freely |
| Delayed pulses | Revisit the two immediate-pulse nonrecoveries after the challenged controls and alternative endpoints have settled |

The six completed component choices are defined by three mechanical contrasts:

| Mechanical arm | $c_\gamma$ | $c_A$ | $\chi$ |
|---|---:|---:|---:|
| Baseline | 0 | 0 | 0 |
| Tension only | 0.25 | 0 | 0 |
| Adhesion only | 0 | 0.35 | 0 |
| Tension plus adhesion | 0.25 | 0.35 | 0 |
| Polarity only | 0 | 0 | 0.35 |
| Full coupling | 0.25 | 0.35 | 0.35 |

Chemistry and polarity dynamics remain active in all six arms. The interventions change which variables affect mechanics. Each component pair starts with exactly the same geometry, volumes, lineage, and polarity. Formation chemistry is perturbed around one with seed 7 and volume-weighted RMS amplitude 0.001. Persistence chemistry is obtained by allowing that perturbation to relax for 240 chemical time units on the frozen starting graph. It is then placed back on the original mechanical state. Reactions can change mean amounts during preparation, so these two starting families differ in more than spatial contrast alone. Within either family, all mechanical arms have identical starting chemistry.

The developed-pattern survival assay avoids that preparation step. It uses chemistry produced during actual no-feedback development to time 90. Enabling full coupling changes only the feedback flag and experimental mode; no fields, concentrations, polarity, timers, or random states are reset. Both copies continue with moving geometry and transport. The assay has completed for all three planned developmental seeds, 7, 8, and 9. These are independent initial random streams, not independent species or biological samples.

The sixth component arm retains tension and adhesion but sets $\chi=0$, testing whether removing polar mechanics restores formation. Both original chemical starting families have completed this intervention. Polarity dynamics remain active. A separate intervention fixing the polarity activation factor at one would be needed to establish whether chemical modulation of polarity is required; that test has not been performed.

## Replication units and conditional starting states

A developmental history is a zygote simulation with its specified initial random stream. The current evidence uses histories 7, 8, and 9. Cloning a checkpoint, freezing a graph, changing a parameter, targeting another cell, adding a pulse, or halving the timestep does not create a new history. Paired branches within a history share their starting state. Histories are distinct simulation realizations, not biological samples, and the three-history pilot does not support a precise population success probability. We do not pool nested measurements into a binomial sample size or combine repeated appearances of the same histories across assays.

| Question | Histories represented | Nested measurements and qualification |
|---|---:|---|
| Zygote comparison and component suppression | 1 | Paired development and twelve mature-geometry component continuations |
| Maintenance after enabling feedback | 3 | Two moving branches per history |
| Endpoint chemical coexistence | 3, reused | Two graphs and forty-two chemical starts per graph |
| Frozen exchange-response comparisons | 3, reused | Six endpoint graphs; forty-eight exchanged-cell reference comparisons |
| Moving exchange-response comparisons | 3, reused | Eight reference comparisons per history |
| Neighbor-dependent state and response | 3, reused | Twelve challenges and twenty-four pulse comparisons per history |
| Delayed-pulse follow-up | 1, selected | Two cases, six endpoint jobs, twenty-four pulses |

Moving exchange starts with the patterned chemical equilibrium on each mature time-150 geometry, rather than the raw developmental chemistry. In each history, the initially minimum- and maximum-activator cells are selected before observing responses. Untouched, freshly exchanged, and frozen pre-relaxed exchange chemistry then coevolve with mechanics to time 210. Pulse/control comparisons run from those physical endpoints to time 270. These are tests of prepared mature basins, not additional zygote-to-identity development. The labels "low" and "high" identify starting concentrations only.

## Context interventions and matched pulse measurements

For a chemical species $c$, exchanging cells $i,j$ with unequal volumes uses $c'_i=sc_j$, $c'_j=sc_i$, where $s=(V_ic_i+V_jc_j)/(V_ic_j+V_jc_i)$. This preserves the pair's amount; the factor is computed separately for each species. Shape, polarity, lineage, and other cells' initial chemistry stay at their destinations. It is a chemical transplant, not physical relocation. Frozen assays also include exact concentration exchange and a sham. A separate exploratory one-history screen tests all 120 pairs; only 96 have sufficient initial contrast for its state classification.

An activator pulse changes only the selected target by a factor $f$. Relative pulses add or remove external chemical amount. Compare every pulse with its own unperturbed continuation, which removes the background's ordinary drift. The signed target response is $\log(c_i^{\mathrm{pulse}}/c_i^{\mathrm{control}})/|\log f|$, for each species. Response distances are time RMS differences over a matched horizon. A response nearer the original donor than the destination is a descriptive comparison, not equivalence of responses, proof of unmodified chemical transfer, or autonomous identity. Reference waveforms must be separated by more than 0.01 for the nearest-reference comparison to be informative.

The neighbor assay uses the original untouched and fresh-exchange time-210 starts, with two targets in each of three histories. Direct neighbors have positive off-diagonal $\Delta_{ij}$. Arms retain the original state, copy it back as a sham, replace neighbor concentrations by their volume-weighted mean, or reset both species of every neighbor to an initial low- or high-activator cell's values. Mixing conserves each species' total amount. Resets change recorded amounts externally. Target and nonneighbor initial chemistry are identical across arms. Interventions occur once; afterward all cells evolve freely on the original frozen graph, without clamping or a reservoir.

Each arm has a control and immediate target pulses $f=0.9,1.1$, followed for 240 chemical time units. Neighbor challenges require log-RMS size at least 0.05 and initial target-pair separation above 0.1. A late state shift means target/control log-RMS, maximized over the final 24 units and divided by the original target-pair separation, exceeds 0.1. Response shifts exceed time RMS 0.01 relative to the same target's untouched waveform. These predeclared thresholds describe changes, not identities. Recovery requires all later samples to stay within 10% of the initial displacement and at least 24 further units to be observed. Nonrecovery is censored at the observation horizon.

The delayed follow-up selects all and only the two immediate-pulse nonrecoveries. For each, copy the full sixteen-cell untouched-control, challenged-control, or alternative negative-pulse endpoint after 240 frozen chemical units. Do not reapply the neighbor intervention or restore earlier concentrations. On each network, compare a control with current $\pm10\%$ pulses and original-amount-matched pulses for another 240 units. For original target activator $a_i^0$ and settled value $a_i^*$, the amount-matched factor is $1+(f-1)a_i^0/a_i^*$. Recovery uses this actual factor. The frozen settling interval is not additional moving developmental age. This planned case follow-up tests perturbation timing and dose, not cross-history replication of switching frequency.

Isolation and common-environment assays provide supporting controls in one history. Removing transport tests reaction-only persistence. Clamped reservoirs replace cell-specific surroundings by common exchange conditions; finite-reservoir assays additionally let the shared reservoir respond and include both formation and release of developed states. These controls examine whether support can come from evolving chemical surroundings rather than an imposed fixed composition. Their numerical exceptions remain reported, and reservoir initial states or perturbation seeds do not count as new developmental histories.

## What counts as a pattern, and what counts as persistence?

Our primary measure of chemical contrast is the standard deviation of log activator across cells:

$$
s_a(t)=\sqrt{\frac1N\sum_i[\log a_i(t)-\overline{\log a}(t)]^2}.
$$

A value near zero indicates similar concentrations; a larger value indicates greater differences. Logs measure proportional differences and prevent the interpretation from depending on an arbitrary common concentration scale. We require $s_a>0.1$ throughout the specified final observation window. For the component study this window is time 63--78; for moving survival it is time 135--150. Developmental descriptive statistics use time 75--90. A threshold crossing is an operational measure of contrast, not a cell-type assignment or a proof of equilibrium.

For moving survival, we also ask whether the differences remain associated with the same cells. We correlate each cell's current log activator with its value at time 90, using fixed initial-volume weights. A minimum correlation of 0.8 throughout the final window is required. An undefined correlation caused by vanishing contrast does not pass. High correlation does not require identical concentrations and does not demonstrate inheritance, because the assay begins after cleavage has stopped. The exact correlation and chemical error measures appear in the appendix.

## Stability tests: can differences grow, and can developed states recover?

An instantaneous graph spectrum answers a restricted question: if this geometry were held fixed and chemistry were nearly uniform, which small spatial perturbations would grow? With homogeneous reaction Jacobian $J_0$, each transport mode is tested using

$$
J_0=\begin{pmatrix}1&-1\\2\beta&-\beta\end{pmatrix},
\qquad
\sigma_k=\max\operatorname{Re}\operatorname{eig}
\left[J_0-\lambda_k\operatorname{diag}(D_a,D_b)\right].
$$

$\lambda_k$ is a nonnegative eigenvalue of $-\Delta_V$; $\sigma_k>0$ predicts initial amplification of that mode. With our parameters, the candidate instability band is $6.492189<\lambda<38.507811$. A finite cell graph must actually contain an eigenvalue in that band. This calculation omits ongoing deformation, dilution, division, and nonlinear chemical changes. In particular, stabilization of uniform chemistry need not destroy a different, already-developed state.

We test that distinction directly on the six frozen time-150 graphs from seeds 7, 8, and 9. Each graph receives 42 chemical starts: the actual developed state, the exact uniform state, twenty small perturbations of the developed state, and twenty of the uniform state. Perturbations use multiplicative log-noise scales 0.01 and 0.001, respectively, with each species rescaled to preserve its starting amount. Each trajectory runs for 240 chemical time units.

A patterned trial must retain $s_a>0.1$ and return within $10^{-4}$ volume-weighted chemical log-RMS distance of the unperturbed patterned continuation. A uniform trial must return within maximum absolute log deviation $10^{-4}$ of $(1,1)$. All endpoints must have maximum chemical derivative below $10^{-6}$ and negative largest real eigenvalue of the full frozen-chemical Jacobian. The latter is evaluated at the actual, possibly nonuniform state, using the formula in the appendix. These checks test local attraction on fixed graphs, not stability of all moving variables. The 252 trajectories are repeated chemical interventions on three developmental histories, not 252 embryos. Seeds 8 and 9 use the original assay horizon, perturbation amplitudes, and acceptance criteria unchanged.

## Numerical implementation and acceptance checks

The main runs use a $72^3$ grid on $[-2.24,2.24]^3$, interface parameter 0.085, and timestep 0.0075. Observations are recorded every 0.6 model-time units. Each step advances chemical reactions and exchange on the current graph, updates polarity, advances mechanics, and applies volume dilution. Progressive cleavage and amount-balanced chemical partition follow. This is a first-order split update, even though the frozen chemical substeps use a second-order method. Division mechanics and substep restrictions are detailed in the appendix.

Numerical quality is checked separately from chemical outcomes: every-step individual volume error must remain below 5%, equivalent cell radius must remain at least four grid spacings, phase-field clipping must remain zero, and sampled boundary occupancy must remain below 0.01. A numerical failure is not counted as biological pattern loss.

For the long survival refinement, both seed-7 branches are repeated from exactly the same time-90 physical state with timestep 0.00375. Matching physical times and cell IDs, we require maximum chemical log-RMS discrepancy below 0.02, relative aggregate axis-ratio discrepancy below 1%, and agreement in the final contrast and association classifications. These practical tolerances were fixed before the repeat. Passing two timesteps bounds sensitivity for this test; it does not establish an asymptotic convergence order or convergence of zygote development.

Frozen endpoint chemistry is integrated with DOP853 at relative/absolute tolerances $10^{-8}/10^{-10}$ and $10^{-11}/10^{-13}$. Maximum absolute log disagreement must be below $10^{-5}$. Radau independently checks the most tolerance-sensitive trajectory on each graph, and the chemical Jacobian is checked by finite differences. Seeds 8 and 9 were selected before their developmental outcomes were known; a seed that fails to form a mature pattern remains reported rather than being replaced. Full spatial convergence and the accuracy of the live conductance approximation remain separate questions.

For the newer neighbor and delayed assays, DOP853 and Radau independently integrate every trajectory at relative/absolute tolerances $10^{-12}/10^{-14}$. Maximum solver log discrepancy must stay below $10^{-5}$. Delayed controls must start at stationary locally stable chemical endpoints and drift by less than $10^{-7}$ in network log RMS. Saved full trajectories, doses, time grids, recovery classifications, and endpoint spectra are checked again during assessment.

Accepted mature moving continuations use resident PyTorch arrays and matrix operations with custom CUDA mechanics and geometry/polarity kernels. The NumPy/SciPy and C++/OpenMP implementations remain numerical references. Four full CPU/GPU replays pass the prescribed waveform, recovery, chemical, and field comparisons; new starting contexts require additional prefix checks. Histories 8 and 9 also pass one coupled timestep halving from 0.00375 to 0.001875 in six mature exchange/retention and twelve selected same-state response continuations. This validates the tested mature regime, not GPU cleavage, full developmental spatial convergence, or a general continuum transport law. Legacy snapshot backends are archived separately from these active scientific paths.

# Results


## Development generates chemical contrast without identity labels, but direct feedback suppresses its formation

Both developmental branches reach sixteen cells by time 90 and pass the specified quality screens. In the late window, the no-feedback branch has mean log-activator dispersion 0.98610, compared with 0.05539 under direct feedback. Its median within-cell log-activator range is 0.01007, indicating comparatively little chemical drift over that window. The attributes remain continuous; these statistics do not establish two populations.

| Late-window measurement | Direct feedback | No feedback |
|---|---:|---:|
| Mean across-cell SD of log activator | 0.05539 | 0.98610 |
| Median within-cell range of log activator | 0.02151 | 0.01007 |
| Mean across-cell SD of log inhibitor | 0.05153 | 0.32426 |
| Mean aggregate axis ratio | 1.33328 | 1.33299 |

The nearly identical aggregate axis ratios do not support feedback-driven elongation in this comparison. Chemical and geometric symmetry breaking therefore cannot be treated as interchangeable outcomes. The result concerns one paired developmental initialization, not the frequency of organization across embryos.

## The feedback branch loses transport-supported homogeneous instability

Both sampled developmental branches support unstable modes between times 12 and 60. The direct branch loses its last unstable sampled mode between times 60 and 66. At time 90, its largest transport eigenvalue is 5.41025, below the instability band's lower edge; the no-feedback value is 6.63631, supporting one unstable mode. Maximum spatial growth rates are respectively -0.09002 and +0.00968.

![Developmental spectra reconstructed from saved phase fields. The shaded band indicates homogeneous-state instability for the fixed reaction and diffusion coefficients. These are conditional frozen-geometry diagnostics, not stability eigenvalues of the full moving model.](figures/developmental-spectrum.png){width=95%}

Total unique-edge conductance is 2.81436 under direct feedback and 3.55781 without it. Positive-edge conductance variation decreases, but weighted-degree variation increases; a general claim that feedback homogenizes the geometry is therefore unsupported. Algebraic graph controls show that reducing the no-feedback graph's overall conductance to the direct total is sufficient to stabilize uniform chemistry. Equalizing positive edge weights at fixed support and total conductance also stabilizes it. Swapping volume capacities has little effect. These operator interventions identify sufficient spectral changes, not physically realizable deformations or a unique mechanical cause.

## Polarity-dependent mechanics reproduces suppression of formation

The twelve matched component continuations have completed through time 78. Baseline, tension-only, adhesion-only, and combined tension-plus-adhesion mechanics permit sustained formation under the declared criterion. Polarity-only and full coupling do not. In contrast, all six mechanical choices retain the prepared pattern.

| Mechanical arm | Formation: late minimum $s_a$ | Persistence: late minimum $s_a$ |
|---|---:|---:|
| Baseline | 0.12322 | 1.14875 |
| Tension | 0.12593 | 1.14786 |
| Adhesion | 0.12289 | 1.16040 |
| Tension plus adhesion | 0.12558 | 1.15970 |
| Polarity | 0.02212 | 1.20381 |
| Full | 0.02243 | 1.21569 |

The values above are minima across time 63--78, not final-time values or averages. The acceptance threshold is 0.1. Polarity-only formation reaches a maximum $s_a$ of only 0.02860 over the entire run, compared with 0.02901 under full coupling. Their final homogeneous-state growth rates are -0.044791 and -0.044806, respectively; both first become negative at sampled time 65.4. Tension and adhesion alone have positive final growth rates and satisfy the formation criterion. The successful formation trajectories remain time dependent, so crossing the contrast threshold does not establish equilibrium.

These interventions show that polarity-dependent mechanical action is sufficient to reproduce suppression in this tested setting. Removing only polar mechanics from full coupling restores formation: the tension-plus-adhesion arm reaches final contrast 0.83349, with late-window minimum 0.12558, and retains prepared contrast with late minimum 1.15970. This supports necessity within the tested combination, in addition to sufficiency of the polarity-only intervention. Nor do they show that suppression is independent of chemical modulation, because activator still enters the polarity equation. All twelve continuations pass the numerical-quality screens, but this component study has not been independently refined in timestep or repeated on additional starting geometries.

![The original ten matched mechanical controls. Formation begins near uniform chemistry; persistence begins with a pattern prepared on the same frozen starting graph. Polarity-only and full coupling suppress formation, while all five displayed mechanical choices retain prepared contrast. The additional tension-plus-adhesion control is reported in the text. Spectral curves refer to instantaneous homogeneous-state stability.](figures/component-controls.png){width=95%}

## An established pattern survives switching on direct mechanical feedback

For developmental seed 7, starting from the actual no-feedback state at time 90, both moving branches retain contrast and initial-cell association through time 150. The switched branch has late minimum log-activator dispersion 1.14762 and minimum correlation 0.998965. The corresponding unchanged-control values are 1.08683 and 0.999485.

| Moving branch | Minimum $s_a$, times 135-150 | Minimum $\rho$, times 135-150 | Final homogeneous growth rate |
|---|---:|---:|---:|
| Switch feedback on | 1.14762 | 0.998965 | -0.18052 |
| Keep feedback off | 1.08683 | 0.999485 | -0.08755 |

![Paired continuations from the same developed state. Chemical contrast and initial-cell association remain high while instantaneous homogeneous-state growth rates become negative. The frozen reference evolves only chemistry on the original time-90 graph.](figures/moving-survival.png){width=100%}

Maximum per-cell volume errors are 1.139% and 1.062%, respectively; the minimum equivalent radius is 5.108 grid spacings, clipping is zero, and sampled boundary occupancy is below $1.5\times10^{-10}$. Both branches pass the quality screens. The frozen reference retains late minimum log-activator dispersion 1.01453.

At time 150, neither moving branch supports an unstable homogeneous mode. The largest transport eigenvalues, 4.64586 and 5.43490, both lie below the instability band. Nevertheless, neither developed pattern is erased over the observed horizon. Thus the developmental suppression result cannot be interpreted as universal destruction of chemical organization by direct feedback.

## Two additional developmental histories reproduce moving-pattern survival

Seed 8 develops independently from a fresh zygote with the same model parameters and reaches time 90 with sixteen cells and no active divisions. Its starting chemical contrast is 1.05438, compared with 1.01091 for seed 7. Its aggregate axis ratio is 1.37398 rather than 1.33501, so the survival challenge begins from a different geometry as well as a different chemical history.

Both seed-8 branches complete time 90--150 and pass the same criteria:

| Seed-8 branch | Minimum $s_a$, times 135--150 | Minimum initial-cell correlation |
|---|---:|---:|
| Switch feedback on | 1.17351 | 0.998701 |
| Keep feedback off | 1.12224 | 0.999202 |

The final homogeneous-state growth rates are -0.18244 and -0.09654, with no unstable homogeneous modes in either branch. Thus both completed developmental histories retain an established pattern while their endpoint graphs stabilize uniform chemistry. Maximum per-cell volume error during seed-8 continuation is below 1.165%, minimum equivalent radius exceeds 5.106 grid spacings, clipping is zero, and sampled boundary occupancy is below $8.53\times10^{-10}$.

Seed 9 also completes development and both moving continuations. Its late-window minimum contrast is 1.25722 with feedback and 1.20476 without it; corresponding minimum initial-cell correlations are 0.999587 and 0.999725. Its final homogeneous growth rates are -0.12015 and -0.04669. All quality checks pass.

All three planned histories therefore retain developed patterns while all six endpoint graphs stabilize uniform chemistry. No seed was replaced or excluded. These are replications of maintenance after no-feedback development, not direct-feedback development from the zygote or inheritance after division. Three histories remain insufficient for a reliable population success probability.

## The long seed-7 survival result passes timestep halving

Both fine-timestep branches have completed through time 150. Each passes the chemical, aggregate-shape, outcome-agreement, and numerical-quality checks against its coarse-timestep reference:

| Branch | Maximum chemical log-RMS error | Maximum relative axis-ratio error |
|---|---:|---:|
| Switch feedback on | 5.54e-05 | 0.000383% |
| Keep feedback off | 5.34e-05 | 0.000578% |

The declared limits are 0.02 for chemical error and 1% for shape error. Fine-timestep late minimum contrasts are 1.14758 with feedback and 1.08679 without it; corresponding minimum correlations are 0.998965 and 0.999485. Both classifications agree with the original runs. These results support temporal robustness of the finite-horizon maintenance conclusion from the fixed seed-7 time-90 state. They do not validate the timestep used to generate that developmental state, the component-screen trajectories, or spatial convergence.

## Uniform and patterned chemical attractors coexist on all six endpoint graphs

On every fixed endpoint graph, uniform and developed families approach distinct locally stable chemical equilibria. All twenty perturbations around each family return to the corresponding reference. All 252 trajectories pass the original numerical and attraction criteria, including 168 added trials on seeds 8 and 9.

| Seed | Feedback | Patterned equilibrium $s_a$ | Largest real eigenvalue, patterned | Largest real eigenvalue, uniform |
|---|---|---:|---:|---:|
| 7 | On | 1.18186 | -0.30409 | -0.18052 |
| 7 | Off | 1.10734 | -0.16378 | -0.08755 |
| 8 | On | 1.20362 | -0.37154 | -0.18244 |
| 8 | Off | 1.14039 | -0.21548 | -0.09654 |
| 9 | On | 1.28962 | -0.58994 | -0.12015 |
| 9 | Off | 1.22481 | -0.41868 | -0.04669 |

![Seed-7 chemical continuations on the two frozen endpoint graphs. Developed and near-uniform starts approach distinct equilibria under identical equations on each graph. The table adds the same assay for seeds 8 and 9. Perturbations probe local attraction, not independent developmental histories.](figures/endpoint-attractors.png){width=95%}

Across all three histories, the maximum tolerance-refinement discrepancy is $2.79\times10^{-6}$, independent Radau discrepancies are below $2.05\times10^{-9}$, and the maximum endpoint derivative is below $7.55\times10^{-10}$. Negative Jacobian eigenvalues together with these residual and perturbation checks support at least two locally stable chemical attractors on each graph. They do not enumerate all attractors or establish stability of the coupled moving model.

## Chemical exchange changes behavior in the same three histories

The frozen exchange-response assay revisits both endpoint graphs in each of histories 7, 8, and 9. All three histories produce donor-nearer exchanged-cell pulse responses. Across their nested graph, target, exchange, and pulse conditions, all 48 reference comparisons favor the donor. However, neither a nearer reference nor reversed activator ordering means the actual transferred concentration vector stays unchanged. Exchange can reorganize the collective chemical state, so this result concerns history-sensitive behavior within a network.

The moving assay preserves each mature geometry, polarity, and random state while preparing untouched, fresh-exchange, or pre-relaxed-exchange chemical backgrounds. It repeats the same design in histories 8 and 9 with history 7 as a historical anchor. Each history produces eight donor-nearer reference comparisons, including both pulse signs and both selected cells in the two exchange backgrounds. No planned history is replaced. These measurements represent **three reused developmental histories**, not 24 independent replicas. Prepared chemistry coevolves with geometry for time 150--210, and responses are assessed from time 210--270. The fresh and pre-relaxed states reorganize rather than meeting the declared likeness criterion for unchanged destination or transferred concentrations.

All new moving continuations pass numerical quality and starting-context backend checks. One timestep halving on histories 8 and 9 preserves the exchange/retention classifications and selected negative-pulse recovery times. Maximum normalized response discrepancy is $2.206\times10^{-4}$ and relative target response-integral discrepancy is $6.156\times10^{-4}$, below the respective 0.01 and 0.02 limits. Responses are repeated from the original time-210 physical states, separately from the refined mature exchange trajectories. Positive and pre-relaxed responses are not covered by that selected refinement. These are mature-basin tests, not fresh development from an independently perturbed zygote.

## Neighbor support changes late states and responses in every history

Changing only direct neighbors on the frozen time-210 graphs alters cell-associated chemistry and behavior while preserving the target's initial concentrations and geometry. All challenges meet the prespecified size and reference-separation criteria. The history-level results are:

| History | Late state shifts / nested challenges | Response shifts / nested comparisons | Conservative recipient state ratio |
|---|---:|---:|---:|
| 7 | 4 / 12 | 17 / 24 | 0.342 |
| 8 | 4 / 12 | 15 / 24 | 0.321 |
| 9 | 3 / 12 | 14 / 24 | 0.583 |

The final column concerns conservative neighbor averaging around the initially low-activator recipient in each fresh-exchange background. Each ratio exceeds the declared 0.1 shift threshold. Thus an amount-preserving change in surrounding chemical organization shifts that recipient's late state in **all three histories**. The target is unchanged at intervention, all other initial changes are confined to its direct neighbors, and the graph is fixed. Subsequent differences arise through transport and nonlinear network dynamics.

For completeness, the nested totals are 11 state shifts among 36 challenges and 46 response shifts among 72 comparisons. These totals describe this intervention set; they are not estimated probabilities across embryos. All 72 nearest-reference descriptions still favor the target's own untouched response, illustrating that a nearest label can remain unchanged despite a substantial waveform change. Continuous distances and recovery therefore provide information that the nearest-reference label discards.

![Neighbor-context measurements grouped by history, background, and target. The first two panels report continuous state and response changes; the lower panels show the two immediate negative-pulse nonrecoveries, both selected from history 8. Rows and pulse conditions are nested interventions, not independent developmental histories.](figures/neighbor-context.png){width=95%}

All unperturbed controls settle to locally stable chemical endpoints. Two immediate negative pulses in history 8 fail the recovery criterion: target 19 after a low-neighbor reset, and target 20 after conservative mixing in the untouched background. Their own untouched counterparts recover. An exploratory inspection of the two nonrecovering pulse endpoints finds residuals below $5.4\times10^{-13}$ and largest real chemical Jacobian eigenvalues -0.601 and -0.336. Each differs from its own control. This supports alternative stable network endpoints selected during reorganization, rather than a numerical integration failure. It does not imply that the control equilibrium, which the pulse never started from, is unstable.

## The selected history-8 failures recover when pulses are delayed

The delayed assay revisits all and only those two cases from history 8. Every pulse applied after settling returns to its own equilibrium, under both the current-percentage and original-amount dose definitions. This holds for the challenged control endpoints, untouched controls, and the alternative endpoints selected by the earlier pulses. The primary challenged-control results are:

| Case in history 8 | Delayed dose | Target recovery, negative / positive | Network recovery, negative / positive |
|---|---|---:|---:|
| Cell 19, low-neighbor reset | Current $\pm10\%$ | 4.0 / 4.0 | 4.5 / 4.6 |
| Cell 19, low-neighbor reset | Same original amount, now $\pm0.397\%$ | 4.0 / 4.0 | 4.5 / 4.5 |
| Cell 20, conservative mixing | Current $\pm10\%$ | 3.0 / 2.9 | 6.3 / 6.5 |
| Cell 20, conservative mixing | Same original amount, now $\pm11.458\%$ | 3.0 / 2.9 | 6.3 / 6.5 |

Times are chemical time since each pulse. The interval used to settle the network is frozen chemistry, not further moving development. There are eight primary challenged-control pulses, eight untouched-control pulses, and eight alternative-endpoint pulses, all within **one selected history**. The two untouched endpoint jobs share a chemical network and test different targets. None of these 24 pulses reaches a different stable endpoint. Maximum independent-solver log discrepancy is $4.05\times10^{-10}$, and maximum unperturbed control network drift is $5.43\times10^{-11}$.

![Immediate and delayed negative pulses in two selected cases from history 8. Each curve is compared with its own control; time is measured since the pulse, so starting chemical networks differ. Red curves select alternative endpoints during reorganization. Blue and green pulses applied to the settled challenged control recover under percentage and absolute-dose matching. The figure does not represent additional histories.](figures/neighbor-context-delayed.png){width=95%}

The timing contrast persists when the removed chemical amount is matched, so a changed absolute dose cannot by itself explain it. Both selected attractors resist the tested later pulses. This establishes sensitivity to the network state at perturbation alongside robustness of the resulting equilibria. It does not identify a critical delay window, the responsible concentrations, a switching probability in other histories, or persistence under moving geometry.

## Supporting release controls distinguish network support from isolation

In the one-history isolation screen, all sampled cell states approach the same positive reaction-only equilibrium after transport is removed. Shared clamped reservoirs can instead support two stable chemical equilibria under identical surroundings. A finite, evolving reservoir also supports established differences at weak and intermediate exchange and formation from small perturbations at strong exchange. These controls demonstrate possibilities within the supplied reaction-exchange system; they do not show metabolism-independent self-sufficiency or infer new cell types.

The finite-reservoir exceptions are retained: one original strong-release trajectory fails its numerical agreement criterion, symmetric stationary-looking releases can have a positive chemical eigenvalue, and one independently checked perturbed release remains evolving at the horizon. Successful formation and selected release outcomes are therefore qualified separately. All reservoir and isolation interventions descend from the same registered history; thousands of sampled chemical starting states do not enlarge the developmental replication unit.

# Discussion

## What the controls establish

The controls establish three different properties. First, on the tested component geometry, polar mechanics prevents initiation from the chosen small perturbation. Second, developed patterns persist after the coupling switch in all three histories, and locally stable uniform and patterned chemistry coexist on both resulting graphs in each history. Third, chemical history and surrounding chemical organization both influence subsequent behavior. Exchange responses are donor-nearer, while a neighbor-only change can shift late chemistry without changing the target's initial state or geometry. These observations support collective maintenance and network-supported context sensitivity. They do not establish autonomous cells or stability of the complete moving system.

The delayed follow-up sharpens the distinction between formation and robustness. During reorganization, a small pulse can redirect the collective state toward another attractor; after settling, the tested percentage-matched and amount-matched pulses return. The observation is conditional on two selected cases in one history. A response can therefore carry chemical history without preserving a transferred concentration vector, and a persistent cell-associated state can still require its network context. This is a dynamical interpretation of the model, not a demonstrated inherited biological identity.

The component controls identify polarity-dependent mechanics as sufficient to reproduce suppression and necessary within the tested full-coupling combination: tension plus adhesion without polar mechanics restores formation. This conclusion applies to one mature starting geometry, one perturbation, one coefficient set, and the tested horizon. Near unit activity, direct tension modulation is small and adhesion modulation multiplies two small responses, whereas exposure drives polarity at baseline strength. This supplies a plausible explanation for early transport changes, but chemical modulation of polarity remains active. Holding that modulation at one is still required to test independence from chemical regulation.

## Relation to previous theory

Our endpoint coexistence result is consistent with established nonlinear network theory [@nakao2010; @wolfrum2012]. It should not be presented as a new form of bistability. Likewise, the relation between network structure, spectral gaps, and pattern diversity has recent explicit treatment [@hajali2025]. Our extension is experimental in the computational sense: isolate mechanical contributions, change coupling after development, follow chemistry and geometry together, and challenge chemical states through exchange, neighbor changes, and pulse timing. The separation between initiation, maintenance, and state-sensitive reorganization is more informative than interpreting every stable concentration difference as a cell type. A continuation analysis could test whether these equilibria connect through a subcritical bifurcation, but the present experiments do not determine that bifurcation structure. A cyclic intervention would be required to demonstrate hysteresis.

Earlier morphogen-curvature and three-dimensional vertex models establish that mechanics can participate in symmetry breaking [@mercker2013; @yu2024]. Our observed suppression is not a contradiction: changing transport can remove an unstable mode as well as create one, and the outcome depends on kinetics, geometry, and initial state. We have not shown that polar mechanics always suppresses organization or that its chemical consequences dominate morphology. Indeed, the original developmental branches have almost equal aggregate axis ratios despite very different chemical contrasts.

## What wet-lab evidence supports, and what it does not

The experimental literature supports several principles reproduced by the model: activator-inhibitor self-organization, mechanical regulation of chemical communication, different requirements for initiation and maintenance, and history- and neighbor-dependent responses. Support for a principle is distinct from validation of our particular equations. No study identified in this targeted review tests the full sequence proposed here: apical-polarity-dependent changes in measured chemical conductances removing unstable graph modes, followed by persistence of an established pattern and transferable pulse-response tendencies on the same moving aggregate.

**Chemical differences without prescribed types.** In zebrafish, fluorescence recovery and pulse-labeling measurements distinguish differential effective transport from differential clearance of Nodal and Lefty [@muller2012]. The engineered mammalian circuit of Sekine and colleagues supplies a causal patterning test: genetically identical cells form signaling domains, and restricting Lefty spread prevents their formation [@sekine2018]. The authors note that their experimental patterns may be localized, initial-condition-dependent structures rather than classical Turing patterns. Their result supports chemical self-organization without supplied cell-type labels, but it does not establish embryonic lineages or our reaction laws. Digit-patterning experiments implicate a BMP-Sox9-Wnt network containing an intracellular regulatory component and external gradients [@raspopovic2014]. Theory further shows that immobile components broaden the possible patterning networks [@marcon2016]; faster inhibitor transport is a property of our selected two-variable mechanism, not a universal developmental requirement.

**A close experimental route from mechanics to transport.** Autorino and colleagues combine genetics, live imaging, and optogenetic interventions in zebrafish to link adhesion-driven tissue rigidity with Nodal-Lefty dynamics [@autorino2026]. Fluidizing tissue extends Nodal signaling; restoring rigidity rescues its spatial and temporal profile. The proposed physical route involves closure of interstitial fluid pathways and altered ligand dispersal, together with earlier Lefty expression. This is strong evidence that changing tissue architecture can negatively feed back on an activator-inhibitor system. It involves planar cell polarity and extracellular porosity, however, rather than our apical vector and symmetric contact-exchange operator. It does not measure a chemical graph's eigenvalues crossing a Turing instability boundary. The agreement therefore supports the physical route while leaving our spectral prediction open to testing.

**Polarity and communication have context-dependent effects.** Mouse embryo manipulations connect apical-domain inheritance and contractility to cell sorting and Hippo/Yap-associated lineage regulation [@maitre2016; @korotkevich2017]. Apical-domain transplantation can initiate lineage segregation, rather than suppress it. In human stem-cell micropatterns, receptor access and BMP-induced Noggin link epithelial organization to signaling boundaries [@etoc2016]. In chicken skin, connexin perturbation disrupts primary feather buds, whereas pharmacological inhibition of gap-junction communication permits additional buds [@tseng2024]. In regenerating Hydra, reducing tissue stretching prevents organizer formation, and Wnt3 overexpression rescues regeneration without normal mechanical stimulation [@ferenc2021]. These interventions establish the importance of mechanics and communication, but their different outcomes rule out a universal suppressive sign. Our observed suppression is conditional on the tested chemistry, transport, geometry, and mechanical laws.

**Initiation and maintenance have experimental precedents.** Photosensitive chemical reaction-diffusion experiments show different illumination thresholds for pattern appearance and disappearance; within the overlap, a locally imprinted pattern persists while the surrounding homogeneous medium does not spontaneously pattern [@kaminaga2006]. This is direct experimental evidence for chemical history dependence, but it is a nonliving microemulsion without cell mechanics. Intestinal organoid experiments provide the biological counterpart: early lumen inflation impairs crypt budding, whereas established buds resist inflation [@xue2025]. Contractility perturbations and washout further support path-dependent morphology. Neither result directly demonstrates uniform and patterned chemical equilibria on the same frozen embryo graph. Our endpoint coexistence is a conditional chemical result, and our model lacks the lumen-pressure mechanism of the organoid experiments.

**History-dependent response is supported; chemistry-only transfer remains a prediction.** WNT priming changes later ACTIVIN-dependent differentiation without changing the subsequent SMAD2 signaling dynamics [@yoney2018]. Sorted epiblast stem-cell populations also retain different differentiation biases after reseeding, with some biases persisting and others relaxing [@menezes2025]. These experiments show that previous state can affect later behavior under common external cues. They do not exchange only two regulator concentrations while initially preserving shape, polarity, position, and lineage, nor do they compare the resulting pulse waveform with donor and destination references. Biological reseeding transfers transcriptional and epigenetic states absent from our model. Donor-nearer response in our assay should therefore be described as a specific computational prediction consistent with biological priming, not as an experimentally validated identity transfer.

**Neighbor support is established; amount-controlled redistribution requires a matching experiment.** Human stem-cell microcolonies show neighbor-dependent sustained signaling needed for differentiation [@nemashkalo2017]. Endogenous-protein imaging and genetic interventions in human gastruloids show that Nodal propagates through a transcriptional relay and that receiving cells help sustain signaling in source cells [@liu2022]. These findings support reciprocal network dependence, but not our passive contact-exchange approximation in every biological setting. A closer experimental analogue of the model's context intervention would alter the chemical arrangement of surrounding cells while preserving the recipient's initial state and separately controlling initial total amounts. The reviewed experiments do not perform that matched redistribution-versus-dosage comparison. They also show why extracellular pools, receptors, uptake, and relay delays may be needed for molecular calibration.

**Persistent differences need not be autonomous identities.** Xenopus muscle-precursor transplantation shows that early cells require neighboring precursors to differentiate, whereas late-gastrula cells can differentiate singly in a different tissue [@katogurdon1993]. Network dependence and autonomy can therefore be stage-specific. In human stem cells, experimentally controlled BMP histories link differentiation to integrated signaling rather than an instantaneous concentration [@teague2024]. In mouse intestinal regeneration and organoids, perturbations of mechanical context, FOXA1, and Delta-Notch signaling constrain a mechanism that stabilizes transient YAP1-associated heterogeneity into spatial fate patterns [@schwayer2026]. The latter is a close biological precedent for combining intrinsic history with neighbor feedback. It includes chromatin regulation and molecular lineage readouts, however, and uniform YAP1 reactivation can disrupt established patterning. None of these results makes our two chemical activities a complete, autonomous, or inherited identity program.

## Experimental predictions that distinguish the proposed mechanism

A useful validation would compare the same mechanical intervention applied before and after a chemical pattern develops, with cell number, culture conditions, and measured mechanical state matched as closely as possible. The prediction is a range in which early intervention prevents contrast but late intervention permits recovery after small chemical disturbances. Viability, cell division, reporter turnover, and irreversible differentiation must be monitored so that apparent persistence is not simply a long-lived fluorescent label or a changed cell population. This is a proposed experiment, not an outcome established by the cited studies.

To test the transport explanation, three-dimensional segmentation should be paired with measurements of signal exchange, receptor access, and reaction response. Cell contact area alone is insufficient to infer a molecular conductance. A fitted volume-weighted operator and measured local reaction Jacobian would then permit a prospective comparison between spatial growth predictions and observed small-perturbation responses. A strong test would hold geometry approximately fixed while perturbing signaling, and separately alter geometry while checking whether reaction responses changed. If suppression occurs with unchanged transport but altered receptor sensitivity, the current proposed route would need revision.

For O1, cells or small groups should also be moved between contexts while measuring multiple molecular markers and function across subsequent divisions. Loss of state after separation would support collective maintenance; persistence under matched external signals could support a more autonomous memory. Neither outcome alone identifies the molecular circuit. Synthetic patterning cultures offer a possible controlled starting point, but a three-dimensional version would require additional engineering and transport validation rather than a direct reuse of our dimensionless parameters.

## Limits of the current explanation

Our chemical states are network attributes, not demonstrated biological cell types. Survival begins after the sixteen-cell cap, so it does not establish inheritance. The model omits a specified gene-regulatory identity program, metabolic constraints, explicit extracellular signaling, and the lumen mechanics needed for cavitation. There are three distinct developmental histories in the central evidence, reused across maintenance, exchange, and neighbor assays. Their nested interventions do not provide a population-level success estimate. Polarity-component necessity is tested on one geometry, and delayed-pulse recovery is a selected one-history follow-up. The component assay still requires timestep and geometry replication; the neighbor challenge and timing comparison still require moving-geometry confirmation.

Nutrients, biosynthetic precursors, and energy availability are assumed permissive throughout the simulated interval. Metabolism and resource limitation are not modeled. Cleavage partitions existing cell volume rather than increasing tissue biomass: daughter target volumes sum to the mother's target volume, and there is no growth between divisions. Actual volumes can deviate under mechanics. Cell-cycle timing, regulator production, polarity dynamics, and contractile-ring forcing do not depend on nutrient, oxygen, or ATP availability. Our results therefore concern organization under permissive conditions; they do not explain nutrient-dependent growth, division arrest, or metabolic self-sufficiency. The finite signaling-reservoir assays do not remove this assumption, because their exchange conserves regulator amounts while reactions still produce and consume regulators without a tracked nutrient or energy budget.

The accepted timestep/backend comparisons cover the tested mature continuations. They do not validate preceding cleavage or replace the failed historical full-development refinement. Curved-interface bias, diffuse-gap exchange, and nonorthogonal flux errors remain in the live conductance approximation. Separate corrected-flux prototypes have not been coupled to these simulations. Conservation of amount is essential but does not make the operator a validated continuum discretization on arbitrary geometry. Spatial refinement, transport-model validation, and biological calibration remain necessary before interpreting the results as a quantitative account of living development.

# Reproducibility and evidence status

The manuscript reports completed developmental/component controls, maintenance and chemical coexistence in three histories, frozen and moving exchange responses in those same histories, neighbor challenges grouped by history, and a selected delayed-pulse follow-up in history 8. Numerical validation and supporting release controls retain their separate scopes and failures. A human-readable and machine-readable results ledger records 22 curated studies, explicit shared history IDs, validation decisions, limitations, evidence hashes, and an inventory of historical output directories. Completion alone is not numerical acceptance. The evidence manifest captures the reports supporting this revision and checks the reported table values.

The literature assessment was extended on 7 October 2026 with primary experimental studies of mechanically regulated transport, chemical and morphological history dependence, cellular priming, network support, and signaling dynamics. It is a targeted narrative review, not an exhaustive systematic search or a claim that no closer experiment exists. The accompanying literature evidence map records sources, applicability, contrary results, and access limitations. Claims about biological identities remain proposals for future validation, not results inferred from visual resemblance.

This literature revision does not replace the numerical results or their provenance. The tables and figures retain the original no-carry moving trajectories or assays inheriting their graphs and states; their original numerical decisions remain recorded in the evidence manifest and claim-level precision audit. Later carry-qualified studies in the repository are separate evidence and are not incorporated as new results in this draft. Experimental support for a biological principle does not resolve a numerical qualification or calibrate the model to a living system.

The manuscript README maps results to their implementations and output files. Retired prototype backends are archived as exact source snapshots, with pinned originals and shared dependencies retained for evidence verification. Building this revision launches no new simulation and changes no accepted model source or outcome report. Author information, funding, and acknowledgements have not yet been supplied. No public data deposit or archival DOI is claimed. An exact source snapshot, configurations, checkpoints, and analysis outputs should accompany a durable public deposit before submission.

# Appendix: implementation details {.unnumbered}

## Supporting common-environment and reservoir controls {.unnumbered}

The common-environment screen releases 2256 chemical starting vectors obtained from the original one-history network state, its pair exchanges, and seeded resets. Repeated cell states remain included for traceability. Isolation uses only $R_a,R_b$. With an infinite fixed reservoir $(A,B)$, cells instead obey $\dot a_i=R_a+k_a(A-a_i)$ and $\dot b_i=R_b+k_b(B-b_i)$. Exchange strengths are $s=0.25,1,4$, with $k_a=sD_ar$ and $k_b=sD_br$, where $r$ is the original graph's median exit rate, approximately 3.27684. Reservoir compositions are either $(1,1)$ or the original volume-weighted mean, approximately $(0.85905,0.93958)$. Every cell sees identical coefficients and reservoir concentrations. Positive steady-state roots and their two-species Jacobian spectra provide independent equilibrium diagnostics; no type labels are imposed.

In the finite-reservoir extension, direct cell-contact exchange is removed, fixed cell volumes are retained, and the reservoir responds through

$$
\begin{aligned}
\dot a_i&=R_a+k_a(A-a_i),&\qquad
\dot b_i&=R_b+k_b(B-b_i),\\
V_R\dot A&=\sum_iV_ik_a(a_i-A),&
V_R\dot B&=\sum_iV_ik_b(b_i-B).
\end{aligned}
$$

Here $V_R$ is reservoir volume. Opposite exchange terms conserve each species' amount across cells and reservoir; intracellular reactions still create or consume amount. There is no external replenishment or reservoir reaction, but reaction sources do not track nutrients or energy. This globally mixed coupling is not an extracellular spatial diffusion model. Reservoir-to-total-cell-volume ratios are 0.1, 1, 10, and 100. Release conditions start from the previous clamped equilibria at their former reservoir compositions. Formation conditions start near unit chemistry with positive amount-preserving cell perturbations of log-noise scale 0.001; five perturbation seeds per volume/strength pair are nested within the same source history.

The original finite study has 24 releases, 60 formation perturbations, twelve exactly uniform controls, six clamped controls, and one isolated control. Horizons account for the slow reservoir response and range from 1200 to approximately 122069 chemical units. Numerical agreement, endpoint derivatives, and the full coupled chemical Jacobian are assessed separately; a stationary-looking symmetric trajectory is not accepted if its largest real eigenvalue is positive. Eight additional strong-release checks add amount-preserving cell perturbations of log-noise scale $10^{-4}$ and retain one unsettled outcome. Independent DOP853/Radau checks cover those eight releases and twenty original strong-formation trajectories, with maximum log disagreement below $10^{-5}$. The original numerical failure is retained rather than relabeled by the follow-up.

The earlier common-environment screen likewise retains one original coarse/tight tolerance failure. All-state further refinement and one selected worst-state Radau check per condition support the reported tight solutions without relaxing its threshold. These distinct levels of validation are recorded in the ledger. Neither the many sampled cell vectors nor reservoir or perturbation conditions count as additional developmental histories.

## Initial cell, interface geometry, and conductance {.unnumbered}

The initial radius is $R_0=0.8$, and the founder field is

$$
\phi_0(\mathbf x,0)=\frac12\left[1-\tanh\left(\frac{|\mathbf x|-R_0}{\sqrt2\epsilon}\right)\right].
$$

Target volume is the initial discretely measured occupancy. Cell centroids use the same occupancy weights. With $\omega_i=\phi_i(1-\phi_i)$, the contact approximation is

$$
W_{ij}=\int_\Omega\omega_i\omega_j\,\mathrm d\mathbf x,\qquad
\widehat{\mathcal A}_{ij}=6\sqrt2 W_{ij}/\epsilon,\qquad
g_{ij}=\widehat{\mathcal A}_{ij}/|\mathbf c_i-\mathbf c_j|.
$$

Diagonal conductances are zero; symmetric filtering removes overlaps below 0.02 times the largest overlap. The area factor comes from integrating complementary planar equilibrium interfaces. A positive contact between coincident centers is rejected. On frozen geometry, the symmetric matrix $M^{-1/2}KM^{-1/2}$ has the same eigenvalues as $-\Delta_V$. Constant concentrations are preserved by exchange, and $\mathbf1^TM\Delta_V=0$ expresses amount conservation. Neither statement proves the geometric flux approximation accurate on general cell contacts.

## Mechanical cost and force {.unnumbered}


For coefficients held fixed during a mechanical update, the effective cost is

$$
\begin{aligned}
E={}&\sum_i\int_\Omega\gamma_i(\mathbf x)
\left[\frac{\epsilon^2}{2}|\nabla\phi_i|^2+Q(\phi_i)\right]\mathrm d\mathbf x\\
&+\sum_i\frac{K_V}{2V_i^\star}(V_i-V_i^\star)^2
+\frac{R}{2}\sum_{i<j}\int_\Omega\phi_i^2\phi_j^2\,\mathrm d\mathbf x\\
&-\sum_{i<j}A_{ij}\int_\Omega Q(\phi_i)Q(\phi_j)\,\mathrm d\mathbf x.
\end{aligned}
$$

The gradient and double-well terms regularize the interface; $K_V$ penalizes departure from target volume; $R$ penalizes overlap; and $A_{ij}$ rewards contact between interfaces. This is a mechanical cost, not ATP expenditure or an assertion that development minimizes a fixed energy. Baseline values are $\gamma_0=1$, $K_V=12$, $R=3$, and $A_0=4$.

With unit mobility, the implemented force is

$$
\begin{aligned}
\partial_t\phi_i={}&\epsilon^2\nabla\cdot(\gamma_i\nabla\phi_i)-\gamma_iQ'(\phi_i)
+K_V\frac{V_i^\star-V_i}{V_i^\star}H'(\phi_i)\\
&-R\phi_i\sum_{j\ne i}\phi_j^2
+Q'(\phi_i)\sum_{j\ne i}A_{ij}Q(\phi_j)+F_i^{\mathrm{furrow}}.
\end{aligned}
$$

Centroids, polarity, chemical coefficients, and ring coordinates are fixed when evaluating this force; their implicit dependence on the phase fields is not differentiated. Zero normal mechanical flux is imposed at the box boundary. Active chemistry, changing coefficients, and cytokinesis preclude interpreting the full evolution as passive relaxation of one time-independent functional.

## Polarity geometry and progressive cytokinesis {.unnumbered}

The exposure cue $\mathbf q_i$ is the outward normal averaged over unoccupied cortex, with total shell weight as denominator. The normal is $-\nabla\phi_i/\max(|\nabla\phi_i|,10^{-10})$ and its weight is $\omega_i[1-\operatorname{clip}(2\sum_{j\ne i}H(\phi_j),0,1)]$. Cues below magnitude $10^{-6}$ are set to zero. After a forward-Euler update, polarity is projected into the unit ball. Directional tension is

$$
\gamma_i(\mathbf x)=\gamma_i^0\left[1-\chi\mathbf p_i\cdot
\frac{\mathbf x-\mathbf c_i}{\sqrt{|\mathbf x-\mathbf c_i|^2+\epsilon^2}}\right].
$$

The spatial coefficient enters the divergence-form mechanical force, including its gradient. Baseline mechanics sets $c_\gamma=c_A=\chi=0$ while retaining polarity dynamics.

Each newborn cell receives a waiting time $2U$, with $U$ uniform on $[0.88,1.12]$. At division onset, the cleavage-plane normal is chosen from the largest-eigenvalue subspace of the occupancy covariance; eigenvalues within relative tolerance 0.03 are treated as degenerate. A random vector within that subspace provides an unbiased fallback for nearly spherical cells. This is a shape-based Hertwig-style rule, not a stress-tensor calculation or an explicit model of spindle machinery [@middelkoop2024].

Cytokinesis progresses through an annular exclusion force over nominal duration $T_{\mathrm{cyt}}=0.9$. With progress $u=\operatorname{clip}[(t-t_0)/T_{\mathrm{cyt}},0,1]$, ramp $s=u^2(3-2u)$, axial distance $z$, and radial distance $r_\perp$,

$$
\begin{aligned}
R_{\mathrm{ring}}&=R_{\mathrm{ring},0}(1-s),\\
F_i^{\mathrm{furrow}}&=-8s\exp[-z^2/(2\epsilon^2)]
\frac{1+\tanh[(r_\perp-R_{\mathrm{ring}})/\epsilon]}{2}H'(\phi_i).
\end{aligned}
$$

The plane initially bisects occupancy; its offset follows the moving centroid. Volume is projected back to the onset value during constriction. Relabeling is postponed until the neck phase value is at most 0.15, relative daughter overlap at most 0.002, and both daughter occupancy fractions exceed 0.1. Daughter regulator amounts are partitioned with balanced noise of scale 0.001, then converted using measured daughter volumes. This avoids creation of chemical amount by a volume mismatch. Division is an imposed biological process, not an outcome derived from activator-inhibitor kinetics.

## Error measures and endpoint Jacobian {.unnumbered}

For fixed initial-volume weights $w_i=V_i(90)/\sum_jV_j(90)$ and log activities $x_i=\log a_i(t)$, $y_i=\log a_i(90)$, the cell-association measure is

$$
\rho(t)=\frac{\sum_iw_i(x_i-\bar x_w)(y_i-\bar y_w)}
{\sqrt{\sum_iw_i(x_i-\bar x_w)^2\sum_iw_i(y_i-\bar y_w)^2}}.
$$

The moving-timestep chemical discrepancy between two states is

$$
E_{\mathrm{chem}}=\sqrt{\frac12\sum_iw_i\left[
\left(\log\frac{a_i^{(1)}}{a_i^{(2)}}\right)^2+
\left(\log\frac{b_i^{(1)}}{b_i^{(2)}}\right)^2\right]}.
$$

The frozen-endpoint return metric uses the same expression with normalized fixed endpoint volumes. Aggregate axis ratio is the square root of the largest-to-smallest eigenvalue ratio of the covariance of the capped occupancy union, $\min(\sum_iH(\phi_i),1)$. Timestep shape error is $|R_{\mathrm{coarse}}/R_{\mathrm{fine}}-1|$ at matched times.

For frozen geometry, the complete chemical Jacobian at a nonuniform state is

$$
J_{\mathrm{chem}}=\begin{pmatrix}
\operatorname{diag}(2a_i/b_i-1)+D_a\Delta_V & -\operatorname{diag}(a_i^2/b_i^2)\\
2\beta\operatorname{diag}(a_i) & -\beta I+D_b\Delta_V
\end{pmatrix}.
$$

It excludes derivatives of mechanics and polarity. For the homogeneous state, the modal determinant is $\beta+(\beta D_a-D_b)\lambda+D_aD_b\lambda^2$, whose roots give the stated instability band when local kinetics are stable.

The chemical substep uses SSP-RK2 with loss-rate bound $\delta t\max(1+D_ar_{\max},\beta+D_br_{\max})\le0.2$, where $r_{\max}=\max_i\sum_jg_{ij}/V_i$. After mechanics, concentrations are multiplied by $V_i^{\mathrm{old}}/V_i^{\mathrm{new}}$. Checkpoints are saved every six model-time units and include state, random streams, observation history, and audit information. The unit tests check exact restart, conservative exchange and division, and the intended coefficient changes in controlled interventions.
