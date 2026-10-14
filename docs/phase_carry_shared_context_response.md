# Chemical responses in reciprocal shared physical contexts

## Question and intervention

The [accepted frozen challenge-time test](phase_carry_frozen_response_assessment.md) retains the selective response differences without continuing physical evolution. Those continuations still inherit different chemistry, transport and volumes. **Do different inherited chemical states respond differently when placed in the same physical context, and does that answer depend on which physical context is chosen?**

Use the actual t=510 sham and recipient-specific redistributed states in histories 7, 8 and 9. The selected recipients remain 22/23, 19/20 and 29/30, with 16 cells per context. For each recipient, cross two chemical starts with two physical contexts:

| Chemical concentrations | Fixed physical context |
| --- | --- |
| Sham | Sham transport and volumes |
| Redistributed | Sham transport and volumes |
| Sham | Redistributed transport and volumes |
| Redistributed | Redistributed transport and volumes |

Transfer the **entire two-species concentration vector**, preserving concentrations and cell IDs exactly. The conservative operator and its measured volume vector travel together. Never attach an operator to a different capacity vector. Challenge immediately at reassignment; do not equilibrate the crossed starts first. Their own unperturbed continuations may relax after reassignment, which must be considered when interpreting response differences.

The chemical factor means **inherited concentration levels and their distribution across cells**, not spatial arrangement alone. Total amounts and recipient concentrations are not matched between chemical starts. Original redistribution matched amounts at t=450, not at this t=510 challenge. Reassignment at fixed concentrations changes per-cell amounts when volumes change. Save original and reassigned per-cell and total activator/inhibitor amounts, their changes, and volume ratios. This is a conditional chemical-state comparison, not an amount-controlled allocation experiment, an independent geometry-versus-volume decomposition, or a causal mediation fraction.

## Fixed equations and pulse preparation

Let `a_i` and `b_i` denote activator and inhibitor concentrations. For each physical context, `V_i` is its measured cell volume, `M=diag(V_i)`, and `Delta=-M^{-1}K` is its exact saved conservative transport operator. Continue unchanged Gierer–Meinhardt reactions and exchange from t=510 to t=570:

$$
\frac{d a_i}{dt}=\frac{a_i^2}{b_i}-a_i+D_a\sum_j\Delta_{ij}a_j,
$$

$$
\frac{d b_i}{dt}=\beta(a_i^2-b_i)+D_b\sum_j\Delta_{ij}b_j.
$$

Self-enhancement `a_i^2/b_i` is opposed by inhibitor and activator loss. Activator produces inhibitor, which also decays. Exchange follows measured contacts. Retain `beta=2`, `D_a=0.02`, `D_b=0.55`, and inherited directional contrast `chi=0.35`. Mechanics, polarity, volumes and transport are not advanced. Geometric dilution is zero because volumes are fixed. Exchange conserves volume-weighted species amount; reactions do not. Check positivity, finiteness, zero row sums, nonnegative off-diagonal exchange, volume-weighted conservation and reciprocal conductances without modifying the saved operator.

Each chemical–physical background has its own unperturbed control. For each target apply activator-only depletion, leaving inhibitor and other cells unchanged:

1. Fractional pulse: multiply target activator by `0.9`.
2. Common removed amount: within each shared physical context remove `deltaN=0.1*V_i*min(a_i^S,a_i^R)` from both chemical starts.

To make the dose factor exactly identical across physical contexts, compute it canonically from concentrations: `f_C=1-0.1*min(a_i^S,a_i^R)/a_i^C`. Compute requested amount separately using the chosen `V_i`; verify measured removal against it within relative error `2e-13`. Record actual factors and absolute initial log amplitudes. This is the same within-context amount design with an explicit rounding convention. Amounts rescale with physical-context volume; factors do not. These common pulses generally differ from the previous native comparison, which paired different volumes. Do not treat those older common-dose paths as identical comparators.

Run every case anew. Old native controls and fractional pulses are numerical comparators only where source chemistry, operator, volumes, dose and timestep match exactly; they are not reused scientific trajectories. There are 21 backgrounds, 69 logical starts and an expected 57 distinct initial-value problems, to be counted from prepared bytes. The sham–sham control is shared across two targets within a history. Two primary timesteps yield 138 production paths; a third yields 69 diagnostic paths, for 207 total. Independently verify all 207 with NumPy/SciPy and solve all 69 logical starts with each of DOP853 and Radau, giving 138 reference trajectories. None of these paths creates a new developmental history.

## Response contrasts and interpretation

Observe every 0.15 units for the full 60-unit horizon. For each pulse define its signed two-species response `R_C,P(t)` as the pulse/own-control log concentration ratio divided by the pulse's actual absolute initial log amplitude. Subscripts `C` and `P` denote chemical and physical starts. The primary contrast on each shared physical context is:

$$
E_P=\left[\frac{1}{60}\int_0^{60}\frac12\sum_{s\in\{a,b\}}
\left(R_{R,P,s}(t)-R_{S,P,s}(t)\right)^2dt\right]^{1/2}.
$$

Use trapezoidal integration on the observation grid. The primary operational screen remains strictly **`E_P>0.01`**, separately for both physical contexts and both doses. Report all 24 chemical contrasts and history counts out of two recipients for each context/dose; report whether each recipient clears both doses in both contexts. A below-screen result is not equivalence. Do not select a favorable physical context or collapse the two outcomes into one.

Also report continuous physical contrasts at fixed chemistry and the signed interaction:

$$
I(t)=R_{R,R}(t)-R_{S,R}(t)-R_{R,S}(t)+R_{S,S}(t).
$$

Compute the interaction waveform before its two-species time RMS. Do not subtract unsigned effect norms or interpret an interaction RMS as a mediated fraction. Physical and interaction contrasts are descriptive continuous measures, with numerical qualification, not additional scientific detection screens.

Report full waveforms, target activator response area, target and volume-weighted network recovery, and chemical contrast. Retain the previous recovery definition: all remaining sampled displacements within 10% of the initial displacement, with at least 24 later observed units; otherwise right-censor at t=570. Recovery is toward the corresponding continuing control. Describe unperturbed control drift from its reassigned start, late (elapsed 36–60) target and volume-weighted network two-species log differences, and final activator/inhibitor concentration ratios on each shared physical context. Do not require recovery, maintained contrast, or a late state difference as a scientific positive. If crossed controls converge while their immediate responses differ, distinguish transient inherited-state sensitivity from maintained state differences.

Chemical differences on a common physical context show that different measured chemistry can produce different responses conditional on that fixed environment. They do not isolate recipient-autonomous memory from the surrounding chemical network. **Recipient-state matching is a separate later test**, not included here. Dependence on physical context shows conditional sensitivity, not that one component alone created biological identity.

## Numerical gates and immutable evidence

Use unchanged float64 PyTorch SSP-RK2 (`gpu_backend.gm_step`) with GPU-resident arrays and matrix operations. Batch independent contexts with an exact block-diagonal matrix and separate trajectory columns; padding duplicates positive controls and is not a scientific path. Require the global positivity substep count to equal every context's count at every timestep. Check both internal Euler stages for positivity and finiteness; do not clip.

Use `dt=0.00375`, `0.001875`, and diagnostic `0.0009375`, with identical starts and observations. Independently solve every logical start using DOP853 and Radau at `rtol=1e-13`, `atol=1e-15`, maximum step `0.15`, retaining both. Raw and normalized dual-reference discrepancies must be at most `1e-8`. Production/reference raw log error must be at most `0.01`; CPU/PyTorch raw discrepancy must be at most `1e-8`. Where finest production error exceeds 100 times reference disagreement, observed convergence order must lie in `1.7–2.3`; report unresolved order otherwise.

Preserve the inherited pulse gates: maximum full-network normalized waveform discrepancy `0.01`, relative target activator-area discrepancy `0.02`, matching recovery/censoring and time difference at most `0.3`. For each chemical, physical and interaction contrast, directly compare signed target contrast waveforms as well as RMS norms: maximum waveform discrepancy at most `0.01` for production comparisons and `1e-8` between independent references; RMS discrepancy at most `0.01` (or `1e-8` for references). Require primary chemical-effect classifications to agree across all timesteps and references. No separate scientific classification is assigned to physical/interaction norms. Retain failures without changing thresholds or excluding targets.

Pin this protocol, runner/tests, all scientific source dependencies, accepted parent assessment and independent review/provenance, and every consumed t=510 source array. Preserve exact prepared chemistry/physical contexts and dose/amount accounting; reference, CPU and production trajectories; checkpoint/resume state and hashes; numerical gates, runtime, timing, status and summary. Audit actual sampled arrays independently of project analysis helpers before accepting a scientific conclusion. Runtime completion alone is not acceptance.

## Execution and scope

Run `python -m embryo.phase_carry_shared_context_response prepare`, then `run`, with `--output output/studies/phase-carry-shared-context-response`. Keep this small study in the writable repository workspace; read parent outputs in place and copy no 3D phase fields. Use physical GPU 1 as logical `cuda:0` with `CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1`, four CPU workers and one BLAS/PyTorch thread per worker. Validate all actual starts on a short CUDA prefix against unchanged CPU stepping before full production.

The independent units are three existing mature histories with recipients and interventions nested within them. These chemical continuations add no new zygotes or divisions, no continuing reservoirs/clamps, no inherited biological identities, and no fresh carry-corrected developmental replication. Retain qualifications concerning historical no-carry formation and prior spatial/conductance validation. This step does not update the results ledger, manuscript numerical Results, shared kernels or baseline simulator, and does not authorize a subsequent experiment.
