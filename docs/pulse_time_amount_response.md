# Challenge-time local-state and surrounding-amount control

## Question and design

The [accepted challenge-time matching experiment](pulse_time_recipient_response_assessment.md) shows different immediate pulse responses despite equal starting recipient activator/inhibitor and a common frozen physical context. Current surrounding amounts remain unequal. Does a difference in **spatial chemical allocation** still alter response when both surrounding species totals are also equal at the pulse?

Use its accepted t=570 preparations on the twelve actual co-evolved endpoint graphs. Six selected recipients remain nested in the same three developmental histories (7–9); this adds zero histories. Physical S/R labels mean endpoint geometry co-evolved with S/R chemical allocation, with both starting from original sham physics at t=510. Freeze each measured conservative operator and volume vector for t=570–630.

Retain all native/no-reset and locally matched/unmatched-budget controls. For both measured local references M=S/R, additionally cross two measured surrounding-total references B=S/R with both chemical allocations C=S/R on both shared physical contexts P=S/R. Native/no-local-reset chemistry is only an unmatched-budget comparator; no extra budget-only native arms are added. Every trajectory is calculated anew; no prior scientific trajectory is reused.

For each species s, recipient i and common capacity vector V, define its surrounding amount

$$
Q_s=\sum_{j\ne i}V_j u_{sj}.
$$

The common target amount is the measured amount of allocation B on that same physical context, excluding the recipient. Scale every surrounding concentration of species s by the same positive factor:

$$
u'_{sj}=\frac{Q_s^{(B)}}{Q_s^{(C)}}u_{sj}\quad(j\ne i),
\qquad u'_{si}=u_{si}.
$$

Here u is activator or inhibitor concentration. The measured recipient pair is already equal across C. Uniform scaling within each species preserves its relative surrounding allocation fractions V_j u_sj / Q_s. It does not preserve absolute concentrations, local reaction rates or the activator/inhibitor ratio at each surrounding cell. These are explicit external preparation interventions, not new biological equations or conservative transport steps. Both measured budget references are tested to avoid relying on one arbitrary dosage choice.

In matched arms, recipient chemistry, surrounding species amounts and whole-network species amounts agree across allocations at the pulse time. All are free to evolve afterward. Nothing is continuously clamped. Record actual per-cell changes, scaling factors, species budgets and fraction-closure errors. Inherited local-reset and physical-reassignment accounting remains pinned in the accepted parent.

## Pulse, solver and qualification

At elapsed zero apply a −10% recipient activator pulse; retain an identically prepared unpulsed control. Measure against that own freely evolving control so chemical preparation relaxation is not counted as pulse response. Equal target concentrations and common volume make the fractional pulse an equal-removed-amount pulse across each matched allocation pair. Verify amounts before and after the pulse. This is an **immediate matching/reset assay**, not a later settled challenge.

Use unchanged Gierer–Meinhardt beta=2, D_a=0.02, D_b=0.55 and the existing conservative volume-weighted operators. Use qualified compensated float64 SSP-RK2 with PyTorch GPU arrays and matrix operations on the GTX 1080 Ti, selected by UUID. Initialize numerical corrections to zero at the explicit rounded-state preparation, then retain them through every substep and checkpoint. Corrections are numerical bookkeeping, not biological memory. Shared CUDA mechanics/geometry/polarity kernels remain unchanged and inactive in this frozen control. No moving dilution is present.

There are 168 labeled contexts, 336 logical starts, 42 allocation/physics factorials, 1,008 GPU paths, 1,008 independent NumPy paths and 672 independent DOP853/Radau trajectories. Report the exact distinct physical initial-value count separately; labels and duplicate starts are not replication. Identical physical matrices/capacities/IDs form twelve independent blocks with 28 chemical columns each, no padding or communication between columns. Use eight single-thread CPU workers and BLAS threads one. Observe every 0.15 units, checkpoint every three, and use dt=0.00375, 0.001875 and 0.0009375.

Before production, verify accepted parent protocol/review/source/input hashes and every actual consumed preparation. Independently reconstruct all local target values, surrounding rescaling, per-species allocation preservation, global/surrounding amount equality, conserved graph structure, IDs, positive chemistry and exact pulses. Require relative amount/fraction errors <=2e-13. Preserve the parent's explicit non-input manuscript-export discrepancy.

Every logical start receives a 0.6-unit GPU/NumPy comparison at all three timesteps and interrupted/restored versus uninterrupted bitwise comparison of concentrations, corrections and histories, with a checkpoint at 0.3. Independently reconstruct the actual saved prefix arrays before full production. Require raw agreement <=1e-8 and normalized pulse response <=0.01. Correct descriptive inherited case-count metadata for this new design without changing shared stepping. Missing/invalid corrections fail continuation.

Generate fresh dual references with rtol=1e-13, atol=1e-15 and maximum step 0.15. Qualify raw and normalized reference/contrast agreement before production. Retain all acceptance rules: dual-reference raw/normalized discrepancy <=1e-8; GPU/NumPy raw discrepancy <=1e-8; production raw discrepancy <=0.01; normalized pulse response, signed contrast and E errors <=0.01; relative response-area discrepancy <=0.02; recovery-time discrepancy <=0.3 with matching censoring. If diagnostic error exceeds 100 times reference noise, both observed orders must lie in 1.7–2.3. Otherwise report unresolved order, not order two. Preserve all failures and do not relax any limit.

## Measurements and decision rules

The primary comparison is C=R versus C=S at fixed recipient reference M, budget reference B and endpoint physical context P. E is the 60-unit RMS difference between signed, two-species recipient log pulse/own-control responses, normalized by each actual pulse's log magnitude. Use the inherited strict E>0.01 screen. Require classifications to agree across timesteps and both independent reference methods; report each history with its two nested recipients. Below-screen is not equivalence, and recipient counts are not population probabilities.

Compare matched-budget versus unmatched-budget signed chemical contrast waveforms before taking their RMS. These budget-change magnitudes, physical contrasts and interactions are descriptive, not additive causal percentages or additional primary screens. Report both budget references, both local references and both contexts even when they disagree. Current dosage is controlled only at initiation of the pulse; subsequent reaction-driven amount divergence is an outcome, not a failed preparation.

Record pulse recovery, own-control relaxation, re-emerging recipient differences, chemical contrast and freely evolving global/surrounding amount curves. Check all numerical and primary decisions before interpretation. If differences survive both measured budget references, this supports allocation-sensitive behavior at matched initial local state and dosage. Loss or reference-dependent survival would constrain that claim; it would not prove that the original effect was solely dosage-driven because rescaling also changes surrounding kinetics.

A coordinator pass remains pending separate independent saved-array review. This tests mature immediate context sensitivity, not moving post-reset maintenance, settled response, nearest-neighbor specificity, autonomous/inherited identity, new zygote development or spatial convergence. Preserve accepted parents, prospective protocols, shared kernels, results ledger and manuscript. Output is `output/studies/pulse-time-amount-response`; no 3D fields are duplicated.
