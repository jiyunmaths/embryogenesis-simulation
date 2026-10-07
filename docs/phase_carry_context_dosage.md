# Separate surrounding spatial chemistry from bulk chemical dosage

## Question and scope

The [completed reset](phase_carry_network_context_assessment.md) changed the recipient's maintained chemical state. The [delayed-response assay](phase_carry_delayed_response_assessment.md) found selective later behavior changes, even with matched pulse doses. Neither controls the amounts added or removed by the **original surrounding-cell reset**. Does changing the distribution of surrounding chemistry still change the recipient when those initial amounts are held fixed? Can changing bulk dosage alone reproduce the reset effect?

This first stage tests **unperturbed chemical states**, not a new pulse-response assay. Start from the original common fine **t=450 exchanged state**, before any surrounding reset. Cell IDs, phase fields, inherited float64 residual carry, rounding counters, polarity, lineage and random streams are identical initially across the four arms. Recipient activator and inhibitor concentrations remain exactly unchanged at preparation. The two chosen recipients in each history are separate interventions. Replication is **three existing mature histories**, not six independent recipients or 42 independent embryos. No fresh carry-corrected zygote history is generated.

## The four controls

For each recipient, let S be all other cells. For species s (activator or inhibitor), x_sj is the original concentration and r_sj is the same-age untouched reference concentration in cell j. V_j is its measured volume in the original exchanged physical state. Define surrounding amounts Q_s = sum over j in S of V_j x_sj and R_s = sum over j in S of V_j r_sj. Reference concentrations are evaluated using these common original volumes; the reference embryo's volumes are not substituted.

| Arm | Relative concentrations among surrounding cells | Initial surrounding amount per species | Concentration assigned to surrounding cell j |
| --- | --- | --- | --- |
| Sham | Original | Q_s | x_sj |
| Raw reset | Reference | R_s | r_sj |
| Amount-preserving redistribution | Reference | Q_s | r_sj times Q_s / R_s |
| Bulk dosage | Original | R_s | x_sj times R_s / Q_s |

Each species has its **own** positive scaling factor. All within-species concentration ratios among surrounding cells retain the assigned distribution. The recipient is excluded from rescaling; preserving surrounding amounts therefore also preserves total embryo amounts at preparation. Sham and raw reset use exact copies for compatibility with existing evidence. Record requested and actual surrounding amounts, both species' scaling factors, added/removed amounts, volume-weighted surrounding log shifts and exact recipient preservation. Relative initial amount error must be at most 2e-14.

Distribution here means a spatial allocation over existing cell IDs, not a change in the contact graph or a molecular identity label. The bulk arm preserves within-species relative concentrations but can change the activator/inhibitor ratio because its two species factors can differ. This is a one-time intervention, with no reservoir or ongoing amount clamp. Reactions, transport, dilution and mechanics can subsequently change amounts and geometry.

## Moving runs and exact reuse

Continue chemistry, mechanics, polarity, conservative volume-weighted transport and geometric dilution together for **60 model units, t=450–510**. Use dt=0.00375 and 0.001875 from the same fine physical state. Model parameters stay at beta=2, D_a=0.02, D_b=0.55, c_gamma=0.25, c_A=0.35, chi=0.35, 72 cubed grid, extent 2.24 and interface width 0.085. Use qualified PyTorch GPU arrays/matrices, custom CUDA mechanics/geometry/polarity, float32 phase fields with float64 residual carry and float64 chemistry. No earlier scientific source or kernel is edited.

There are **42 distinct paths: 24 new and 18 exact reused**. Per history and timestep, one sham is shared by both recipients; each recipient has raw-reset, redistribution and bulk arms. Reuse the completed sham and raw-reset paths only after regenerating their preparation temporarily and requiring identical array dtypes/values, all metadata, carry, rounding and random streams. Verify parent protocol, accepted independent review, checkpoint/history/prefix evidence and original pilot/full field snapshots. Reuse is not inferred from parameter agreement. Reused field snapshots are read in place. Estimated additional storage is 8.4 GiB, with a 2 GiB reserve; stop before a new continuation if the reserve falls below 2 GiB, retaining committed checkpoints.

## Prespecified comparisons

For each recipient and species, use log concentration differences between arms. Primary magnitude is the time RMS over the final **24 units**, with both species equally weighted. The effect screen is **0.01** in this log-state metric. Preserve below-screen outcomes; they do not demonstrate equivalence or define a biological cell type. Report signed final log differences and final species ratios, plus descriptive full-window RMS.

- Spatial effect at original dosage: redistribution minus sham.
- Spatial effect at reference dosage: raw reset minus bulk.
- Dosage effect at original distribution: bulk minus sham.
- Dosage effect at reference distribution: raw reset minus redistribution.
- Total historical reset effect: raw reset minus sham.
- Interaction: raw-reset log state minus bulk log state minus redistribution log state plus sham log state. Form this **signed** difference first, then take its time RMS. A nonzero interaction indicates nonadditivity in log-state effects at this tested point; it does not identify a unique biological mediator.

Report both selected recipients and aggregate outcomes by history. A spatial effect at matched initial amounts supports sensitivity to initial spatial chemical context. A dosage-only effect narrows the prior reset explanation. Effects in both factors or an interaction imply a coupled explanation. Later geometry remains part of the response; this design does not separate direct chemical action from mechanical mediation or test neighbor-only topology, network necessity, autonomous/inherited identity, pulse behavior or fresh developmental differentiation.

## Numerical acceptance and execution

Run CPU implementation/regression tests, four small GPU tests (both new arms at both timesteps), and four full-size native/carry prefixes for recipient 19 in history 8 (both new arms/timesteps). Every remaining new preparation requires its own prefix gate. Each prefix checks eight independent chemical steps (maximum log error 1e-11), amount conversion (2e-14), four exact restart steps and the unchanged 0.6-unit native/carry tolerances.

Within each history all 14 six-unit paths must pass its paired gate before new long continuations. Sample every 0.15 units; checkpoint every three units. Require raw log chemistry, polarity, relative axis and transport errors at most 0.01, relative volume at most 0.005, final phase fields at most 0.02, agreeing contrast classifications and loss-time differences at most 0.3. All state-effect classifications must agree at both timesteps; full-window and late-window RMS discrepancies must be at most 0.01, including interaction. The six-unit pilot uses its entire available window and retains the same screen. Full-window refinement is still mandatory.

Every-step target-volume error must stay below 5%, equivalent radius at least four grid spacings, clipping zero, boundary occupancy below 0.01, and dilution amount-conversion error at most 2e-14. No temporal alignment or post hoc relaxation is permitted. Any paired failure blocks that history's long continuation and stops the coordinator for review. The existing exact time-zero restart wrapper is retained without editing its pinned source.

Use `python -m embryo.phase_carry_context_dosage prepare`, then `preflight` and `run`, with output `outputs/phase-carry-context-dosage`. Use the dedicated qualified GTX 1080 Ti with `CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1`. The protocol pins the latest accepted delayed review, original reset review, inherited sources/inputs, this document, the new runner/tests and prepared starts. Launch success is separate from scientific acceptance. A later response assay should be considered only after reviewing this state factorial.
