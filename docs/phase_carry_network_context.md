# Moving chemical-network context reset with phase carry

## Question and design

The completed [carry exchange/response study](phase_carry_exchange_response_assessment.md) finds donor-nearer responses in both recipients within each of three existing histories. This test asks whether changing the surrounding chemical state changes a recipient's maintained state or perturbation response while the tissue continues to move.

Use the qualified **fine-level exchanged endpoints at t=450**, chi=0.35, in histories 7, 8 and 9. Their previously selected pairs remain fixed: 22/23, 19/20 and 29/30. Preserve the actual exchanged geometry, polarity, volumes, cell IDs, lineage, random streams and inherited float64 phase carry. Both new timestep levels start from this identical fine physical preparation. No new history, fresh zygote development, division or supplied identity class is introduced.

Compare three chemical contexts per history:

- **Sham:** retain all exchanged concentrations exactly; share its unperturbed control between the two pulse targets.
- **Reset surrounding chemistry for the first recipient:** keep that recipient's activator and inhibitor exactly unchanged; set every other cell to its same-ID concentration in the qualified fine untouched t=450 endpoint.
- **Reset surrounding chemistry for the second recipient:** apply the same rule while preserving the other recipient instead.

For recipient i and either species s, the reset is

$$
c_{s,j}^{\mathrm{reset}} =
\begin{cases}
c_{s,i}^{\mathrm{exchanged}}, & j=i,\\
c_{s,j}^{\mathrm{untouched}}, & j\ne i.
\end{cases}
$$

The reset changes **all other cells**, not just immediate contact neighbors. Record its initial surrounding-cell volume-weighted log RMS; require more than 0.05 so an ineffective intervention is explicit. It externally adds/removes recorded amounts,

$$
\delta N_s=\sum_{j\ne i}V_j^{\mathrm{exchanged}}
\left(c_{s,j}^{\mathrm{untouched}}-c_{s,j}^{\mathrm{exchanged}}\right).
$$

Use the recipient background's actual saved GPU-measured volumes for this accounting. It is not an amount-preserving exchange. Changes in total amounts and their spatial distribution are both part of this collective context intervention; it cannot isolate topology or distinguish spatial redistribution from bulk dosage by itself.

Each reset receives its own unperturbed moving control and a separate **immediate -10% activator pulse in the preserved recipient**. The sham background has its control and separate pulses in both selected cells. Pulse amount removal is recorded separately from the reset. There is no relaxation before the pulse, clamping, reservoir, fixed chemical graph or permanent target constraint: all cells and both species subsequently evolve freely.

## Numerical implementation and reuse

The model and kernels stay unchanged: resident PyTorch arrays/matrix operations plus custom CUDA mechanics and spatial geometry/polarity, float32 phase fields with inherited float64 phase-update residuals, float64 chemistry/kinetics, live conservative transport and geometric dilution. Parameters remain beta=2, D_a=0.02, D_b=0.55, c_gamma=0.25, c_A=0.35, chi=0.35, grid 72 cubed, extent 2.24 and interface width 0.085.

Use dt=0.00375/0.001875, with observations every 0.15 and atomic carry-preserving checkpoints every three units. Retiming changes only the discrete clock, dt, steps and save interval. Preserve the recipient's initial state and all nonchemical arrays exactly.

The nine exact fine sham controls/pulses are reused from the completed parent response study. Reuse requires **identical complete initial checkpoint arrays, dtypes and metadata**, plus the parent's original protocol, completed history/checkpoint, native-prefix and pilot/full-field hashes. Matching parameters or approximate fields alone cannot authorize reuse. Their full 60-unit numerical/physical qualification remains applicable. Coarse sham paths are new because they now start from the fine physical state, rather than the older coarse t=450 state.

There are **42 paths nested within three histories: 33 new and nine reused**. Nine new coarse sham paths and 24 reset control/pulse paths are new. Targets, chemical contexts, timesteps, controls and pulses are nested interventions, not independent developmental replicas.

Every new context requires the inherited strict native/carry prefix, eight independent NumPy chemical steps with dilution accounting, and four exact GPU carry-restart steps. Six-unit coarse/fine pilots gate long continuation separately within each history. A failed gate prevents that history's long runs and remains unresolved; scientific pattern loss or context-sensitive behavior is not a numerical failure.

Before the initial launch, an additional full-size preflight checks the four reset/pulse contexts in history 7: both preserved recipients at both timesteps. These 72-cubed checks supplement the small implementation/GPU tests. They do not qualify the remaining contexts or replace the paired pilots.

Both pilots and full runs preserve their exact field/history snapshots. The full horizon is 60 units, t=450 to 510. Completion does not automatically constitute acceptance.

## Measurements and decisions

Keep chemical-state persistence and response sensitivity distinct:

1. **Contrast:** S is across-cell SD of natural-log activator. Late contrast requires S>0.1 at every saved observation in the final 24 units. Whole-window retention and the first downward crossing are separate diagnostics.
2. **Recipient chemical-state shift:** compare the reset control with the sham moving control at the same times, using the two-species target log RMS. Report the final shift and maximum over the final 24 units. A late maximum above 0.01 is a context-effect screen, not a cell-type label, stationarity test or proof of a sustained plateau.
3. **Response shift:** use signed target log(pulse/control)/|ln(0.9)| for both species, with each context's own moving control. Compare reset versus sham using the full-window two-species time RMS; a shift above 0.01 is a response-effect screen. An immediate-pulse effect can reflect a transient context challenge even if the background state later recovers.
4. **Donor/destination references:** retain the qualified fine untouched t=450 response references from the parent. Both timestep levels compare with these same fixed reference waveforms. References must be separated by more than 0.01 to classify the nearer one; otherwise report unresolved. Nearest-reference status is descriptive, not donor equivalence or identity autonomy.
5. **Response/recovery:** record activator/inhibitor peak gains, activator area, other-cell response, and target/whole-network recovery to within 10% of initial displacement, requiring at least 24 later observed units. Nonrecovery is right-censored. Pilot recovery uses the same rule and is therefore censored at six units.

Require raw matched-time coarse/fine errors no larger than the inherited gates: log chemistry 0.01, polarity absolute 0.01, relative axis 0.01, relative volumes 0.005, relative transport 0.01, final phase absolute 0.02. Contrast/retention/loss presence must agree and observed loss-time differences must be at most 0.3.

Additionally require pulse-normalized response error <=0.01 across **all cells, both species and all saved times**, relative target activator-area error <=0.02, matching recovery presence and times within 0.3, and agreement of effect/reference classifications. Absolute discrepancies in late target-state shift and reset-versus-sham response shift must each be <=0.01. Do not align trajectories or relax thresholds after observing results.

Every-step physical gates remain target-volume error <5%, minimum equivalent radius >=4 grid spacings, no clipping and dilution amount-conversion error <=2e-14. Sampled boundary occupancy must remain <0.01. Chemistry and geometric diagnostics must be positive and finite where required.

## Interpretation limits

A qualified state or response shift is evidence that the chemical context affects this recipient in the coupled model. Retained donor-nearer behavior despite this particular reset shows resilience to the tested challenge, not autonomy from every environment. A transient pulse effect is separate from changed maintenance. The intervention also changes chemical feedback to mechanics, so it does not isolate signaling from later mechanical mediation.

Report outcomes by three existing histories. Both selected cells are highly exposed in these preparations; this is not an inner/outer transplantation test. No new endpoint bistability assay, inherited fate, broad parameter map, fresh carry developmental trajectory, spatial convergence or general geometric transport closure is qualified. This t=450–510 test does not retroactively accept older neighbor/reset, reservoir or response-amplitude claims. Preserve the parent evidence and original kernels unchanged.

The runner is `python -m embryo.phase_carry_network_context prepare`, followed by `preflight` and `run` on the qualified GPU. The run validates context gates and pilot comparisons before allowing long trajectories, saves resumable carry checkpoints, and records scientific outcomes separately from numerical acceptance.
