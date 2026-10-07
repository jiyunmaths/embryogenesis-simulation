# Carry-corrected moving state exchange and response

## Question and starting evidence

Do established chemical states carry response behavior into a different geometric context, or does the recipient's surroundings determine its response? The [completed maintenance test](phase_carry_maintenance_assessment.md) qualifies mature developed states in histories 7–9. This experiment uses their **fine-timestep, chi=0.35 endpoints at t=390**. It does not transplant an equilibrium from an older no-carry assay or generate a new developmental history.

Each history supplies one identical physical starting state for both new timesteps, 0.00375 and 0.001875. Preserve the phase-field rounding residual, geometry, polarity, lineage, chemical state, and random streams. Only the numerical clock/configuration is retimed. PyTorch owns GPU arrays and matrix operations; the existing custom CUDA mechanics and spatial geometry/polarity kernels remain unchanged. Chemistry, current cell volumes, dilution and conservative contact conductances evolve in every scientific trajectory. The parameters remain beta=2, D_a=0.02, D_b=0.55, c_gamma=0.25, c_A=0.35 and chi=0.35 on the 72-cubed grid with extent 2.24.

## Selection and intervention

Select one pair per history before observing any intervention outcome. Standardize initial cell exposure, conservative transport exit rate -Delta_ii, and measured volume across the sixteen cells. A contrasting pair has Euclidean distance in these three standardized features at least the median across all unordered pairs. Among contrasting pairs, choose the greatest absolute log-activator difference, with deterministic cell-ID tie breaking. Require initial exchanged-pair chemical log distance above 0.1. Record all features, context distances, and selected cell IDs. This deliberately selects contrasting chemistry and context; it does not estimate how often such pairs arise.

Compare **unexchanged** and **fresh conservative exchange**. Exchange both activator and inhibitor concentrations between the pair, then rescale each species within that pair to conserve its measured initial amount. For species s and selected cells i,j:

```text
alpha_s = (V_i c_i^s + V_j c_j^s) / (V_i c_j^s + V_j c_i^s)
new c_i^s = alpha_s c_j^s
new c_j^s = alpha_s c_i^s
```

Other cells remain unchanged. V denotes the actual saved GPU-measured initial volume; c denotes positive concentration/activity. The common scale is necessary when cell volumes differ. This is a species-amount-conserving chemical transplant, rather than an exact concentration swap. The exchange does not move cell surfaces or polarity and does not introduce identity labels. Verify a sham assignment is identical to the unexchanged preparation; verify that chemical editing and retiming leave every other checkpoint field, including the float64 phase carry, unchanged. Subsequent reactions can change chemical amounts.

## Two moving stages

1. **Exchange relaxation:** co-evolve both backgrounds from t=390 to 450. Record every 0.15 units. Compare pair concentrations with the unexchanged moving destination and with the conservatively transferred unexchanged reference. Use fixed initial measured volumes for descriptive distances; actual evolving volumes govern transport and dilution. Report late-window contrast, continuous retention, first observed contrast loss, and continuous pair-reference distances separately.
2. **Response:** use each background's own completed t=450 geometry, chemistry, polarity and carry. Continue an unperturbed control and a separate -10% activator pulse in each selected cell to t=510. Each pulse is compared with its own background's control. The pulse removes activator externally, with its amount explicitly recorded; it is not an amount-conserving exchange. No chemical re-equilibration or carry reset occurs between stages.

This gives **12 exchange-stage paths and 36 response-stage paths**, 48 scientific paths nested within **three existing developmental histories**, with no new histories. Negative pulses and fresh exchange are the prespecified scope. Positive pulses, frozen pre-relaxed exchange, other polarity contrasts, neighbor reset, and reservoir tests remain separate work.

## Measurements and scientific decisions

Chemical contrast is S=SD(ln a) across cells. Sustained contrast means S>0.1 at every saved observation over the final 24 units; whole-window retention is reported separately. A loss is an admissible scientific outcome.

Chemical distance is the volume-weighted RMS difference in log concentrations over the two species. Pair distances use only the selected cells. Divide the maximum late-window distance by the initial exchange separation. A ratio below 0.1 is described as destination-like or transferred-like; otherwise report reorganized. These are moving similarity criteria, not stationarity or cell-type classifications. The donor reference at each observation is the amount-preserving exchange of the unexchanged control using fixed initial masses; it is a diagnostic reference, not another simulated trajectory.

The signed target-cell response is:

```text
R_s(t) = ln(c_pulse^s(t) / c_control^s(t)) / abs(ln(0.9))
```

Compare the two-species response waveform with the same-age unexchanged donor and destination waveforms using time-RMS distance. Report both distances and their signed difference. A nearest reference is informative only when donor/destination reference separation exceeds 0.01. Donor-nearer behavior indicates transferable behavior within the interacting tissue; it does not imply donor equivalence, cell autonomy, biological commitment or preservation of an unchanged concentration vector. If references become indistinguishable, report unresolved rather than force a label.

Also record response integrals, target and other-cell gains, and sustained target/network recovery. Recovery requires all subsequent samples within 10% of the initial displacement and at least 24 units of subsequent observation. Nonrecovery is right-censored at 60 units.

## Numerical gates fixed before launch

Every exact starting context, including all response controls and pulses, requires a 0.6-unit native CPU/carry-GPU comparison before its long continuation. Retain the existing strict limits: chemical log error 1e-5, polarity absolute error 1e-5, relative volume/transport error 1e-5, relative axis error 1e-4, and final phase-field absolute error 2e-5. Independently evaluate eight chemical steps with NumPy, requiring log error below 1e-11 and dilution amount error at most 2e-14. Check exact carry-preserving restart for four steps. This is a short context check; the existing full-horizon backend comparisons and float64 mechanics references remain supporting evidence at their tested scopes.

Run paired six-unit pilots before the full 60-unit stages. An exchange-stage pilot must pass its timestep gates before advancing that history to long exchange runs. Response pilots include both matched controls and both targets in both backgrounds; their waveform checks must pass before long responses start. A failed numerical gate prevents the corresponding scientific continuation and remains recorded.

Exchange-stage timestep checks use raw matched times: maximum log-concentration difference 0.01, polarity difference 0.01, relative axis difference 0.01, relative volume difference 0.005, relative transport difference 0.01, final phase-field absolute difference 0.02, pair-reference ratio difference 0.05, matching contrast/retention/state-similarity outcomes, and observed loss times within 0.30. Pilots use all six units for their contrast/ratio window.

Response-stage checks retain the existing response-specific limits: maximum normalized pulse-minus-control waveform discrepancy across all cells/species/samples 0.01, relative target-activator response-integral discrepancy 0.02, matching recovery presence and recovery times within 0.30, raw chemical log discrepancy 0.01, and matching informative/nearest-reference decisions. Also require the same mechanical/transport/field limits as exchange. Pilot recovery is censored by its short window; a pilot pass does not qualify full response measurements.

Every step checks positive finite chemistry and geometry, per-cell target-volume error below 5%, minimum equivalent radius at least four grid spacings, zero clipping, and dilution accounting at most 2e-14. Sampled boundary occupancy must stay below 0.01. Checkpoints are atomic every three model units and include float64 phase carry; interruption/resume must not repeat an exchange or pulse. Numerical acceptance and scientific donor/destination outcomes are separate decisions.

## Provenance, scheduling and interpretation

Preparation verifies the completed maintenance protocol, summary, assessment, trajectories, endpoints and checkpoint hashes. Freeze this protocol document, new runner/tests, unchanged inherited scientific sources, selected sources and prepared inputs before launch. Response handoffs are materialized only from complete, passing exchange results, then independently hashed. Preserve all earlier studies and unresolved old-method measurements.

A single worker uses the dedicated GTX 1080 Ti, which was the measured efficient scheduling choice. Histories and timestep pairs are scheduled in stages. Root/child statuses record elapsed time, completed scientific paths and gate decisions. An outcome is assessed by history, with backgrounds, targets, pulses and timesteps treated as nested interventions. This is a carry-corrected mature exchange/response experiment on three existing histories, not a new zygote-to-identity result or a broad promotion of the old ledger.

```bash
CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python -m embryo.phase_carry_exchange_response prepare

CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python -m embryo.phase_carry_exchange_response run

python -m embryo.phase_carry_exchange_response assess
```

The new output directory is `outputs/phase-carry-exchange-response/`. Preparation and scientific execution are distinct. Runtime estimates will use the completed carry maintenance throughput and then the new paired pilots; estimates are not acceptance criteria.
