# Moving-geometry cell-response pilot

Current scope: the seed-7 pilot, its response refinement, fifteen moving exchange-response runs, and six targeted refinement runs are complete. The replication on seeds 8 and 9 is running through the accepted PyTorch/custom-CUDA backend after all six new-context CPU/GPU checks passed. The following protocols retain their distinct starts, horizons, and inference limits; [the current replication](#replication-across-developmental-histories) is conditional on mature patterned basins.

## Question

Do the response differences found in the frozen assay persist when cell shape, polarity, transport, and dilution coevolve? The experiment compares a local chemical pulse with a matched unperturbed **moving** control, not with its initial state. This subtracts the background's ongoing development.

## Prespecified pilot

- Source: seed 7, full-coupling time-150 endpoint (`outputs/feedback-survival/switch_on/latest_state.npz`).
- Two chemical starting backgrounds: the validated patterned and uniform equilibria on that graph, taken from the preceding frozen bistability experiment. The equilibrated chemistry replaces raw time-150 chemistry; geometry, polarity, cell IDs, clocks, material parameters, and random streams retain their source values.
- Two target cells: minimum and maximum initial patterned activator, with ties resolved by cell ID. The same IDs are targeted on the uniform background. This selection uses initial state, not response outcomes, and is not a classifier assigning two cell types.
- Interventions: multiply only the target activator by 0.9 or 1.1 at elapsed time zero. Inhibitor and other concentrations are unchanged. These are external chemical boluses with unequal absolute amounts across backgrounds; they are not conserved redistribution.
- Four pulse runs and one unperturbed control per background: ten moving continuations total.
- Duration: 60 model-time units, from source time 150 to 210. Grid 72³, source timestep 0.0075, observation interval 0.15. Two processes run concurrently. No further cleavage occurs at the sixteen-cell cap.
- Matched frozen references start from precisely the same initial chemistry and operator, sampled over the same 60-unit interval. Their DOP853/Radau comparisons must pass maximum log error 1e-5 before the moving study is launched.

The pilot focuses on one geometry with full mechanical coupling. It neither repeats the entire three-history/all-cell frozen screen nor isolates polarity from tension and adhesion. Uniform versus patterned backgrounds change the entire network's chemistry. Controls are scheduled first so later pulse comparisons have completed reference trajectories.

## Outcomes and safeguards

Target activator response integral, induced inhibitor peak, and response in other cells use the same log-ratio definitions as `cell_response.py`. Moving changes are measured relative to the matching moving control; frozen changes use the matching frozen control. Initial measured volumes provide fixed weights so changes in weighting do not masquerade as chemical response. Recovery requires deviation to remain below 10% of the initial log displacement for all remaining samples, with at least 24 subsequent time units observed. This preserves the previous hold criterion, but changes the total observation horizon from 240 to 60 and the sampling interval to 0.15. Recovery is therefore censored beyond the supported observation window.

We also report pulse-control centroid displacement, polarity difference, relative aggregate axis-ratio difference, and relative transport-operator difference at the endpoint. Shape and chemical effects are reported separately. Static-equilibrium classification is deliberately not applied to the moving trajectories: an endpoint difference does not prove another attractor or cell identity.

Every-step per-cell volume error must remain below 5%, minimum equivalent radius at least four grid spacings, and phase-field clipping zero. Chemistry must remain positive and finite. Sampled boundary occupancy must remain below 0.01. Quality failures remain failures and are not interpreted as biological state loss. The earlier survival timestep check does not validate this new pulse experiment; response-specific refinement remains necessary before a strong quantitative conclusion.

## Reproducibility and operation

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.cell_response_moving prepare
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.cell_response_moving run --workers 2
```

Preparation requires a fresh output directory. The default is `outputs/cell-response-moving/`. Protocol, source code, starting checkpoints, initial chemistry, and frozen references are hashed before launch. Each job retains a status, history, and atomic restart checkpoint every three model-time units. Re-running `run` verifies immutable inputs and resumes the matching checkpoint; it does not repeat the pulse. Completed jobs are reused only if their protocol hash matches. An aggregate comparison is generated only after all jobs finish with passing quality checks.

Eight focused tests passed before preparation, including preservation of nonchemical state under the pulse, subtraction of a drifting matched control, and exact restart agreement on a small-grid continuation. These are implementation checks, not convergence evidence for the production runs.

## Completed pilot assessment

All ten jobs completed through elapsed time 60, with eight pulse trials and two matched controls. Input and code hashes match the prepared protocol. All quality checks passed: maximum volume error 1.143%, minimum equivalent radius 5.108 grid spacings, zero clipping, and maximum sampled boundary occupancy 3.004e-10.

| Background | Target recovery time | Target response integral, normalized by log pulse |
| --- | --- | --- |
| Patterned | 3.30–3.75 | 1.38–1.64 |
| Uniform | 7.80–9.00 | 3.29–4.06 |

Every pulse recovered under the prespecified sustained-recovery criterion. Moving versus frozen target response integrals differed by at most 2.26%. The patterned control retained chemical contrast (standard deviation of log activator 1.182 to 1.274); the uniform control ended near uniform (7.03e-8). The high-activator patterned target induced substantially more inhibitor and neighbor response than the low-activator target. These are state-dependent responses within a collective chemical background, not evidence of autonomous cell types: only one history and two selected cells were tested, and relative pulses deliver different absolute amounts.

Final pulse-control geometric differences were small, with centroid RMS differences no larger than 3.48e-5 model-length units and relative operator differences below 7.08e-5. Such residuals do not establish mechanical memory. The two chemical backgrounds were preconditioned equilibria on a mature geometry, not independently developed cell identities.

Results: `outputs/cell-response-moving/comparison.json`; complete trajectories and numerical audits reside in each job directory. The original protocol and results remain unchanged.

## Completed response-specific timestep refinement

Repeat all ten jobs at timestep 0.00375, preserving the same physical starting state, random streams, target IDs, pulses, 60-unit horizon, and 0.15 observation interval. Compare each pulse against its own matched moving control at each timestep. Halving the timestep changes the numerical integration; it does not add a developmental replicate.

Before launch, the acceptance thresholds are fixed at: maximum difference between coarse and fine pulse-minus-control log responses below 1% of the initial log pulse; target response-integral relative difference below 2%; target and network recovery classification unchanged and recovery times within 0.30 units (two sample intervals); maximum raw chemical log discrepancy below 0.01. The waveform condition checks every sampled cell and both chemical species. These are practical agreement criteria for one timestep halving, not a proof of asymptotic convergence. Tiny mechanical residuals are excluded from the validated claims.

```bash
python -m embryo.cell_response_moving_refinement prepare
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.cell_response_moving_refinement run --workers 2
```

The refinement batch completed and passed all declared agreement criteria. The maximum normalized response discrepancy was 5.03e-5, the maximum relative response-integral discrepancy was 2.37e-4, and sampled target/network recovery times were unchanged. The maximum raw chemical log discrepancy was 1.16e-4. The separate directory `outputs/cell-response-moving-refined/` preserves the pilot; `refinement.json` records the decision. This validates one timestep halving for the original t=150 pulse assay, not every later starting state.

## Completed moving exchange-response study

All fifteen seed-7 continuations completed from t=210 to 270 at dt=0.00375, on the 72³ grid with extent 2.24 and sixteen cells. Each of three backgrounds—unexchanged, fresh exchange, and pre-relaxed exchange—has its own moving control and ±10% activator pulses in cells 27 and 20. Pre-relaxed chemistry was transplanted onto the source t=150 geometry and then co-evolved to t=210. These are matched interventions in one developmental history.

We compare the signed, two-species target log-response waveform against the matching control, normalized by the absolute log pulse. Untouched same-age donor and destination responses provide references. All eight exchange comparisons favor the donor:

| Recipient | Background | Pulse | Distance to donor | Distance to destination |
|---|---|---|---:|---:|
| 27 (high-state recipient) | Fresh | −10% | 0.11234 | 0.18710 |
| 27 | Fresh | +10% | 0.08664 | 0.17379 |
| 27 | Pre-relaxed | −10% | 0.11236 | 0.18711 |
| 27 | Pre-relaxed | +10% | 0.08665 | 0.17379 |
| 20 (low-state recipient) | Fresh | −10% | 0.00698 | 0.16624 |
| 20 | Fresh | +10% | 0.00791 | 0.16571 |
| 20 | Pre-relaxed | −10% | 0.00698 | 0.16624 |
| 20 | Pre-relaxed | +10% | 0.00791 | 0.16571 |

The original reference separation is approximately 0.165 for both pulse signs, exceeding the declared 0.01 resolution threshold. Fresh and pre-relaxed responses are very similar. The low-state recipient closely resembles its donor; the high-state recipient is closer to its donor but retains appreciable differences. For negative pulses, target recovery shifts from 3.15 to 6.45 units in cell 27 and from 4.05 to 2.85 units in cell 20 after exchange. These observations support transferable response behavior within the coupled system; nearest-reference classification does not establish donor equivalence, autonomous identity, or independent biological replication.

Every run passes the numerical-quality screen: maximum per-cell volume error 1.132%, minimum equivalent radius 5.1099 grid spacings, zero clipping, and maximum sampled boundary occupancy 5.03e-10. Controls preserve strong chemical contrast through t=270. Results and complete histories are in `outputs/cell-exchange-response-moving/`; `comparison.json` contains response metrics and donor/destination distances. The native-backend transition manifest records the accepted backend and source provenance.

## Targeted exchange-response timestep refinement

Repeat the unexchanged and fresh-exchange controls and −10% pulses in both cells at **dt=0.001875**: six runs, three workers with two native threads each. Start from exactly the same t=210 physical states and chemistry, retain the 60-unit horizon and 0.15 observation spacing, and compare with the completed dt=0.00375 continuations. This holds exchange formation fixed; it does not refine the earlier t=150–210 formation trajectory.

Acceptance criteria are unchanged from the original response refinement: normalized pulse-minus-control waveform discrepancy ≤0.01 across every cell/species/sample; response-integral relative discrepancy ≤0.02; recovery classification unchanged with time discrepancy ≤0.30; and maximum raw chemical log discrepancy ≤0.01. Also require resolvable donor/destination references and unchanged nearest-reference classification at each timestep, using that timestep's own untouched references. The assessment fails if any criterion fails.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.cell_exchange_response_refinement prepare
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.cell_exchange_response_refinement run
python -m embryo.cell_exchange_response_refinement assess
```

All six fine-timestep runs completed and the targeted agreement check passed. Maximum normalized response discrepancy was 6.25e-5 (limit 0.01), maximum relative response-integral discrepancy was 3.22e-4 (limit 0.02), and raw chemical log discrepancy was 1.64e-4 (limit 0.01). All sampled target/network recovery times were unchanged. Both exchanged cells remained closer to the donor reference. The separate directory `outputs/cell-exchange-response-refined/` preserves the original study; `refinement.json` records the decision.

The targeted check excludes positive pulses, pre-relaxed exchange, spatial refinement, and independent developmental histories. Fifteen focused tests passed before preparation, including changed transfer classification, loss of recovery, invalid/misaligned evidence, retiming, and native-worker integration. The [full-horizon GPU-backend validation](model.md#resident-gpu-backend-validation) subsequently passed; the following replication uses the gated mature GPU backend.

## Replication across developmental histories

The running experiment repeats the moving exchange-and-response assay on **seeds 8 and 9**, with the completed seed-7 study retained as a historical reference. It asks whether transferable response behavior occurs in more than one independently developed geometry. Results are pending; this protocol does not assume that either new history will reproduce the seed-7 outcome.

Both sources are existing sixteen-cell, direct-feedback, polarity-enabled checkpoints at t=150 from the completed survival study. Their frozen-endpoint bistability and conservative-exchange assays must have passed. Their physical parameters match the accepted GPU regime. Retiming from dt=0.0075 to **0.00375** preserves geometry, chemical state, polarity, cell IDs, physical age, and random streams; replacing chemistry with the specified equilibrium/exchange state is a separate experimental intervention.

| Stage | Starting state | Runs per history | Time interval |
|---|---|---:|---|
| Moving formation/retention | Untouched patterned equilibrium, fresh conservative exchange, pre-relaxed exchange, all on the same source geometry | 3 | t=150–210 |
| Moving pulse response | Each formation endpoint, preserving its own evolved geometry and chemistry; one control and ±10% activator pulses in both selected cells | 15 | t=210–270 |

This gives **36 GPU continuations across two additional histories**, plus six short backend checks. The minimum/maximum initial patterned-activator pair is selected before observing responses, with cell-ID tie breaking: seed 8 uses IDs **19 and 20**, and seed 9 uses **26 and 30**. Both species are exchanged, with a species-specific rescaling of the pair that conserves each species' amount despite unequal cell volumes. The pre-relaxed state is the endpoint of the earlier 1200-unit conservative frozen exchange assay. No prescribed identity labels or fate attractors are introduced.

All six new history/background starting contexts first undergo matched **0.6-unit native CPU/GPU continuations**, sampled every 0.15 units. They use the accepted backend criteria: maximum chemical log error ≤1e-5, polarity absolute error ≤1e-5, relative volume/operator error ≤1e-5, relative axis-ratio error ≤1e-4, and endpoint phase-field error ≤2e-5. **All six checks passed before the formation runs started.** A failed check would prevent the entire scientific batch from starting. This supplements the completed full-horizon seed-7 validation; it does not establish full-horizon equivalence on every new history. PyTorch owns resident arrays and matrix operations; custom CUDA computes mechanics and spatial geometry/polarity. A single worker uses the dedicated GTX 1080 Ti.

Formation outcomes and response outcomes are assessed separately. Formation reports late chemical contrast, pair ordering, and normalized distance to the moving untouched destination and conservatively exchanged donor references. Initial pair log distance must exceed 0.1 to be informative; late destination/donor ratios below 0.1 indicate corresponding state likeness. A state can reorganize and still retain donor-nearer pulse behavior. Fixed initial measured volumes weight chemical distances, while actual evolving volumes continue to govern transport and dilution. Independently checked DOP853/Radau frozen references provide a separate comparison for every formation background.

For initial measured-volume weights, state distance is

$$
d_{\mathrm{state}}(c,\widetilde c)=
\sqrt{\frac{\sum_i V_i^0\sum_{s\in\{a,b\}}
(\log c_i^s-\log\widetilde c_i^s)^2}{2\sum_i V_i^0}}.
$$

The sum is restricted to the pair for pair outcomes and uses all cells for global outcomes. Late pair ratios use the maximum distance over the last 24 formation units, divided by initial fresh-exchange pair separation. The donor reference is the conservative exchange of the untouched moving control at each observation, using the fixed initial masses; it is a diagnostic reference, not a second simulated trajectory. Chemical amount conservation at the initial intervention does not imply conservation during subsequent reactions or equivalence of the pre-relaxed transplant's total amount.

Response analysis retains the seed-7 definitions: each pulse is compared with its own unperturbed moving control, using signed two-species log-response waveforms normalized by the absolute log pulse. Untouched same-age donor/destination responses provide references; their separation must exceed 0.01 before nearest-reference classification. Recovery requires a sustained deviation below 10% of the initial displacement with at least 24 subsequent units observed. All moving runs retain every-step volume, radius, clipping, positivity, and sampled boundary checks. Observations remain 0.15 units apart, with restart checkpoints every three units.

The [response equations](cell_response.md#current-response-definition-and-moving-history-replication) explain signed normalization and waveform distance. Each continuation uses the first-order coupled splitting described in [the model method](model.md#coupled-integration-and-computing-paths). A passing backend comparison therefore does not turn this into a higher-order solver or supply missing spatial/time-convergence evidence.

The coordinator hashes sources and inputs, preserves prior studies, rejects incomplete/misaligned evidence, and materializes pulse sources only from completed formation endpoints whose checkpoints match their histories. It aggregates the two new histories separately and retains the seed-7 reference. Fifty-one focused checks passed before preparation, including conservative exchange with unequal volumes, changed formation/response outcomes, uninformative pairs, malformed evidence, endpoint handoff, failure gating, and a real GPU adapter restart. All six frozen-reference integrations passed independent solver checks; the largest chemical log discrepancy was 5.06e-10, below the 1e-5 limit.

```bash
CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python -m embryo.exchange_response_histories prepare

CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python -m embryo.exchange_response_histories run

python -m embryo.exchange_response_histories assess
```

The default output is `outputs/exchange-response-histories/`. The root status identifies the current stage; child statuses report individual trajectory progress. Saved passing prefixes and completed jobs are verified on restart, and unfinished GPU jobs resume without repeating a pulse. `formation_comparison.json` is written for each history after formation, and the aggregate `comparison.json` is generated only after all scientific jobs pass their numerical screens.

**Interpretation limits:** this is conditional on previously formed mature patterned basins. Equilibrated or exchanged chemistry is transplanted at t=150; it does not test spontaneous identity formation from a new zygote or estimate pattern-formation frequency. Three developmental histories including seed 7 are a pilot. Backgrounds, cells, and pulse signs are nested interventions, not independent developmental replicas. Nearest donor behavior does not demonstrate donor equivalence or autonomous identity. Exchange-formation timestep refinement, new-history response refinement, spatial convergence, and broader developmental replication remain separate tests.
