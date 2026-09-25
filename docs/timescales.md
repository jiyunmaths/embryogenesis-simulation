# Signaling growth versus cleavage timescale

This experiment asks whether the weak signals in the initial time-15 run reflect insufficient amplification time, loss of growing modes as contacts change, or a failure of the default signaling mechanism to sustain nonlinear patterns. It keeps signaling symmetry breaking separate from downstream fate labels and embryo shape.

## Results of the first timescale screen

**The weak signals at time 15 were not evidence that the default mechanism cannot pattern.** With the same parameters and seed, the fully coupled simulation develops large signaling differences when followed to time 60. The frozen-graph comparisons show why cell count and available amplification time must both be considered.

| Frozen actual graph | Fastest spatial growth rate | Growth e-folding time | Additional time to activator SD 0.1 | Activator SD after 120 units |
|---|---:|---:|---:|---:|
| 4 cells | −0.202513 | No growing mode | Not reached; perturbations decay | Approximately zero |
| 8 cells | 0.069621 | 14.36 | 88.20 | 0.3236 |
| 16 cells | 0.214884 | 4.65 | 23.30 | 0.6546 |

The independently measured modal growth rates agree with their predictions within $6\times10^{-7}$. The 8-cell graph is linearly unstable but requires much longer than one cell cycle to reach the larger amplitude marker from its inherited state. In the baseline run, the interval at exactly eight cells is only 2.61 units, before asynchronous division changes the graph again. Stability at a snapshot does not imply enough time to form a mature pattern there.

| Fully coupled run | First reaches 16 cells | Activator SD first reaches 0.1 | Delay after reaching 16 cells | Activator SD at time 60 |
|---|---:|---:|---:|---:|
| Normal cycle interval, 2 | 12.375 | 36.135 | 23.760 | 0.6781 |
| Doubled cycle interval, 4 | 20.865 | 44.205 | 23.340 | 0.6642 |

Slower cleavage delays large signaling differences in absolute time, mainly because it delays access to the larger contact graph. The delays measured from first reaching 16 cells are nearly equal. At equal follow-up of 39.135 units after reaching 16 cells, activator SD is 0.6417 for normal timing and 0.6642 for slower cleavage. Thus this screen does not support a strong suppression of signaling by the default cleavage schedule. It also does not establish that slower cleavage is generally beneficial or irrelevant across parameters and seeds.

The prescribed-geometry control reproduces the complete baseline signal history to a maximum absolute error of $2.33\times10^{-13}$ for activator and $1.06\times10^{-13}$ for inhibitor. Stretching that same graph history by a factor of two through its final cleavage gives:

- First attainment of 16 cells at time 24.750, compared with 12.375 in native replay.
- Activator SD 0.1 at time 46.095: 21.345 units after reaching 16 cells, compared with 23.760 in native replay.
- At equal post-cleavage follow-up of 47.625 units, activator SD 0.6756, compared with 0.6781 in native replay.

Extra time along the earlier graphs modestly advances onset relative to the last cleavage in this control, but does not qualitatively rescue an otherwise absent signal pattern. Replay stretches all recorded pre-final-cleavage geometry, including furrowing; the coupled slower-cycle experiment changes cycle intervals only. Their birth times are consequently different, and they should not be treated as identical interventions.

At time 60, thresholded A/B/uncommitted counts are 6/10/0 for normal timing and 4/9/3 for slower cleavage. These are observations of the downstream fate variable, not persistence or irreversible-commitment tests. This experiment also does not establish a developmental shape axis.

### Numerical checks and interpretation limits

Reducing the frozen signaling outer time step from 0.05 to 0.005, which also reduces the actual SSP-RK2 substep, changes activator by at most $1.03\times10^{-6}$ and inhibitor by at most $4.94\times10^{-6}$ over the 120-unit trajectories. The graph snapshots and initial signals exactly match the corresponding division events in the earlier time-15 baseline experiment.

The largest sampled individual-cell volume errors in the coupled runs are 1.65% and 1.68%; there are no sampled overdue divisions. Maximum sampled boundary occupancies are 0.00211 and 0.00364. These diagnostics reveal no gross volume failure, but do not replace mechanical time/grid refinement or a larger-domain check. The test suite passes 51 tests, including the new replay and timing checks.

The result supports **insufficient available growth time as the main explanation of the original weak time-15 signals in this particular run**. Growing contacts do not prevent later amplification here. It leaves open robustness across seeds, sensitivity to contact pruning and grid resolution, nonlinear mode selection, and whether the resulting identities persist under perturbations.

Machine-readable evidence is generated in `outputs/signaling-timescales/analysis/analysis.json`; the comparison figure is `outputs/signaling-timescales/analysis/timescales.png`. The finer signaling analysis is in `outputs/signaling-timescales/analysis-fine`. These directories are local generated artifacts; the commands below reproduce them.

## Protocol

All primary runs use seed 7, the default $40^3$ mechanical grid, time step 0.015, a 16-cell cap, and the default Gierer–Meinhardt coefficients $\beta=2$, $D_a=1$, $D_b=20$. Partition-noise scale remains 0.001. This is a single-seed timescale screen; it is not a robustness or numerical-convergence study.

| Comparison | What changes | What it tests |
|---|---|---|
| Frozen 4-, 8-, and 16-cell graphs | Hold the actual graph at first attainment of each count; evolve the inherited signals for 120 additional time units | Whether a supported unstable mode has enough time to grow on that particular graph |
| Independent small-eigenmode check | Apply a tiny perturbation aligned with the fastest spatial mode and its chemical eigenvector | Whether measured early exponential growth agrees with linear stability |
| Coupled baseline | Continue full 3D signaling, fate, polarity, and mechanics to time 60 | Whether the original weak pattern develops with longer follow-up |
| Coupled slower cleavage | Double the mean cell-cycle interval from 2 to 4; run to time 60 | How more time between cleavage events changes the coupled trajectory |
| Prescribed-geometry replay | Double the dwell time of every recorded graph through the final cleavage, then retain the original post-cleavage durations | Isolate timing along an identical sequence of contacts, without altered signaling changing mechanics |

The slower coupled run changes the cell-cycle interval only. Cytokinesis duration, ring force, signal rates, noise scale, and mechanical relaxation rates are unchanged. Thus it is not a uniform slowing of every process. Because the mechanics responds to signaling, its shapes, division axes, and actual abscission times can differ from the baseline.

Frozen graphs are followed for longer than the coupled runs because the 8-cell graph has a much slower predicted growth rate: one e-fold takes approximately 14 time units. A 120-unit observation window can resolve growth that remains small after 60 units. This extended window changes observation time only, not the kinetics.

The baseline and slower run are compared both at the same absolute time and at the same elapsed time after reaching 16 cells. First crossings of activator standard deviations 0.01 and 0.1 are descriptive amplitude markers. They are not bifurcation points, persistence criteria, or fate-commitment thresholds. Crossings are measured on the recorded sampling grid without interpolation; equal-follow-up amplitudes use linear interpolation between samples.

## Frozen graphs and linear predictions

Stage graphs are captured immediately after first reaching exactly 4, 8, or 16 cells. These are actual weighted contacts at exact division events, rather than the later, regularly sampled graph frames in earlier reports. Other cells may still be undergoing cytokinesis. A frozen graph is consequently a controlled snapshot, not an equilibrated packing or an assertion that those contacts remain fixed biologically.

The nonlinear trajectory begins from the activator and inhibitor values present at that snapshot. It receives no further partition noise and has no subsequent cleavage. Its reference linear trajectory starts from the same departures from the homogeneous equilibrium:

$$
\begin{aligned}
\delta\mathbf{z}(0)&=(\mathbf{a}(0)-\mathbf{1},\mathbf{b}(0)-\mathbf{1}),\\
\frac{\mathrm{d}\delta\mathbf{z}}{\mathrm{d}t}
&=\begin{pmatrix}
I+D_a\Delta&-I\\
2\beta I&-\beta I+D_b\Delta
\end{pmatrix}\delta\mathbf{z}.
\end{aligned}
$$

The linear trajectory uses a matrix exponential; the nonlinear trajectory uses the existing positive SSP-RK2 integrator, with an outer time step of 0.05 and its internal rate-dependent substeps. Once perturbations become large, linear continuation is a comparison curve rather than a valid nonlinear prediction.

Separately, `verify_mode` uses a small modal perturbation and fits its logarithmic amplitude over five time units. Its growth rate is compared with the eigenvalue of $J-\lambda_k\operatorname{diag}(D_a,D_b)$. For positive rate $r$, the e-folding time is $1/r$. The initial signal standard deviation alone does not specify the amplitude projected onto that growing mode; therefore $\log(A_{\mathrm{target}}/A_{\mathrm{initial}})/r$ is only a rough timescale estimate.

## Why replay is needed

Full coupled runs alone cannot isolate the effect of spending longer on the same graphs: their geometries respond differently. The capture records the contact matrix, cell order, and both signal variables at every mechanical step, plus every division's prolongation matrix and partition jumps.

At a replay division, the parent-to-daughter map $P$ is applied in the original cell order. A fractional partition perturbation is reconstructed from each original species event:

$$
\boldsymbol\epsilon=
\frac{\mathbf{x}_{\mathrm{after}}-P\mathbf{x}_{\mathrm{before}}}
{P\mathbf{x}_{\mathrm{before}}},
\qquad
\mathbf{x}'_{\mathrm{after}}=
(P\mathbf{x}'_{\mathrm{before}})\odot(\mathbf{1}+\boldsymbol\epsilon).
$$

The division draws are reused as relative perturbations, not as the old additive jumps. This retains the original partitioning rule when the replayed parent activity differs. The division order and lobe fractions are prescribed from the recorded path. Both signals are inherited; neither is reset to equilibrium at a new cell count.

A replay with unchanged timing must reproduce the recorded coupled signals to within $10^{-8}$ maximum absolute error over the entire trace. This validates the contact timing, event order, inheritance, and numerical update before interpreting the stretched replay. Stretched replay does not predict a physically realizable mechanical trajectory and does not evolve a separate fate or polarity state.

## Reproduce the experiment

Each destination must be new. The two captures are independent and can run concurrently if resources permit:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.timescales capture \
  --output outputs/signaling-timescales/baseline --duration 60

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.timescales capture \
  --output outputs/signaling-timescales/slow --duration 60 --cycle-factor 2

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m embryo.timescales analyze \
  --baseline outputs/signaling-timescales/baseline \
  --slow outputs/signaling-timescales/slow \
  --output outputs/signaling-timescales/analysis
```

For a signaling time-step check, repeat the analysis into a new directory with `--frozen-dt 0.005`. This reduces the actual SSP-RK2 substep as well as the outer sampling step; it does not rerun 3D mechanics.

Each capture writes `config.json`, `graph_preflight.json`, `trace.npz`, `cleavage_spectra.json`, `metrics.json`, `metrics.csv`, `run.json`, and a full `final_state.npz` checkpoint. It omits sampled surface movies to keep the experiment focused on signaling and timing. The trace contains raw contact weights at every step; the configured graph cutoff is reapplied consistently during analysis.

Analysis writes frozen trajectories and their linear references, native and stretched replays, a machine-readable `analysis.json`, and `timescales.png`. First-attainment times, signal amplitude crossings, matched follow-up, modal checks, and mechanical diagnostics are available without rerunning the 3D simulations. Generated output directories remain excluded from Git; the protocol and results interpretation are versioned here.

The relevant implementation is [timescales.py](../embryo/timescales.py); mathematical assumptions are in [graph_signaling.md](graph_signaling.md). Tests verify full signal reconstruction across actual divisions, volume-weighted inheritance under reused fractional noise, analytic eigenmode growth, unchanged post-cleavage follow-up under time stretching, and refusal to overwrite experiments.
