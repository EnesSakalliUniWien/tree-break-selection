# Eigenspace benchmark and calibration reassessment — 2026-09-28

The complete-eigenspace repair removes the demonstrated basis dependence, but
the current methods remain miscalibrated in several tested settings. The matched
benchmark shows both improved and worsened recovery. Independent null simulation
shows large false-partition rates even when the optional adaptive Wald rank rule
is disabled.

## Matched benchmark after the repair

The same 11 `method_proof` datasets, seeds, 47 spectral blocks, eight tree
inference methods and two geometries were rerun: **752 attempts**. Inputs,
retained spectra, blocks, package versions and settings match the earlier run.

| Quantity | Result |
|---|---:|
| Successful partitions, before and after | 525 |
| Failed tree construction, before and after | 170 |
| Unsupported calibration, before and after | 57 |
| Changed rooted trees among available pairs | 18/582 |
| Changed partitions among successful pairs | 33/525 |
| Changed partitions on unchanged trees | 22 |
| Single-cluster outputs | 436 → 459 |
| Generator-label ARI increased / decreased | 15 / 15 runs |

All raw block geometries and 44/47 diffusion geometries are bitwise identical.
The three changed diffusion geometries have disconnected graphs with leading
eigenvalue multiplicities 4, 2 and 6. Retaining those complete leading groups
changes the selected directions. Changes on unchanged trees occur downstream
of tree inference.

The one-population below-BBP example changed from 10/36 multicluster outputs to
0/36. Those are repeated analyses of **one dataset**, not 36 independent null
experiments. Above-BBP outputs also moved toward one cluster; collapse alone is
not improved recovery. See the [matched benchmark report](/Users/berksakalli/Projects/kl-te-cluster/reports/eigenspace_reassessment_20260928/spectral_blocks/README.md)
and [comparison summary](/Users/berksakalli/Projects/kl-te-cluster/reports/eigenspace_reassessment_20260928/spectral_blocks/comparison_summary.json).

## Complete-pipeline null calibration

The frozen evaluation generated **4,200 independent datasets**: 500 null and
200 signal datasets in each of six families. Five methods share each exact
dataset; iid Bernoulli also has a prespecified alpha sensitivity comparison.
Every replicate refits geometry, tree, configured NNLS branch lengths,
covariance/PCA/rank, edge selection, empirical inflation, sibling BH and
traversal. This is an evaluation of the **current implementation**; it is not a
matched before/after estimate of full-pipeline calibration change.

At edge alpha .001 and sibling alpha .01, the table gives final `K > 1` events
per completed one-population null dataset:

| Null family | Average linkage | Complete linkage | Neighbor joining | Hamming diffusion | Adaptive diffusion + NNLS |
|---|---:|---:|---:|---:|---:|
| Independent Bernoulli | 404/499 | 123/499 | 295/500 | 165/500 | **316/500** |
| Dependent Bernoulli | 496/500 | 490/500 | 432/500 | 439/500 | 491/500 |
| Categorical | 357/500 | 96/500 | 238/500 | 226/500 | **263/500** |
| Gaussian, identity covariance | 8/500 | 193/500 | 315/500 | unavailable | unavailable |
| Gaussian, AR1 covariance | 54/500 | 46/500 | 427/500 | unavailable | unavailable |
| Gaussian, 64 features / 40 samples | 0/500 | 0/500 | 0/500 | unavailable | unavailable |

The dependent Bernoulli family copies features and intentionally stresses the
independent-feature covariance model. Independent Bernoulli and categorical
results already demonstrate failure without that stressor. Native continuous
inputs are incompatible with the registered Hamming diffusion methods; their
unavailable results are not zeros.

For adaptive diffusion + NNLS, the independent Bernoulli rate is **.632** with
pointwise Wilson 95% interval **[.589, .673]**; categorical is **.526**, interval
**[.482, .569]**. Twenty of the 24 primary null cells with available outputs
exceed the .01 diagnostic target after Holm adjustment across all 34 planned
global-null method/alpha cells. This compares a declared whole-pipeline
criterion; it does not assume that per-depth BH already promises whole-tree
error control. An unflagged cell does not establish general calibration.

The selected root p-values also have excess lower-tail mass. In independent
Bernoulli data, adaptive NNLS has finite root-null p-values on 482/500 datasets;
280/482 are at most .01 (rate .581, Wilson interval [.536, .624]). In categorical
data the count is 240/481 (rate .499, interval [.454, .543]). These conditional
root tails use one root per independent dataset; they are not a pooled-node
binomial calculation. The root-tail Holm family includes all 272 planned cells.

![Complete-pipeline null partition rates](/Users/berksakalli/Projects/kl-te-cluster/reports/eigenspace_reassessment_20260928/null_partition_rates.png)

## Signal recovery and calibration-record contamination

Signal datasets contain two known populations. Mean ARI for adaptive diffusion
+ NNLS is .877 on independent Bernoulli data and .649 on categorical data,
versus .763 and .400 for average linkage on the same datasets. Mean numbers of
clusters are 2.295 and 2.230, respectively. These recovery improvements coexist
with the null failures above.

On the 200 independent Bernoulli signal datasets, adaptive NNLS has mean
generator-null split FDP .0628 (descriptive dataset-bootstrap 95% interval
[.0380, .0902]). Generator-null truth is equality of the two known population
mixtures, computed from exact integer population counts in the children.
Only reached `split` decisions count as sibling rejections; `pass_through`
decisions are recorded separately.

Nonnull calibration records occur in 118/200 independent Bernoulli signal runs
and 194/200 categorical signal runs for adaptive NNLS. Their counts are
298/10,846 and 1,320/11,111 calibration records. These overlapping records are
not independent trials. Saved weights are base support weights, not the
focal-context kernel weights, so this identifies contamination without proving
its causal contribution to the fitted tail.

The high-dimensional Gaussian null has zero false partitions, but every one of
the 200 signal datasets also returns a single cluster for all three supported
methods: mean ARI zero. Their root test dimension is 39. Zero observed null
errors here cannot be presented as successful calibration **and** recovery.

Prespecified alpha sensitivity also fails to resolve the Bernoulli null issue:
adaptive NNLS gives 256/500, 316/500 and 436/500 false partitions at sibling
alpha .005, .01 and .05. Corresponding signal mean ARIs are .892, .877 and .777.
These evaluation outcomes were not used to tune or replace the current rule.

## Isolated adaptive-Wald calibration

Six million independent Gaussian contrasts were generated, including four
million held-out draws. Comparisons within a seed/dimension share their draws.
At alpha .01 and adaptive energy fraction .90:

| Spectrum | Current held-out rate, seed 2026092801 | Current held-out rate, seed 2026092802 |
|---|---:|---:|
| Fully tied `(1,1)` | .009949 | .009993 |
| Distinct `(2,1)` | .014464 | .014595 |
| Partially tied `(4,2,2,1)` | .015531 | .015342 |

Complete groups fix the fully tied example. Selecting rank from the tested
contrast across distinct groups still invalidates the ordinary fixed-rank
chi-square tail in these examples. All 16 held-out adaptive cells for distinct
or partially tied spectra show excess rejection after Holm correction over
216 cells, across alpha .001/.005/.01/.05. No held-out current fixed-rank control
is flagged.

All 2,688 vectorized-versus-kernel checks pass. Within repeated eigenspaces,
100,000-contrast rotation checks produce zero current-rule rank or rejection
changes; the largest p-value difference is 2.78e-15. The generic complete-pipeline
methods above disable adaptive energy rank, so this isolated defect does not
by itself explain their false-partition rates.

![Isolated adaptive-Wald calibration](/Users/berksakalli/Projects/kl-te-cluster/reports/eigenspace_reassessment_20260928/wald_rank_calibration.png)

See [Wald experiment details](/Users/berksakalli/Projects/kl-te-cluster/reports/eigenspace_reassessment_20260928/wald_rank/README.md)
for exact binomial intervals, the matched historical replay and compressed
per-contrast outcomes.

## Coverage, verification and remaining questions

The complete-pipeline evaluation contains **23,800 method/alpha attempts**:
19,595 successful outputs, 4,201 explicit skips and four unsupported outcomes;
zero execution errors and zero audit errors. The skips comprise 4,200 native
continuous/Hamming incompatibilities and one duplicate-geometry case. The four
unsupported rows lack internal empirical-null support: one iid-null dataset
under average linkage at three alphas and complete linkage at .01. Unavailable
outputs remain unavailable, including in uncertainty bounds.

The current 47-case baseline passed before report construction. All matched
benchmark input/status/assignment checks pass; all 216 isolated rejection-count
rows were independently decoded from saved outcomes. The independent pipeline
audit regenerated all 4,200 input hashes, checked 1,184,337 node records,
reconstructed all 19,595 completed-run aggregates and 4,140 summary/tail
calculations, and passed 39 deterministic production reruns covering successful,
skipped and unsupported outcomes. Every full-pipeline traversal replay matches
its dispatcher labels. Complete per-edge spectral ranks were not retained;
their aggregate fields were checked in the selected reruns. Source hashes remained fixed
throughout both runs. Timing is descriptive: the matched spectral panel took
915.9 seconds and the concurrent six-worker pipeline evaluation 481.8 seconds;
these are not controlled before/after speed comparisons.

The finite panel does not establish transfer to all dimensions, phylogenetic
dependence, real scRNA inputs, every tree method, every support policy or a
selected-Gaussian prototype. It does identify the next calculation problems:
the reference law after contrast-selected rank, the selected root/sibling tail
after the whole fitted pipeline, and the effect of contaminated/dependent
calibration records. Whole-tree control needs separate justification from
geometric invariance and signal recovery.

No production method, default or calibration rule was changed in this
reassessment. No commit or push was made. The prior repair and historical
reports were preserved.

## Reproducible artifacts

- [Frozen complete-pipeline protocol](/Users/berksakalli/Projects/kl-te-cluster/reports/eigenspace_reassessment_20260928/full_pipeline/protocol.md),
  [runner](/Users/berksakalli/Projects/kl-te-cluster/reports/eigenspace_reassessment_20260928/full_pipeline/run_reassessment.py),
  [manifest](/Users/berksakalli/Projects/kl-te-cluster/reports/eigenspace_reassessment_20260928/full_pipeline/evaluation/manifest.json).
- [Null partition summaries](/Users/berksakalli/Projects/kl-te-cluster/reports/eigenspace_reassessment_20260928/null_partition_summary.csv),
  [root tails](/Users/berksakalli/Projects/kl-te-cluster/reports/eigenspace_reassessment_20260928/full_pipeline/evaluation/root_tails.csv),
  [per-run outcomes](/Users/berksakalli/Projects/kl-te-cluster/reports/eigenspace_reassessment_20260928/full_pipeline/evaluation/runs.csv),
  [compressed node records](/Users/berksakalli/Projects/kl-te-cluster/reports/eigenspace_reassessment_20260928/full_pipeline/evaluation/nodes.csv.gz).
- [Signal recovery](/Users/berksakalli/Projects/kl-te-cluster/reports/eigenspace_reassessment_20260928/signal_recovery_summary.csv),
  [paired signal comparisons](/Users/berksakalli/Projects/kl-te-cluster/reports/eigenspace_reassessment_20260928/paired_signal_comparisons.csv),
  [alpha sensitivity](/Users/berksakalli/Projects/kl-te-cluster/reports/eigenspace_reassessment_20260928/alpha_sensitivity.csv).
- [Summary and plotting script](/Users/berksakalli/Projects/kl-te-cluster/reports/eigenspace_reassessment_20260928/summarize_and_plot.py),
  [null-rate PDF](/Users/berksakalli/Projects/kl-te-cluster/reports/eigenspace_reassessment_20260928/null_partition_rates.pdf),
  [Wald-rate PDF](/Users/berksakalli/Projects/kl-te-cluster/reports/eigenspace_reassessment_20260928/wald_rank_calibration.pdf).
- [Independent pipeline verification](/Users/berksakalli/Projects/kl-te-cluster/reports/eigenspace_reassessment_20260928/full_pipeline/evaluation/verification.json),
  [verification script](/Users/berksakalli/Projects/kl-te-cluster/reports/eigenspace_reassessment_20260928/full_pipeline/verify_pipeline_outputs.py).

Evaluation command, from the repository root (choose a fresh output directory
to preserve the recorded run):

```bash
TBS_N_JOBS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg \
  .venv/bin/python reports/eigenspace_reassessment_20260928/full_pipeline/run_reassessment.py \
  --phase evaluation --null-replicates 500 --signal-replicates 200 --workers 6 \
  --output-dir reports/eigenspace_reassessment_20260928/full_pipeline/evaluation
```
