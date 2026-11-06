# Diffusion topology, NNLS branch lengths and contrast-energy rank

Executed 2026-11-02 in the current working tree. This is a diagnostic comparison,
not a production-method change or a selective-calibration repair.

## Protocol and reproduction

```sh
uv run --no-sync python -m benchmarks.diagnostics.calibration.diffusion_factorial
uv run --no-sync python -m benchmarks.diagnostics.calibration.diffusion_factorial_trace
uv run --no-sync pytest -q tests/core/test_diffusion_factorial.py tests/core/test_adaptive_diffusion_runner.py tests/tree/test_optimized_branch_lengths.py
```

The complete paired 2×2×2 design uses 180 independent datasets and 1,440 pipeline
runs. Each dataset is evaluated in all eight conditions. Nodes and treatment
runs are not independent replicates. The 72-run pilot checked execution and
invariants; its overlapping seeds are not additional observations.

| Factor | Off | On |
| --- | --- | --- |
| Diffusion topology | Average linkage on canonical original-feature distance: standardized Euclidean for continuous; Hamming for binary/one-hot | Average linkage on existing feature-family adaptive diffusion distances |
| NNLS branch lengths | Every branch length zero; edge variance policy `none`; sibling branch factor exactly one | Fixed-topology NNLS against original-feature mean squared standardized Euclidean distance; normalized branch factors for both edge and sibling tests |
| Contrast-energy rank | Existing nominal spectral rank rule, minimum dimension 2 | Existing smallest PCA-prefix rank retaining 0.9 of observed contrast energy |

Diffusion uses existing family settings: continuous Euclidean k=15 and at most
50 components capped by raw feature count; Bernoulli/categorical Hamming k=10
and at most 30 components; time=3, bandwidth `-1/(d+2)`, epsilon=`median`.
This contrasts the complete existing topology constructors, not individual
kernel, bandwidth, normalization or diffusion-time choices.

NNLS uses all unordered pairs in this panel (cap=50,000), seed=0, tolerance
1e-5, max iterations=1,000, and rejects nonconverged fits. Rank comparisons
reuse identical topology and branch lengths. NNLS comparisons reuse identical
topology. Topology-only effects are the comparisons with NNLS off: NNLS-on
topology comparisons necessarily include the induced change in fitted lengths.
Assertions and stored hashes check these invariants for every dataset.

All tests use the original feature matrix. Alpha remains .001 for edges and
.01 for siblings; covariance is `guarded_within_child` with minimum child
size 8; sibling method is `projected_wald_inflation`; support-threshold
enforcement is false; traversal uses passthrough=true. Numerical inflation
estimates and gate support are recomputed downstream of each intervention:
the estimator policy is held constant, not its fitted values. The experiment
therefore estimates pipeline effects, not effects with covariance matrices or
fitted inflation factors artificially fixed.

Null panel: 50 datasets per family, 80 rows each; 24 independent standard
Gaussian coordinates, 24 independent Bernoulli(.5) coordinates, or 24
independent uniform three-category variables (72 one-hot columns). Seeds are
2026110203 + family_index×10000 + replicate, with family order continuous,
Bernoulli, categorical. These are explicit iid nulls, not fitted nulls.

Signal panel: five consecutive seeds starting at each registered seed for
`gauss_clear_medium_continuous`, `gauss_moderate_3c_continuous`,
`binary_low_noise_4c`, `binary_moderate_4c`, `cat_clear_3cat_4c`, and
`cat_clear_4cat_4c`. Existing generators and case sizes are unchanged.

## Results

All 1,440 conditions returned supported output. All 360 NNLS fits converged
(each fit is reused for the two rank conditions). No diffusion fallback was
used. Relative NNLS RMSE ranged from .03214 to .19024; successful optimization
does not establish that fitted lengths are stochastic variance parameters.

False final splits under the global iid null, each count out of 50:

| Topology | NNLS | Energy rank | Gaussian | Bernoulli | Categorical |
| --- | --- | --- | ---: | ---: | ---: |
| Original distance | Off | Off | 0 | 0 | 0 |
| Original distance | Off | On | 0 | 0 | 0 |
| Original distance | On | Off | 0 | 0 | 0 |
| Original distance | On | On | 0 | 2 | 0 |
| Diffusion | Off | Off | 29 | 0 | 0 |
| Diffusion | Off | On | 31 | 4 | 0 |
| Diffusion | On | Off | 30 | 4 | 2 |
| Diffusion | On | On | 29 | 18 | 0 |

The Gaussian topology-only difference is 29 additional false splits on the
same 50 datasets with branch time off and nominal rank. Its diffusion-arm
split fraction is .58, with exact two-sided binomial .95 interval
[.4321, .7181]. **Zero splits in the original-distance arm is not proof of
calibration**: its smaller root child has 1–2 leaves in all 50 datasets
(median 1), within the known small-node covariance regime. Diffusion's smaller
root child ranges from 7 to 40 (median 24.5); only one falls below 8. The root
edge gate opens in 0 original-distance versus 34 diffusion datasets for that
comparison. The topology intervention changes both selection and node sizes;
this experiment does not separate those two pathways.

On Bernoulli diffusion trees the four (NNLS, energy) false-split counts are
(0,0)=0, (0,1)=4, (1,0)=4, (1,1)=18. The dataset-averaged difference of
differences is (18−4−4+0)/50=.20, so adding the two separate count changes
does not reproduce the combined result. These are descriptive panel effects,
not independent binomial samples or a multiplicity-adjusted discovery claim.
The root edge gate actually opens less often with NNLS: 17→13 for nominal
rank and 41→36 for energy rank. More final splits cannot be explained by more
root edge openings alone.

Categorical null results also show no universal monotone rank effect: the
diffusion/NNLS nominal-rank arm splits 2 datasets, while energy rank splits 0.
Root gates open frequently even in conditions with zero final splits; final
partitions alone conceal individual gate behavior.

Signal recovery is insensitive in five of the six selected case families:
every condition has mean ARI=1 there. These relatively easy examples do not
establish equivalence on difficult signals. In `cat_clear_3cat_4c`, mean ARI
is .80668 for original distance / no NNLS / nominal rank; either NNLS or energy
rank raises it to .89487; diffusion reaches .97340 in all four NNLS/rank
conditions. Each mean uses five datasets. The topology-only paired ARI gain is
.16672. Therefore diffusion improves one measured recovery example while
showing a serious Gaussian-null failure in this panel.

## Why NNLS can increase significance

The trace script replays Bernoulli replicate 3, diffusion topology, energy rank.
It reconstructs raw sibling records and the empirical-null fit using the same
gate pipeline, then asserts agreement with the full diagnostic output.

| Root quantity | No branch time | NNLS branch time |
| --- | ---: | ---: |
| Raw sibling Wald statistic | 109.87885697 | 99.59535746 |
| Fitted inflation c_hat | 8.60026112 | 5.98048832 |
| Adjusted sibling p | .01242258 | .00225679 |
| Final clusters | 1 | 2 |

The branch factor is 1.103253. The raw statistic falls by this factor, as
expected, but the fitted inflation falls more. With unit reference scale,
the adjusted statistic is Q/c_hat; its increase explains the smaller p-value.
The raw-statistic ratio assertion passes at relative tolerance 1e-10; the
reconstructed adjusted p agrees at 1e-12. The trace supports this specific
mechanism; it does not establish that every NNLS-induced split has the same
cause. All raw records and the root calibration decision are retained.

## Evidence and limits

Outputs are in `raw/assets/diffusion-factorial-20261102/`:

- `manifest.json`: source/lock hashes, versions, seeds, configuration and Git head.
- `runs.csv`: all 1,440 outcomes, labels, root metrics and data/topology/branch hashes.
- `root_annotations.csv.gz`: root and both root-child annotations for each run.
- `summary.csv`: condition means, counts and exact binomial intervals.
- `paired_contrasts.csv`: within-dataset differences at each fixed pair of other
  factors; t intervals are exploratory, unadjusted and unreliable for rare
  binary events or five-replicate nonnormal outcomes.
- `nnls_reversal_trace.json`, `nnls_reversal_records.csv`: the traced reversal.

Fifty null datasets cannot certify .001 or .01 tails. Even 0/50 has an exact
two-sided .95 upper bound of .0711. The internal `ok` status means the program
had support under its current policy, not that its p-values are calibrated.
The working-tree empirical-inflation cap already present before this task is
included and identified by source hash. No production statistics, alpha,
method registry, dependency or preset was changed. The restricted analytic
selected-root chi/F derivations do not calibrate adaptive diffusion selection.

The four new parity tests match the existing runner's root statistics,
p-values, degrees of freedom and ARI for both topologies and rank modes with
NNLS on. The full experiment checks unchanged topology and lengths within
the relevant comparisons. Together with the 18 existing diffusion/NNLS tests,
these establish execution fidelity, not the validity of the statistical null.
