# Fixed-partition Gaussian attribution and separate projection/inflation restoration

Executed 2026-11-02 with three Sol subagents: baseline implementation,
restoration implementation, and independent mathematical/code/result review.
Wolfram evaluator inputs and outputs are preserved in
`raw/inbox/diffusion-attribution-wolfram-20261102.md`.

## Design

Reuse the 50 diffusion trees selected from iid Gaussian N(0,I24) matrices with
80 rows in the preceding factorial. Selected-data and topology hashes match
that experiment exactly. Freeze the entire tree and leaf labels, not just root
child counts. Repopulate each tree with either its original selected matrix
or an independent iid N(0,I24) matrix. Independent seeds are 2026110303+r;
selected seeds are 2026110203+r, r=0,...,49. One pair is one independent
replicate; node records and treatment arms are dependent observations.

The baseline crosses these data modes with known identity versus the actual
production `guarded_within_child` covariance. Use all 24 features, identity
projection, zero lengths, inflation factor 1 and no contrast-energy rank
selection. The primary result is the unadjusted root sibling test at alpha=.01.
Recorded full-row child-parent algebra checks at alpha=.001 are explicitly
different from production's split-row edge annotations. No final traversal or
FDR claim is made from these unadjusted root results.

For each fixed binary root partition, write s=1/nL+1/nR and
z=(meanL-meanR)/sqrt(s). Under independent testing rows and known identity
covariance, Q=z'z has chi-square 24 law conditional on the entire selected tree.
This argument does not apply when the same observations selected the groups.

## Baseline findings

| Testing observations | Covariance | Root rejections at .01 | Mean Q |
| --- | --- | ---: | ---: |
| Selected | Known identity | 50/50 | 98.9434 |
| Independent | Known identity | 0/50 | 24.2202 |
| Selected | Production guarded | 50/50 | 169.5356 |
| Independent | Production guarded | 8/50 | 30.4840 |

The known-covariance comparison isolates selection as sufficient for the
observed failure: PCA, covariance estimation, NNLS and fitted inflation are
absent. The independent-data covariance comparison identifies a separate
reference-law mismatch. It is not evidence of independence failure.

Live covariance representation matters, not just the policy name:

| Actual root covariance | Number of pairs | Independent rejections | Independent mean Q |
| --- | ---: | ---: | ---: |
| Dense pooled within children | 25 | 7/25 | 35.6872 |
| Diagonal pooled within children | 24 | 1/24 | 24.9085 |
| Dense parent-total fallback | 1 | 0/1 | 34.2147 |

At p=24, children with fewer than 25 leaves store diagonal variances. Thus the
dense pooled regime requires both children at least 25, not merely passing the
minimum-child-size 8 guard. A child below 8 triggers parent-total fallback in
this root panel. These regimes are recorded for every row.

Under a fixed independent partition, ordinary dense unregularized pooled
covariance with nu=n-2=78 gives

\[
\frac{55}{78\cdot 24}Q\sim F_{24,55}.
\]

Wolfram computes rejection probability .2342701084 if this Q is instead
compared to the nominal .01 chi-square 24 threshold 42.97982014. The ideal mean
is 35.32075, close to the observed 35.6872 in the 25 independent dense cases.
Using this F diagnostic gives 0/25 independent rejections but still 25/25
selected rejections. It is not a selected-data repair.

The independent diagonal-pooled reference is a sum of 24 independent F(1,78)
variables. Dense parent-total fallback instead gives Q/79~Beta(12,27.5).
Wolfram checks the total/pooled rank-one identity and Beta transformation.
Neither is the dense pooled F reference. Production adds numerical ridge;
the largest measured dense Q difference from the unregularized calculation
was 3.22e-8, so the ideal formulas are not claimed as exact runtime laws.

## Restoration design and findings

For each testing matrix, obtain the native production split-row edge
annotations once at alpha=.001 with zero lengths and no contrast-energy rank
selection. Freeze all resulting null-like, blocked, support-weight and
dependency-group fields across covariance, projection and inflation arms.
Edge metadata can differ between selected and independent matrices; within
each matrix it is held fixed. No tree is rebuilt on independent observations.

Cross full-feature versus exact native production projection with inflation 1
versus the existing empirical-null fit, for both covariance modes and both
data modes. The same native all-leaf parent-total spectral basis/rank is used
across covariance modes; identity covariance changes the contrast covariance
only. No identity-chart refit or forced reduced rank is substituted. The
native MP rank and tie expansion remain as in production, while adaptive
contrast-energy rank remains off.

All root projections span 24 dimensions. Therefore restoring native projection
has no direct raw-root-statistic effect in these data. Nonroot ranks differ,
so restoring projection changes the records and context used by the empirical
inflation model. The role assignments remain fixed, but fitted values and
dimension-dependent support summaries are allowed to respond.

Counts of descriptive root p<.01 across all 50 datasets:

| Data and covariance | Full feature, c=1 | Native projection, c=1 | Full feature, empirical c | Native projection, empirical c |
| --- | ---: | ---: | ---: | ---: |
| Selected, known identity | 50 | 50 | 34 | 13 |
| Selected, guarded | 50 | 50 | 34 | 28 |
| Independent, known identity | 0 | 0 | 0 | 0 |
| Independent, guarded | 8 | 8 | 7 | 7 |

All root empirical predictions have available internal support under the
current non-enforcing support policy. Availability is not calibration.
Only 34 selected roots and zero independent roots are focal/eligible according
to the frozen native edge roles. Predictions on skipped roots are descriptive
and explicitly withheld from `production_gate_p`; in particular, 7/50 fresh
adjusted predictions are not seven executed production gate rejections.
Among the 34 eligible selected roots, empirical rejection counts are 34 for
full features,13 for native projection with identity contrast covariance,
and 28 for native projection with guarded covariance.

Native guarded/projected/empirical eligible-root p-values match the original
factorial exactly, including 34 eligible roots and 28 root rejections. That
factorial had 29 fragmented datasets: replicate 38 has root p=.0816779 and an
open edge path, with three final clusters from downstream traversal. Root
rejection and final fragmentation are different measured outcomes.

## Implications for the analytic repair

1. A known-covariance selected-region law is already necessary for diffusion
   selection. Improving covariance estimation alone cannot remove 50/50
   selected identity-covariance rejections.
2. Estimated-covariance references must distinguish dense pooled, diagonal
   pooled and parent-total regimes. The standard chi-square reference is
   insufficient even without tree selection in the dense finite-sample case.
3. Production projection is a root-span no-op here, but it changes empirical
   inflation through lower nodes. An analytic repair must account for the
   selected calibration population and fitted scale; fitting a scale does not
   establish a conditional selected-null law.

No diffusion selection region has been derived here. The earlier quadratic
region solver for average squared-Euclidean linkage does not apply unchanged
to data-dependent diffusion construction. No thresholds or production methods
were modified.

## Reproduction and artifacts

```sh
uv run --no-sync python -m benchmarks.diagnostics.calibration.diffusion_attribution --replicates 50
uv run --no-sync python -m benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_attribution_restoration --replicates 50 --output-dir raw/assets/diffusion-attribution-20261102
uv run --no-sync pytest -q tests/core/test_diffusion_attribution.py tests/core/test_diffusion_attribution_restoration.py tests/core/test_diffusion_factorial.py tests/core/test_adaptive_diffusion_runner.py tests/tree/test_optimized_branch_lengths.py
```

The output directory contains 200 baseline root rows,800 restoration root
rows and 63,200 compressed all-node restoration records; summaries, regime
counts, contracts and provenance accompany them. The restoration table
includes the baseline arms, so these are not 1,000 independent experiments.
Independent review verified every baseline/restoration root row, exact
agreement of baseline arms, frozen roles/bases and native factorial parity.

Fifty pairs provide attribution evidence, not precise tail validation. A0/50
count has exact two-sided .95 upper bound .0711. There is no general claim
about other dimensions, covariance models, diffusion settings or difficult
signal recovery. Full traversal/BH calibration remains untested here.
