# Root sampler validation, cut geometry and clustering performance

2026-11-09, Europe/Vienna. Three Sol agents completed the mathematical,
calibration and performance analyses. The coordinator added a fresh count-depth
panel and independently checked the combined evidence. Production code, defaults,
dependencies and the Git index were not changed.

## Finding

The evidence supports the maintained root sampler within its stated model and
shows substantial recovery gains from B on the CHOIR Markov/Ward tree. It does
not establish general calibration for A, B, C or the final partition. The gains
depend on separation, group balance, dimension and topology. B and C can improve
ARI while losing a small group, and a stable one-cluster answer can be wrong for
recovery. A scalar correction alone cannot generally repair the reference shape
changed by selection. A Gaussian conditional reference can preserve the sibling
score in principle; its implementation for the complete TBS procedure is open.

The original Claude conversation is session
`37455824-398b-44be-86b4-f026b44f5614`, titled *Root sampler validation framework*.
Its five stages and later user constraints are preserved in
[selected source excerpts](root-validation-claude-context-20261109.md).

## Which question each result answers

| Target | Evidence in this study | Boundary |
| --- | --- | --- |
| Equality of the true means averaged over two selected groups | Known-sigma root/ancestor-path chi chain | Independent Gaussian rows with known common spherical noise; native root and immediate children |
| One common Gaussian population at a fixed root | Exact conditional Beta/all-label derivation and numerical controls | Same Gram-measurable subspace for every candidate; fixed union; global rather than selected-pair inference |
| Selected-pair test with unknown common covariance | Chen–Gao polar/Haar literature and the rank argument | Stronger common-mean Gaussian null; TBS adaptation derived, not implemented |
| One expression profile under varying sequencing depth | Refitted count simulations | Biological robustness target; it does not satisfy the root chain's Gaussian contract |
| Recovery of cell-line or simulated group labels | Partition metrics, rare-group matching and perturbation stability | Agreement with labels is not a statistical Type I error guarantee |

The known-noise contrast null permits different row means within either group
when the two average means coincide. This is broader than the common-mean null
needed for the unknown-covariance row-orientation argument. A model class must
define “one population”; an unrestricted iid law can itself have a mixture
representation. Neither mean equality nor rejection identifies a component count.

## Root and immediate-child validation

All 24 historical 999-step root chains reproduced exactly; seven sampled
immediate-child records matched. The new main panel used 119 steps, fixed seeds
and fixed proposal scale, followed by a separately fixed 999-step paired subset.

| Known-noise null | Datasets | Root rejections at .05 | Family false rejections at .05 | One-sided 95% family upper limit |
| --- | ---: | ---: | ---: | ---: |
| n=9, p=3 | 1200 | 38 | 19 | .02315 |
| n=12, p=8 | 120 | 6 | 3 | .06334 |
| n=24, p=5 | 120 | 6 | 7 | .10677 |

At n=9,p=3 the root one-sided upper limit is .04130. These primary endpoints
are conservative at this fixed budget. The larger/wider strata remain too
imprecise to substantiate the same nominal claim. The complete main run used
1,664 datasets including signals and fitted-scale sensitivity: 4,441 moved
chains, five stalled chains and one unsupported observed numerical tie. All
attempts remain in the relevant denominators.

On 20 preselected balanced signal inputs, raising the budget from 119 to 999
steps increased root rejection at .05 from 6 to 12, and any true family rejection
from 2 to 8. Median family runtime increased from 1.42 to 11.65 seconds. The
30 paired null inputs are a sensitivity check, not a precise rate estimate.
In the strong four-group panel the tree recovered the intended top split on
32/32 datasets, but selected child cuts often mixed the intended subgroups;
immediate-child rejection counts were 2/28 and 2/32. Budget and topology both
limit detection. Movement is not a certificate of the selection event or model.

[Calibration report and protocols](root-validation-calibration-20261109.md)
record all counts, intervals, truth definitions and unresolved cases.

## Mathematical findings and useful generalisations

Four new Wolfram evaluations independently checked the counterexample, uniform
boundary moments, R transformation, Beta/F relationship, dimensional obstruction,
tail expansion and group-imbalance thresholds. The preceding A/B/C audit's
limitations remain: ceiling 3 needs a mode-cut restriction for the unimodal
argument; B's uniform-boundary leading variance is 12 rather than the correct
96/5; C lacks its cited C2 regularity and has no stated uniform-marginal model
above dimension 3. A also needs to account for covariance and rank selection.

A new exact root construction is available under a restricted Gaussian null.
Let Y be the centered n-row data, choose the same rank-k feature subspace from
G=Y'Y for every split, and let u be a normalized fixed label contrast. Then

    R² = u' Proj(col(YV)) u,
    R² | G ~ Beta(k/2, (n-1-k)/2),    0 < k < n-1.

This holds for arbitrary unknown common feature covariance because the
conditional row orientation is Haar. An exact global-root bound is

    p_all = min(1, (2^(n-1)-1) * BetaSF(R²)).

The family consists of all deterministic nonempty unordered label splits.
Replacing it with the sample's realized flat-cut count is not justified by
this argument. For fixed k the all-label critical value tends to R²=3/4,
equivalently rho=3. This is a conservative Gaussian threshold limit, not the
true Gaussian best-cut limit or a universal unimodality bound. Selected
descendant unions require additional conditioning.

The spectrum-only-rank control rejected 120/12,000 at .01; substituting a rank
chosen from observed contrast energy gave 256/12,000, with exact interval
[.01882,.02408]. Exhaustive maxima over all 127 nonempty unordered splits of
eight points gave 131/12,000 with the all-label correction. Full row-span
saturation at k=n-1 gives R²=1 for every split; dimension is therefore a
calibration and power issue as well as an implementation limit.

For two true Gaussian components in fixed retained coordinates, let pi be a
component share and Delta their Mahalanobis mean distance in the common within
covariance. The true-label ratio is exactly `rho=pi(1-pi)Delta²`. Exceeding 3
requires Delta above 3.464, 4.330 and 7.947 at pi=.5,.2,.05. This explains a
structural imbalance penalty in the population criterion; it is not by itself
a causal explanation of every four-group tree result.

The October 2026 [Chen–Gao preprint](https://arxiv.org/html/2610.17458v1)
provides a related selected-pair construction under a common-mean Gaussian null
with unknown shared covariance. Its reversible random-index rank argument
permits a deterministic recomputed sibling score. Our extension to the full TBS
score is a mathematical inference; no TBS implementation of that rotation
reference was made. Global root replay and conditioning on the observed pair
are distinct options and answer different hypotheses.

[Full derivations](root-validation-mathematics-20261109.md) and
[actual Wolfram requests/results](root-validation-wolfram-20261109.json)
preserve assumptions and source references.

## Partition quality beyond ARI

All eight saved CHOIR partitions were checked against the same 48,879 cell IDs
and 190 cell-line labels. Their old ARI, NMI, cluster counts and 90/90 recovery
reproduced. Macro F1 uses one-to-one maximum-total-F1 matching; unmatched lines
receive zero, and each line has equal weight. Rare F1 concerns the 48 smallest
lines, with lexical tie breaking.

| Markov/Ward test | Clusters | ARI | AMI | Macro F1 | Rare F1 | Lines recovered 90/90 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| A: production | 954 | .5488 | .8313 | .6173 | .7007 | 41 |
| B: ceiling 3/count | 217 | .8329 | .9346 | .7996 | .7148 | 119 |
| C: best-cut EC | 173 | .8152 | .9330 | .7463 | .5565 | 111 |
| Production rule-out | 664 | .7652 | .8839 | .7408 | .8037 | 80 |

B leads overall recovery here, while rule-out has higher mean F1 among the
smallest lines. C's cluster count is closer to 190 but its recovery is worse
than B's. B reduces same-line pair splitting from .572 to .085; C has more
different-line merging within predicted clusters (.270 versus B's .233).
The default symmetric/average A/B/C partitions retain ARI near .028 and
macro F1 near .18–.19; changing the sibling tail does not recover the boundaries
missing from that tree.

The fresh paired Gaussian-mixture study completed 448 full fits and 224
same-cell jitter comparisons across seven regimes, two trees, four arms and
four seeds. On imbalanced Markov/Ward inputs, B/C raise mean ARI from .332
to .616 while lowering macro F1 from .717 to .510; both leave the 12-cell group
unmatched in every seed. At lower separation they return one cluster on every
original and perturbed input: truth ARI 0, jitter ARI 1. In the larger balanced
sample they recover four clusters, whereas A averages 18.75. These differing
outcomes establish a trade-off, not one universally best construction.

Saved node relationships connect acceptance to size, rank, balance, depth,
reference weight and subspace-distance ties. For 167 matched Ward nodes with
spectral/distance records, test-space tie share is smaller among raw passes.
The nodes are selected, nested and dependent; these associations are descriptive.
Spectral multiplicity identifies a subspace, while individual eigenvectors
within a tied block can rotate. Crossing separations cannot both appear as
branches of one hierarchy. Neither a tie count nor multiplicity alone certifies
another population. [Full metrics and relationship tables](root-validation-metrics-20261109.md)

## Full count preprocessing and cut-level stress tests

The coordinator reran 100 fresh count datasets through normalization, log1p,
Seurat HVG selection, scaling, PCA, both diffusion geometries, NNLS and both
gates, giving 800 matched fits with no failed attempt. Six original count/tree
replays matched before these runs. All 800 outputs were independently rescored;
200 dataset/tree groups have identical input, PCA, distance and linkage hashes
across arms, and all 1,298 B/C split records match on their shared geometry.

| One-profile count null: false partitions out of 30 | Default A | Ward A | Default/Ward B | Default/Ward C | Default/Ward rule-out |
| --- | ---: | ---: | ---: | ---: | ---: |
| Log-depth SD .25 | 0 | 30 | 0 / 0 | 0 / 0 | 0 / 0 |
| Log-depth SD .42 | 4 | 30 | 0 / 0 | 0 / 0 | 0 / 5 |
| Log-depth SD 1.0 | 25 | 30 | 0 / 0 | 0 / 0 | 0 / 0 |

Zero of 30 has a one-sided exact 95% upper limit .09503. These results cannot
substantiate .01 error control. Five balanced and five imbalanced strong-profile
datasets were also run: B/C/rule-out recovered all four groups exactly on both
trees. Default A did too; Ward A averaged 20.4 and 26.6 clusters, with ARI .292
and .074. Their large signal differs from the weak/imbalanced Gaussian-mixture
panel above.

Separately, 10,000 datasets supplied 20,000 paired fixed-label and searched-cut
diagnostic records. On Gaussian n=32,k=20, B rejected 818/1,000 searched
coordinate-threshold cuts at .01. C returns 1 for every attainable statistic
in that configuration because its envelope exceeds one throughout the range.
That is structural loss of power. With shared-composition counts at two
randomly assigned depth levels, the searched cut rejects 977/1,000 for B and
143/1,000 for C after simplified fitted preprocessing. These are isolated
full-coordinate tail checks, not the production traversal or its estimated
invariant rank. Their generators also differ from the lognormal-depth count
panel. Agreement on the latter does not remove these reference-law failures.

## Reproducibility and remaining limits

The run completed 1,248 fresh end-to-end fits, 1,664 fresh chain families plus
50 paired longer-budget runs, historical replays, 20,000 cut records and the
Gaussian geometry controls. Frozen protocols, manifests, chains, assignments,
plots and scripts are in `results/analyses/root_validation_followup_20261109/`.
The main image is `comparison.png`, with PDF and SVG exports. Detailed execution
and archive fingerprints are recorded in `verification.json` and the evidence
manifest next to this report.

The prior audit's 119 recorded source/input hashes matched before this study.
The coordinator recorded 315 source/configuration hashes. Focused tests passed:
40 coordinator baseline tests, the agent's 45 sampler-focused tests, 13
mathematical-arm focused tests and the metrics integration check; these overlap
and should not be added as unique-test counts. The wiki lint passed on 262 pages
after synthesis. The original shorter-budget analysis source was recovered byte for
byte from its manifest after the runner gained a fixed longer-budget option.

Count execution used two processes and one-thread BLAS/OpenMP/Numba settings.
The existing neighbor search retains explicit `n_jobs=-1`, so its actual internal
thread count was not certified by the requested one-job backend context. All
replays and paired geometries matched. The unused exception-handling limitations
are recorded in the run's `counts/deviations.md`; all 800 attempts succeeded.

The CHOIR partitions were reused and rescored; full-cell preprocessing and
clustering were not repeated in this follow-up. The historical 935-dataset panel
was not rerun. No rotation-based TBS test was implemented. Covariance estimation,
count preprocessing, arbitrary descendants, adaptive traversal/pass-through and
final-partition guarantees remain outside the validated chain scope. The
experimental B/C references retain mathematical failures, despite recovery gains.

Evidence: [combined numerical summaries](root-validation-followup-20261109-evidence.json),
[coordinator/math/count sources](root-validation-followup-20261109-sources.py.txt),
[calibration sources](root-validation-calibration-20261109-sources.py.txt), and
[metric sources](root-validation-metrics-20261109-scripts.py.txt).
