# The ceiling-3 rule on the CHOIR cell lines and the continuous benchmark

Run on 2026-11-09 (user choice: test the ceiling-3 rule of the single-population split bound on the real cell
lines and the continuous benchmark). Local folder: `results/analyses/ceiling_rule_real_20261109` (git-ignored);
protocol, deviations, code and outputs in `ceiling-rule-real-and-benchmark-scripts-20261109.py.txt`.

## Setup

- Rule (unchanged from `single-population-split-bound-20261108.md`): for each tested sibling split, k invariant
  directions of its subtree, rho = pi (1 - pi) D^2 along them, p = min{1, sum_{i<=k} binom(n-1, i)
  P(chi2_k(3 n) >= n rho)}, p = 1 when k = 0; then the usual sibling BH step.
- Both parts use the production opt-in path (`sibling_split_count_rule_out=True`). Arm `rule_out` runs it as is;
  arm `ceiling3` swaps `apply_split_count_rule_out` for the ceiling-3 p-value on the same leaf data, leaf rows
  and statistic. No production change.
- Real: 48,879 cells, 30 PCs, saved diffusion coordinates, registry parameters, edge alpha 1e-3, sibling alpha
  1e-2, both registered methods. Benchmark: the 9 continuous cases, both methods, default alphas.
- First run failed: the label `ceiling3` is rejected by the sibling FDR step; relabelled
  `split_count_rule_out:ceiling3` and rerun from the start (deviations.md).

## Replay checks

- Real: `rule_out` through the production path gives the 2026-11-08 harness result exactly (69 and 664
  clusters, ARI 0.0182 and 0.7652).
- Benchmark: today and `rule_out` rows equal the 2026-11-08 run in found clusters and ARI, 36/36.

## Real cell lines (190 lines)

| Method, arm | Clusters | ARI | NMI | Lines 90/90 | ARI vs CHOIR_1 | Clusters >= 90% pure | Clusters < 20 cells |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Markov/Ward, today (2026-11-08) | 954 | 0.549 | 0.863 | 41 | 0.445 | 761 | 374 |
| Markov/Ward, rule_out | 664 | 0.765 | 0.901 | 80 | 0.514 | 494 | 253 |
| Markov/Ward, ceiling3 | 217 | 0.833 | 0.940 | 119 | 0.538 | 164 | 25 |
| Default, today (2026-11-08) | 75 | 0.028 | 0.548 | 27 | 0.014 | 48 | 20 |
| Default, rule_out | 69 | 0.018 | 0.497 | 29 | 0.009 | 49 | 12 |
| Default, ceiling3 | 64 | 0.028 | 0.549 | 29 | 0.014 | 46 | 9 |

Why lines miss 90/90 (added after the endpoints): Markov/Ward split only 72 to 16, merged with other lines only
11 to 33, both 27 to 22; default tree unchanged (145 merged only; largest cluster 28,611 cells under ceiling3,
32,680 under rule_out).

Tested splits kept before BH (Markov/Ward, rule_out to ceiling3): 32-159 cells 17.4% to 1.6%, 160-999 14.3% to
12.9%, 1,000+ 34.0% to 40.8%; at most 31 cells none in either. Of the same 1,838 splits, rule_out keeps 258 and
ceiling3 121; 81 are kept by both, 177 only by rule_out, 40 only by ceiling3. Default tree: 22 vs 14 kept, 12 by
both. Median rho of Markov/Ward splits with invariant directions: 0.97 (32-159 cells), 1.00 (160-999), 3.13
(1,000+).

## Continuous benchmark

`ceiling3` gives the same found clusters and ARI as `rule_out` on all 9 cases and both trees (mean ARI 0.695,
today 0.880 default and 0.852 Markov/Ward; exact cluster count 4 of 9, today 6 and 5). As predicted, the losses
are the subtrees with n <= p + 1 cells (60 cells in 40 features, 45 in 30, 40 in 20,000), which have no
invariant directions under either rule; in the 180-cell, 136-feature case 2 of 13 tested splits have
directions and their rho is at most 1.62. The null case stays at one cluster (largest rho 0.28).

## Reading

- On the real data the ceiling-3 rule is stricter than separability for small subtrees and looser for large
  ones. Markov/Ward moves toward the 190 lines on every measure: 664 to 217 clusters, ARI 0.765 to 0.833, lines
  recovered 80 to 119, small clusters 253 to 25. The cost is more merged lines (11 to 33 merged only).
- The default tree does not change materially: its giant cluster comes from the topology, as before.
- On the continuous benchmark the ceiling-3 rule changes nothing; the losses there are set by having no
  invariant directions in small, wide subtrees, which neither rule addresses.
- Not measured: error control on real data (no known no-group parts); populations with long flat tails, where
  3 is not a ceiling. No production change.

## Merge check (2026-11-09, later; `merge_protocol.md`, `merges.py`)

Markov/Ward assignments only, no new fits. Lines merged with other lines under ceiling3 (merged only or both):
55 of 190 (29%).

- Size: merged lines are not small overall (median 206 merged only, 238 both, 226 recovered), but the very
  smallest are over-represented: 5 of the 10 smallest lines, 6 of the 19 smallest, 16 of the 38 smallest.
- 56 ceiling3 clusters hold two or more lines at >= 5% of their cells. For the two largest lines of each, the
  true-label split was tested with the same rule:
  - 17 clusters have at most 49 cells and no invariant directions; neither rule can split them.
  - 27 of 39 true-label splits pass ceiling 3 easily (rho 5.5 to 27.7, median smaller share 0.32). There the
    tree did not offer the line boundary as a split, or the gates stopped before it; the ceiling is not the
    cause. Rule-out separated 16 of these 27 pairs.
  - 12 fail ceiling 3 (rho 1.5 to 5.8, median smaller share 0.13); 9 of them pass the Gaussian ceiling, where
    g(pi) is 0.3 to 0.7 for these small shares. These are the merges the constant ceiling causes: small groups
    beside a large one. Rule-out separated 7 of the 12.
  - The published CHOIR_1 clustering separates every one of the 56 pairs.
- Expectation of the merge protocol (merged lines are small; their true splits have small pi and rho between
  g(pi) and 3) held only for the 12 failing pairs; most merges come from the tree, not the ceiling.

## Why the tree hides the boundary (2026-11-09, later; `trace_protocol.md`, `trace_fit.py`, `trace.py`)

Markov/Ward refit with ceiling3, tree and annotations kept. Replay: the same 1,838 tested splits and 120 BH
acceptances as the earlier ceiling3 run (p-values equal up to float rounding). For the 27 merged pairs whose
true-label split passes ceiling 3, the traversal stopped at the cluster's own top node in every case:

- 8 pairs: the top split is the line boundary (>= 85% of each line in its own child), but its rho is 1.6 to 4.5,
  below what a split of that size and dimension needs for p <= 0.01 (4.3 to 7.4; the true-label rho is 5.8 to
  27.7). Cells put on the wrong side, and third lines in the cluster, pull rho down; the count of flat cuts in
  k = 6 to 12 invariant directions sets the needed rho well above 3.
- 2 pairs: the top node's children fail the edge gate; no sibling test is run.
- 17 pairs: the top split cuts elsewhere (both lines on one side or each line split), median rho 1.22, max
  5.18; it is rejected and the top-down traversal never reaches deeper nodes. A deeper node isolates both
  lines (F1 >= 0.9) in 4 of them; in the other 13 no node of the tree isolates both lines.
- Every top split's rho is below its needed rho; 8 are above 3. So the binding cost is the margin from the
  flat-cut count, not the ceiling itself, together with a Ward tree on diffusion coordinates whose splits
  inside these clusters do not follow the line boundary in the 30-PC test space.
