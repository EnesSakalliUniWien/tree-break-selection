# What decides the edge test's projection dimension k on the cancer_pc cases

Run on 2026-11-12, Europe/Vienna. Local folder: `results/analyses/cancer_pc_edge_k_20261112` (git-ignored);
protocol, deviations, both scripts and the complete summary output are in
`cancer-pc-edge-k-scripts-20261112.py.txt`. Follow-up of `cancer-pc-edge-test-20261112.md`. Working tree on top
of `bc044e07`; registered cell of `tbs_family_diffusion` with the operator and linkage overrides, today's test,
`TBS_N_JOBS=1`. No production code was changed. Environment: Python 3.12.11, NumPy 2.3.4, SciPy 1.16.3,
uv 0.9.5.

## Question

At nodes of one true group the edge test's k has median 28 with operator "markov" and linkage "ward" and 6, 9
and 3 in the other three cells. What decides k?

## The rule, as read from the code

- The eigenvalues are those of the node's basis rows (one half of each child's rows) after whitening.
- Whitening covariance: the pooled within-child covariance of the basis rows themselves when both basis
  children have at least 8 rows and basis rows minus 2 is at least the number of features (30) ("within
  route"); otherwise the node's total covariance over all its rows ("total route").
- k = max(Marchenko-Pastur signal count, minimum 2), then extended through the numerical eigenvalue group at
  its boundary (`extend_rank_through_ties`, relative tolerance 1e-10).

## Setup

The 16 fits of the 4 cells x 4 cases. The per-node spectral worker and the Marchenko-Pastur estimator were
wrapped at run time while the edge test's basis is computed; 3,241 nodes recorded with all eigenvalues, the
signal count, k, the route, basis rows per child and test rows. Self-checks passed: found clusters and ARI as
before; k equals the k of the earlier run at 146 of 146 nodes; k recomputed from the recorded eigenvalues with
the two rules equals k at 3,241 of 3,241 nodes. The protocol with three predictions derived from the code was
written before the run; the k values of the 146 nodes were already known then.

## Results

Nodes of one true group with at least 40 cells (146):

| Route | Nodes | Signal count | k |
| --- | ---: | --- | --- |
| within | 51 | 0 at 18, 1 at 33 | 30 at 46, 29 at 4, 28 at 1 |
| total, fewer than 30 test rows | 46 | 0 | 30 minus test rows (3 to 10) at 40, 2 at 6 |
| total, 30 or more test rows | 49 | 0 | 2 at 49 |

- P1 supported: on the within route 29 eigenvalues below the first are numerically tied at 46 of 51 nodes (28
  or 27 at the other 5), so the minimum 2 is extended to 30 (29, 28). Example: eigenvalues 149.46, then 0.9375
  twenty-eight times.
- P2 supported: on the total route with fewer than 30 test rows the tied group starts at the first eigenvalue
  and has 30 minus test rows members at 43 of 46 nodes, and k = max(2, 30 minus test rows) at 46 of 46; with
  30 or more test rows there is no tie at the boundary and k = 2 at 49 of 49.
- All nodes with at least 40 cells (383): within route 180 (k = 30 at 168, 28 or 29 at 12); total route 203
  (k = max(2, 30 minus test rows) at 100 of 100 with fewer than 30 test rows; k = 2 at 103 of 103 otherwise).
- Why a node is on the total route (95 of the 146): the smaller basis child has fewer than 8 rows at 69, basis
  rows minus 2 is under 30 at 66.

Per cell, nodes of one true group with at least 40 cells:

| | symmetric/average | symmetric/ward | markov/average | markov/ward |
| --- | ---: | ---: | ---: | ---: |
| Nodes | 48 | 29 | 45 | 24 |
| Within route | 14 (29%) | 12 (41%) | 12 (27%) | 13 (54%) |
| k = 2 | 20 | 8 | 21 | 6 |
| k from 3 to 10 | 14 | 9 | 12 | 5 |
| k from 28 to 30 | 14 | 12 | 12 | 13 |
| k, median | 6 | 9 | 3 | 28 |

- P3 supported: markov/ward has the largest share of nodes on the within route (54%), so its median k is 28.
  The number of nodes with k from 28 to 30 is nearly the same in all four cells (14, 12, 12, 13).

Rejections by route, without branch time (edges of the 146 nodes, from the earlier `edges.jsonl`):

| | Edges | Tree's split, raw p <= 0.001 | Random splits, raw p <= 0.001 | Random splits, raw p <= 0.01 | Statistic per k, tree's split, median |
| --- | ---: | ---: | ---: | ---: | ---: |
| within route, k 28 to 30 | 102 | 84 (82.4%) | 0.66% | 2.89% | 2.74 |
| total route, k = 30 minus test rows | 80 | 2 (2.5%) | 0.01% | 0.06% | 1.16 |
| total route, k = 2 | 110 | 38 (34.5%) | 0.07% | 0.76% | 1.59 |

The within route per cell:

| | symmetric/average | symmetric/ward | markov/average | markov/ward |
| --- | ---: | ---: | ---: | ---: |
| Edges on the within route | 28 | 24 | 24 | 26 |
| Tree's split, raw p <= 0.001 | 16 | 18 | 24 | 26 |
| Tested by Tree-BH | 0 | 4 | 0 | 26 |
| Significant | 0 | 0 | 0 | 26 |

## Reading

- k is not estimated from signal at these nodes; the signal count is 0 or 1 everywhere. k is set by the
  whitening route and by the tie extension: 30 where the basis rows are whitened by their own within-child
  covariance, 30 minus test rows where the node has fewer than 30 test rows, 2 otherwise.
- On the within route the edge test is a test in all 30 directions. There the tree's split is rejected at
  82% of the edges, in every cell (16 of 28 to 26 of 26). With markov/ward these edges are reached by
  Tree-BH and all 26 are significant; in the other cells 0 to 4 of them are reached.
- With k = 30 minus test rows the projection consists of the directions the test rows do not span; random
  splits reject at 0.06% at the 0.01 level and the tree's split at 2.5% at 0.001. The test has almost no
  power and no excess there.
- The docstring of the estimator says a pure-noise node runs a low-rank test. With the tie extension that
  holds only on the total route with 30 or more test rows (k = 2).

## Deviations

One check line of the summary was restricted to nodes with at least 40 cells, and the per-cell table of the
within route was added after the first output. `measure.py` was not changed after its run.

## Not checked

Whether another k rule would remove the rejections of the tree's split (with k = 2 they are still 34.5%). Why
29 or 28 instead of 30 at 5 nodes with 32 to 34 basis rows. Cases with more features than 30, with correlated
noise or with real data. One seed per case.
