# Why the edge test is significant inside one Gaussian group on the cancer_pc cases

Run on 2026-11-12, Europe/Vienna. Local folder: `results/analyses/cancer_pc_edge_test_20261112` (git-ignored);
protocol, deviations, both scripts and the complete summary output are in
`cancer-pc-edge-test-scripts-20261112.py.txt`. Follow-up of `cancer-pc-operator-linkage-20261112.md`. Working
tree on top of `bc044e07`; registered cell of `tbs_family_diffusion` with the operator and linkage overrides,
today's test. No production code was changed. Environment: Python 3.12.11, NumPy 2.3.4, SciPy 1.16.3, uv 0.9.5.

## Question

At nodes of one true group (one Gaussian population) the edge gate opens for both children at 17 of 24 nodes
with operator "markov" and linkage "ward" and at none in the other three cells. Why is the edge test
significant there?

## The edge test, as read from the code

The contrast is the mean of the child's test rows minus the mean of the parent's test rows (each child's rows
alternate between a test half and a basis half). It is whitened with the parent's guarded within-child
covariance, divided by the square root of (1/n_child - 1/n_parent) x (1 + L / mean L), with L the edge's branch
length, and projected on k directions of the parent's PCA basis fitted on the basis half. The p-value is
chi-square(k). Tree-BH at edge alpha 0.001 tests an edge only when the edges above it were rejected.

## Setup

Both child edges of every node with at least 40 cells, all from one true group, in the 4 cells x 4 cases: 292
edges at 146 nodes (symmetric/average 96, symmetric/ward 58, markov/average 90, markov/ward 48). Four variants
of the same test per edge: (a) the tree's split with branch time (the recorded test); (b) the tree's split
without branch time; (c) 200 random splits with branch time; (d) the same 200 random splits without. A random
split divides the parent's test rows and basis rows at random into the tree split's sizes; the within-child
covariance is recomputed for it; the parent's basis and k stay. Self-checks passed: variant (a) recomputed
equals the recorded statistic, k and p at 292 of 292 edges; the recomputed covariance of the tree's split
equals the pipeline's; the 16 fits reproduce found clusters and ARI. The protocol with three hypotheses was
written before the run.

## Results

| | symmetric/average | symmetric/ward | markov/average | markov/ward |
| --- | ---: | ---: | ---: | ---: |
| Edges | 96 | 58 | 90 | 48 |
| k, median (quartiles) | 6 (2 to 30) | 9 (2 to 30) | 3 (2 to 30) | 28 (6 to 30) |
| Multiplier 1 + L / mean L, median (quartiles) | 1.00 (1.00 to 1.14) | 1.00 (1.00 to 1.11) | 1.02 (1.00 to 1.18) | 1.03 (1.00 to 1.11) |
| Statistic per k, (b) tree split, median | 1.16 | 1.82 | 2.02 | 2.60 |
| Statistic per k, (d) random splits, median | 0.99 | 0.99 | 0.99 | 1.01 |
| Raw p <= 0.001, (a) tree split with branch time | 17 (17.7%) | 18 (31.0%) | 35 (38.9%) | 36 (75.0%) |
| Raw p <= 0.001, (b) tree split without | 26 (27.1%) | 22 (37.9%) | 38 (42.2%) | 38 (79.2%) |
| Raw p <= 0.001, (c) random splits with branch time | 6 of 19,200 | 42 of 11,600 | 16 of 18,000 | 35 of 9,600 |
| Raw p <= 0.001, (d) random splits without | 9 of 19,200 (0.05%) | 59 of 11,600 (0.5%) | 32 of 18,000 (0.2%) | 52 of 9,600 (0.5%) |
| Raw p <= 0.01, (d) random splits without | 0.8% | 1.9% | 1.0% | 2.3% |
| Edges tested by Tree-BH | 12 | 12 | 24 | 42 |
| Edges significant | 0 | 0 | 2 | 34 |
| Not tested by Tree-BH, of these raw p <= 0.001 | 17 of 84 | 18 of 46 | 33 of 66 | 1 of 6 |

The statistic with branch time equals the statistic without it divided by the multiplier at 292 of 292 edges
(to 1e-11).

The largest node of each true group per case (the first node inside the group that Tree-BH reaches):

| | symmetric/average | symmetric/ward | markov/average | markov/ward |
| --- | ---: | ---: | ---: | ---: |
| Nodes | 10 | 6 | 10 | 7 |
| Share of the smaller side, median | 0.09 | 0.18 | 0.08 | 0.44 |
| Smallest raw edge p of the node, range | 0.72 to 1.0 | 0.052 to 0.98 | 3.7e-05 (2 nodes), else 0.016 to 0.21 | 2e-20 to 0.00095 |
| Both edges raw p <= 0.001 | 0 | 0 | 0 | 6 |
| Both edges significant | 0 | 0 | 0 | 6 |

Random splits by the form of the covariance (variant (d), raw p <= 0.01 and <= 0.001): where the pipeline pools
variances only (a child with fewer cells than features; 186 of the 292 edges) 0.7% to 1.0% and 0.04% to 0.1% in
all four cells; where it uses the full matrix (106 edges) symmetric/average 0.4% and 0.03%, markov/average
1.2% and 0.2%, symmetric/ward 4.2% and 1.3%, markov/ward 7.0% and 2.0%.

## Reading

- H2 supported: the cause is the selection of the split. The edge test treats the split as fixed, but the
  tree chooses it from the same cells. For the tree's split the test rejects at raw p <= 0.001 in 27% to 79%
  of the edges of one Gaussian group; for random splits of the same sizes in 0.05% to 0.5%.
- H3 not supported: the branch time is not what separates the cells. The multiplier is 1.00 to 1.03 at the
  median in all four cells and removes few rejections (26 to 17, 22 to 18, 38 to 35, 38 to 36).
- What separates the cells is the first split inside a group together with Tree-BH. In the three cells where
  the gate stays closed the first split inside a group is unbalanced (smaller side 8% to 18%) and not
  significant, so the edges below are never tested, although 20% to 50% of them have raw p <= 0.001. With
  markov/ward the first split is balanced (44%) and significant at 6 of 7 nodes, so the edges below are tested
  and rejected (34 of 42).
- So the three cells do not show that the edge test holds its level inside one Gaussian group; they are
  protected by one weak split at the top of each group.
- H1 supported in part: for a fixed split the chi-square reference holds where variances only are pooled, and
  is exceeded 4 to 7 times at the 0.01 level with the full covariance in the two ward cells (20 and 12 edges).

## Deviations

The covariance self-check stopped the first run: where a child has fewer cells than features the pipeline pools
variances only, so the recomputed covariance was reduced to its diagonal at those nodes (93 of 146), for the
tree's split and the random splits. The block on Tree-BH blocking and the largest node per group was added to
the summary after its first output.

## Not checked

Why k is 28 at the median with markov/ward and 3 to 9 in the other cells. The random splits keep the parent's
basis, k and row halves, so they test the law for a fixed split at this node, not the whole procedure. Raw
p-values in the variants, not the Tree-BH decision. One seed per case, 4 synthetic cases with identity noise
covariance, nodes with at least 40 cells, 200 draws per edge. No real data.
