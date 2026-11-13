# Draft definition of k that does not count variance inside one population

Written and checked on 2026-11-12, Europe/Vienna. Local folder:
`results/analyses/k_split_direction_draft_20261112` (git-ignored); the protocol with the definition, the script
and its output are in `k-split-direction-draft-scripts-20261112.py.txt`. Follow-up of
`line-gene-programs-20261112.md` and step A2 of `plan-non-heuristical-statistic-20261110.md`. The draft exists in
the analysis folder only; no production code was changed.

## Why a different definition

k today is a count of eigenvalues of the node's covariance above an edge. On real cells one population has
real directions of variation, so that count is 5 to 12 inside one cell line and cannot tell one population from
several. The reference does not need such a count: for one Gaussian population with any covariance, `rho` of the
best flat cut in k fixed directions depends only on (n, k), and any split's `rho` is at most the best flat cut's.
What the reference needs is that the directions do not depend on the cells that are tested.

## The draft definition

For a node with cells X and the tree's split into left and right:

1. Halves. The node's cells are divided into a basis half and a test half, alternating within each side in the
   order of the cells (first cell of a side to the test half), as the edge test does. The division ignores the
   data values.
2. Direction, from the basis half only; one direction w, so k = 1.
   - Variant A uses the tree's split: w = S^-1 (mean of the left basis cells minus mean of the right basis
     cells), S the pooled within covariance of the two sides in the basis half.
   - Variant B does not use the tree's split: the same formula for the two groups of a 2-means clustering of the
     basis half.
   - With no more than p + 2 basis cells, w is the difference of the two means without S^-1.
3. Test, on the test half only: `rho` of the tree's split of the test cells along w, against the reference at
   (number of test cells, 1).
4. Fail closed (p = 1): fewer than 15 test cells, or a side without basis cells.

Under one Gaussian population variant B is exact up to the Monte Carlo error of the table: w is fixed given the
basis half, the test cells are independent of it, and the tree's split is one of the splits the best flat cut
covers. In variant A the test cells influence w through the tree's split, so its level is measured, not derived.

What the definition does not use: no eigenvalue count, no noise band, no assumption that one population is
spherical. Variation inside one population enters only if it separates the two sides along w.

## Check on simulations (protocol written before the run)

Tree: Ward linkage on the data itself. 1,000 datasets per setting. Datasets with p <= 0.01:

| Setting | Node | Variant A | Variant B | Current rule |
| --- | --- | ---: | ---: | ---: |
| One population, identity, 64 cells, 30 features | root / larger child | 1 / 3 | 0 / 0 | 0 / 0 |
| One population, identity, 256 cells | root / larger child | 0 / 0 | 0 / 0 | 0 / 0 |
| One population, identity, 1,000 cells | root / larger child | 0 / 0 | 0 / 0 | 0 / 0 |
| One population, variances of the CHOIR root, 256 cells | root / larger child | 0 / 0 | 0 / 0 | 0 / 0 |
| One population, variances of the CHOIR root, 1,000 cells | root / larger child | 0 / 0 | 0 / 0 | 0 / 0 |
| One population, 30 cells in 40 features | root | 13 | 0 | 0 |
| Two groups, 256 cells, 30 features, `rho` 2 | root | 0 | 0 | 0 |
| Two groups, `rho` 3 | root | 63 | 39 | 164 |
| Two groups, `rho` 5 | root | 914 | 882 | 996 |
| Two groups, `rho` 10 | root | 1,000 | 1,000 | 1,000 |
| Two groups, unequal within-group variances, `rho` 5 | root | 0 | 0 | 0 |
| Two groups, unequal within-group variances, `rho` 10 | root | 535 | 249 | 543 |
| Two groups of 15 cells in 40 features, `rho` 50 | root | 1,000 | 1,000 | 0 |

In the settings with unequal within-group variances the tree's root split equals the two groups in 0 and 35 of
1,000 datasets, so those rows measure the tree as much as the test.

Predictions:

| Prediction | Result |
| --- | --- |
| Q1 variant B: p <= 0.01 in at most 2.0% of the datasets in every setting with one population | held (0 of 11,000) |
| Q2 variant A exceeds the level in at least one setting with one population | held by the count only: 13 of 1,000 at 30 cells in 40 features, within sampling error of 1%; at most 3 of 1,000 elsewhere |
| Q3 at `rho` 5 and 10 both variants reach p <= 0.01 in at least 90% | not held: variant B 882 of 1,000 at `rho` 5 |
| Q4 with unequal within-group variances both variants beat the current rule | not held: 535 and 249 against 543 |

## Reading

- The draft does what it is meant to do: k is 1 by definition and no variance inside one population is counted;
  with one population it rejects in 0 of 11,000 datasets (variant B) and in at most 13 of 1,000 (variant A).
- It works with more features than cells (1,000 of 1,000 against 0 for the current rule).
- It costs power with more cells than features, because half of the cells are not tested: at `rho` 5 it accepts
  88% to 91% against 99.6%, at `rho` 3 4% to 6% against 16%.
- It gives no advantage when the within-group variances are unequal; there the tree's split is the limit.
- Variant A is closer to the split under test and has more power than variant B; its level is not derived.

## Not checked

Real data and the benchmark cases; the method's own diffusion tree; nodes deeper than the root's larger child;
more than one direction; another way to form the halves. The reference's population class (one Gaussian
population) is unchanged, so skewed or bimodal variation inside one population, for example along the cell
cycle, is not addressed by this definition.
