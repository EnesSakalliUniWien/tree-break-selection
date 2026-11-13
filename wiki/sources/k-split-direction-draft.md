---
title: Draft Definition of k that Does Not Count Variance inside One Population
type: source
status: draft
updated: 2026-11-13
sources:
  - raw/inbox/choir-cell-lines/k-split-direction-draft-20261112.md
  - raw/inbox/choir-cell-lines/k-split-direction-draft-scripts-20261112.py.txt
  - raw/inbox/choir-cell-lines/k-split-direction-runs-20261112.md
  - raw/inbox/choir-cell-lines/k-split-direction-runs-scripts-20261112.py.txt
  - raw/inbox/choir-cell-lines/k-split-direction-level-20261112.md
  - raw/inbox/choir-cell-lines/k-split-direction-level-scripts-20261112.py.txt
  - raw/inbox/choir-cell-lines/variant-a-check-20261112.md
  - raw/inbox/choir-cell-lines/variant-a-check-scripts-20261112.py.txt
  - raw/inbox/choir-cell-lines/variant-a-fixed-20261112.md
  - raw/inbox/choir-cell-lines/variant-a-fixed-scripts-20261112.py.txt
  - raw/inbox/choir-cell-lines/default-tree-third-version-20261112.md
  - raw/inbox/choir-cell-lines/default-tree-third-version-scripts-20261112.py.txt
  - raw/inbox/choir-cell-lines/peel-pass-20261112.md
  - raw/inbox/choir-cell-lines/peel-pass-scripts-20261112.py.txt
  - raw/inbox/choir-cell-lines/peel-pass-row-orders-20261113.md
  - raw/inbox/choir-cell-lines/peel-pass-row-orders-scripts-20261113.py.txt
  - raw/inbox/choir-cell-lines/peel-pass-default-test-20261113.md
  - raw/inbox/choir-cell-lines/peel-pass-default-test-scripts-20261113.py.txt
  - raw/inbox/choir-cell-lines/peel-pass-opt-in-20261113.md
  - raw/inbox/choir-cell-lines/peel-pass-opt-in-scripts-20261113.py.txt
  - raw/inbox/choir-cell-lines/peel-pass-production-20261113.md
  - raw/inbox/choir-cell-lines/peel-pass-production-scripts-20261113.py.txt
  - raw/inbox/choir-cell-lines/table-rebuild-20261113.md
  - tree_break_selection/tree/topology.py
  - tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/best_flat_cut_reference.py
tags:
  - invariant-directions
  - sibling-gate
  - calibration
  - draft
---

# Draft Definition of k that Does Not Count Variance inside One Population

## Summary

A draft, outside production: k is 1. The node's cells are divided into a
basis half and a test half; one direction is computed from the basis half
only; `rho` of the tree's split of the test half along that direction is
compared with the reference at (test cells, 1). No eigenvalue count and no
spherical population are used. With one Gaussian population it rejects in 0
of 11,000 simulated datasets when the direction ignores the tree's split.
With the direction from the tree's split (variant A) the continuous benchmark
suite reaches 0.8302 against 0.6950 (9 cases), with the null cases at 1
cluster. On the CHOIR cells variant A depends on the order of the rows: ARI
0.8043 to 0.8397 over five orders against 0.8219 for the current rule. An
independent check found one defect (an uncaught error) and this dependence.
A third version fixes both: 20 random divisions per node, twice the median
p-value and an exact null of 99,999 draws for one direction. On the CHOIR
cells it gives ARI 0.8349 to 0.8411 over five orders and splits 0 of 800
fits of one Gaussian population.

## Key Points

- Reason: the reference needs directions that do not depend on the tested
  cells, not a count of variance directions; on real cells that count is 5 to
  12 inside one cell line ([[line-non-sphericity]]).
- Halves alternate within each side in the order of the cells, as in the edge
  test. Variant A takes the direction from the tree's split on the basis
  half, variant B from a 2-means clustering of the basis half.
- Variant B is exact under one Gaussian population up to the table's Monte
  Carlo error: the direction is fixed given the basis half and the best flat
  cut covers every split of the test cells.
- One population, 1,000 datasets per setting, p <= 0.01: variant B 0 in all
  11 settings; variant A at most 3, and 13 with 30 cells in 40 features;
  current rule 0.
- Two groups, 256 cells, 30 features: `rho` 3 gives 63 (A), 39 (B), 164
  (current); `rho` 5 gives 914, 882, 996; `rho` 10 gives 1,000 in all.
- Two groups of 15 cells in 40 features: 1,000 and 1,000 against 0 for the
  current rule.
- Unequal within-group variances: no advantage (535 and 249 against 543);
  the tree's split is the limit there.
- Predictions: 2 of 4 held; power at `rho` 5 for variant B (88.2%) and the
  advantage with unequal variances did not hold.
- Benchmark with the opt-in, 14 cases (same day): variant A finds 4, 3 and
  2 clusters in gauss_clear_medium, gauss_moderate_3c and
  gauss_dense_signal_highd (current rule 3, 2, 1); continuous mean ARI 0.8302
  (current 0.6950, prototype rule 0.8880); variant B 0.7768; no case above
  its true groups; one-population grid 80 of 80 runs at 1 cluster for both
  variants, smallest p-value 0.062 and 0.0735.
- CHOIR cells, Markov/Ward tree, same 1,838 calls: variant A 218 clusters,
  ARI 0.8311, 120 lines recovered, 151 splits accepted (current rule 230,
  0.8219, 116, 134); ARI between the two assignments 0.9817. Variant B 153
  clusters, ARI 0.5830.
- Variant B loses power where a node holds more than two groups: its
  direction does not follow the tree's split (`rho` 0.31 against 1,213 at one
  benchmark node).
- The working variant A has a measured, not a derived level. Measured on
  1,960 fits of one Gaussian population with the method's own trees (11
  settings, identity and unequal variances): p <= 0.01 at 109 of 31,482
  tested nodes with at least 30 cells (0.35%) and at 34 of 12,302 tested
  calls (0.28%); the current rule at 25 of 24,135 nodes.
- By test cells: 0.86% at 15 to 19, 0.49% at 20 to 29, 0.12% at 30 to 49,
  none at 100 and more. Identity 0.34%, unequal variances 0.36%.
- Complete method: more than 1 cluster in 12 of 1,960 fits; those fits find
  3 to 9 clusters from one accepted node of 29 to 65 cells, because the
  nodes above it are passed through.
- Independent check (same day), defect: a feature that is constant in the
  basis half makes the within covariance singular and the draft raises
  `LinAlgError`, which stops the fit; the current rule gives p = 1.0 there.
- Independent check, row order: the halves alternate in row order. CHOIR,
  five orders of the same cells: 196 to 225 clusters, ARI 0.8043 to 0.8397
  (current rule 230, 0.8219, no dependence); 39 of 170 ever-accepted nodes
  are accepted in only some of four orders. Benchmark: found clusters and
  ARI unchanged in 28 of 28 rows under 5 orders.
- One population that is not Gaussian, fits with more than 1 cluster of 160:
  skewed 3 (current rule 1), heavy tails 0 (0), uniform 1 (2), closed curve
  143 (158).
- CHOIR, accepted splits at nodes with at least 90% of the cells from one
  line: variant A 3 of 151, current rule 8 of 134.
- The benchmark mean 0.8302 against 0.6950 covers the 9 continuous cases of
  11 November; over all 10 continuous cases it is 0.8472 against 0.7255.
- Fixed versions (same day). Second version: pseudo-inverse for the
  direction (no error), 20 random divisions per node, twice the median of
  their p-values (Meinshausen, Meier and Buehlmann 2009,
  https://doi.org/10.1198/jasa.2009.tm08647). CHOIR over five orders: 133
  to 135 clusters, ARI 0.5371 to 0.5374.
- Cause of that loss: the table's smallest p-value. 151 calls sit at
  0.001 = 2 x 0.0005 and Benjamini-Hochberg accepts 97 of them; 46 of the 54
  splits lost against variant A have p = 0.001.
- Third version: the same with an exact null of the best cut of one
  direction, 99,999 Gaussian datasets per grid n (smallest node p-value
  0.00002). CHOIR over five orders: 206 to 215 clusters, ARI 0.8349 to
  0.8411, 126 to 127 lines recovered, 29 to 31 merged only (current rule
  230, 0.8219, 116, 38); 149 nodes accepted in every order, 9 in only some.
- Third version, level: 0 of 800 fits of one Gaussian population with more
  than 1 cluster and 0 of 7,249 tested calls with p <= 0.01 (variant A: 8
  fits, 20 calls); benchmark equal to variant A in 28 of 28 rows.
- Exact 99% points for one direction against the table: 5.343 against 4.888
  at n = 15, 2.531 against 2.439 at n = 100.
- Default tree of the CHOIR cells: current rule 64 clusters, third version
  63 in three row orders, ARI 0.0275 for both, largest cluster 28,611 cells
  from 144 lines. The tree is a chain of small groups peeled from one
  remainder; the sibling test is called at 90 nodes (1,838 on Markov/Ward)
  and not at all below a peel of 4 cells above 28,606. The statistic is
  not the limit there.
- The stop, read from the edge columns: at 4 against 28,606 both edges
  have statistic 44.4 with k = 30 (raw p 0.044, edge alpha 0.001); Tree-BH
  then closes the subtree, where edge statistics are 1,044 to 22,769. Of
  155 nodes with at least 1,000 cells, 5 fail at peels of 1 to 14 cells and
  101 are never tested below them, 94 with raw p <= 0.001.
- Peel pass, at run time only: a node whose smaller child has fewer than 15
  cells passes the edge gate like a peel of one cell. Default tree: 402
  clusters (third version) and 432 (current rule), ARI 0.6172 and 0.6210
  against 0.0275, largest cluster 1,990 against 28,611, 42 lines recovered
  against 29, sibling test called at 808 nodes against 90. Markov/Ward tree,
  third version, five row orders: 210 to 217 clusters and ARI 0.8642 to
  0.8706 against 206 to 215 and 0.8349 to 0.8411; the gain within one order
  is 0.0283 to 0.0305. Benchmark equal in 28 of 28 rows; one Gaussian population at
  1 cluster in 880 of 880 fits, 0 of 8,723 tested calls with p <= 0.01. The
  edge gate itself is unchanged.
- The peel pass with today's default sibling test (no opt-in) on one
  Gaussian population in 30 features: more than 1 cluster in 825 of 860
  fits against 703 without it; mean found clusters 9.68 against 7.61. The
  increase is in the cells with average linkage: the default cell at 496
  cells goes from 0 of 50 fits split to 28 of 50. The peel pass is safe only
  with a sibling rule that holds its level on the method's own trees
  ([[selected-tree-null-over-splitting]]).
- The peel pass with the opt-in's current rule (table, k up to 15) on the
  same 880 fits: found clusters equal with and without it in 880 of 880; 1
  fit is split in both arms (496 cells, unequal variances, 4 clusters). The
  over-splitting belongs to today's default test, not to the peel pass.
- In production since 13 November: with `sibling_best_flat_cut_reference` the
  gate pipeline sets `small_peel_leaf_limit` to the reference's `MIN_CELLS`
  (15) and `peel_continuations` in `tree_break_selection/tree/topology.py`
  maps small peels too; without the opt-in nothing changes. Production
  equals the trial on the default tree (48,879 of 48,879 cells, 432 clusters,
  ARI 0.6210), on 880 of 880 one-population fits and on 26 of 26 benchmark
  rows. On the Markov/Ward tree the opt-in's current rule goes from 230
  clusters, ARI 0.8219 and 116 recovered lines to 223, 0.8214 and 110: of
  the calls at the table's smallest p-value 129 of 144 are accepted against
  133 of 140. Call by call, the 1,819 calls in both fits have the same
  p-value; 7 nodes at p = 0.0005 (84 to 1,022 cells) are refused only with
  the peel pass, and 380 of the 391 added calls have p = 1. The loss is in
  the Benjamini-Hochberg step: with a smallest p-value of 0.0005 a call can
  be rejected at 0.01 only where 5% of its family is at that value. The
  third version is still outside production.
- Fixed on 13 November by a table with 19,999 draws per cell (smallest
  p-value 1/20,000): Markov/Ward 244 clusters, ARI 0.8648, 124 lines
  recovered; default tree ARI 0.6255; 2 of 880 one-population fits split;
  benchmark unchanged ([[gaussian-best-flat-cut-reference]]).
- Sprs under the third version ([[spr-refit]]): the tree after the sprs of
  11 November gives ARI 0.7691 against 0.8359; only 32 of 1,302 refused
  nodes have a passing best cut along their direction.
- Limits: the level of variant A is measured for Gaussian populations only;
  a node needs about 30 cells; one direction; half of the cells untested; one real dataset and
  tree; the population class of the reference is unchanged.

## Evidence

`raw/inbox/choir-cell-lines/k-split-direction-draft-20261112.md` holds the
reason, the definition, the table of the check, the four predictions, the
reading and the limits; the protocol, the script and its output are in
`raw/inbox/choir-cell-lines/k-split-direction-draft-scripts-20261112.py.txt`.
The runs on the benchmark cases and the CHOIR cells are in
`raw/inbox/choir-cell-lines/k-split-direction-runs-20261112.md` with protocol,
plug-in, scripts and outputs in
`raw/inbox/choir-cell-lines/k-split-direction-runs-scripts-20261112.py.txt`.
The level of variant A on one-population data is in
`raw/inbox/choir-cell-lines/k-split-direction-level-20261112.md` with protocol,
scripts, summary and spot check in
`raw/inbox/choir-cell-lines/k-split-direction-level-scripts-20261112.py.txt`.
The independent check and its verification are in
`raw/inbox/choir-cell-lines/variant-a-check-20261112.md` with brief, summaries
and scripts in
`raw/inbox/choir-cell-lines/variant-a-check-scripts-20261112.py.txt`.
The fixed versions, their looked-up sources and results are in
`raw/inbox/choir-cell-lines/variant-a-fixed-20261112.md` with protocol,
rules, the builder of the exact null and outputs in
`raw/inbox/choir-cell-lines/variant-a-fixed-scripts-20261112.py.txt`.
The run on the default tree is in
`raw/inbox/choir-cell-lines/default-tree-third-version-20261112.md` with
protocol, scripts and outputs in
`raw/inbox/choir-cell-lines/default-tree-third-version-scripts-20261112.py.txt`.
The peel pass is in `raw/inbox/choir-cell-lines/peel-pass-20261112.md` with
protocol, scripts and outputs in
`raw/inbox/choir-cell-lines/peel-pass-scripts-20261112.py.txt`; its five
row orders are in
`raw/inbox/choir-cell-lines/peel-pass-row-orders-20261113.md` with
`raw/inbox/choir-cell-lines/peel-pass-row-orders-scripts-20261113.py.txt`;
today's default test with the peel pass is in
`raw/inbox/choir-cell-lines/peel-pass-default-test-20261113.md` with
`raw/inbox/choir-cell-lines/peel-pass-default-test-scripts-20261113.py.txt`;
the opt-in's current rule with the peel pass is in
`raw/inbox/choir-cell-lines/peel-pass-opt-in-20261113.md` with
`raw/inbox/choir-cell-lines/peel-pass-opt-in-scripts-20261113.py.txt`.
The production change and its verification are in
`raw/inbox/choir-cell-lines/peel-pass-production-20261113.md` with
`raw/inbox/choir-cell-lines/peel-pass-production-scripts-20261113.py.txt`.

## Links

- [[spr-refit]]
- [[non-heuristical-statistic-plan]]
- [[selected-tree-null-over-splitting]]
- [[line-non-sphericity]]
- [[k-margin-prototype]]
- [[gaussian-best-flat-cut-reference]]
- [[exact-gate-sibling-test]]
