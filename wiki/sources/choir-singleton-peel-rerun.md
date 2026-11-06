---
title: Cancer Cell-Line Rerun After the Singleton-Peel Merge
type: source
status: reviewed
updated: 2026-11-06
sources:
  - raw/inbox/choir-cell-lines/tbs-merged-rerun-20261106.md
  - raw/inbox/choir-cell-lines/tbs-merged-rerun-summary-20261106.json
  - raw/inbox/choir-cell-lines/tbs-merged-rerun-composition-20261106.csv
  - raw/inbox/choir-cell-lines/tbs-merged-rerun-fragmentation-20261106.csv
  - raw/inbox/choir-cell-lines/tbs-merged-null-rerun-20261106.md
  - raw/inbox/choir-cell-lines/tbs-merged-null-summary-20261106.csv
  - raw/inbox/choir-cell-lines/tbs-recovery-validation-null.csv
  - raw/inbox/choir-cell-lines/tbs-merged-cell-lines-pvalues-20261106.md
  - raw/inbox/choir-cell-lines/tbs-merged-cell-lines-pvalues-script.py.txt
  - raw/inbox/choir-cell-lines/tbs-merged-cell-line-breakdown-20261106.csv
  - raw/inbox/choir-cell-lines/tbs-merged-sibling-tests-20261106.csv
  - raw/inbox/choir-cell-lines/tbs-merged-cluster-stops-20261106.csv
  - raw/inbox/choir-cell-lines/tbs-inflation-debug-20261106.md
  - raw/inbox/choir-cell-lines/tbs-inflation-factors-20261106.csv
  - tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflation_correction/inflation_adjusted_sibling_tests.py
  - raw/inbox/choir-cell-lines/tbs-run.md
  - raw/inbox/choir-cell-lines/tbs-run-script.py.txt
  - tree_break_selection/hierarchy_analysis/statistics/multiple_testing/tree_bh/correction.py
  - tree_break_selection/hierarchy_analysis/tree_decomposition.py
  - tree_break_selection/hierarchy_analysis/decomposition/gates/gate_evaluator.py
  - tests/core/test_cancer_pc_surrogate.py
  - CHANGELOG.md
tags:
  - cancer
  - traversal
  - multiple-testing
---

# Cancer Cell-Line Rerun After the Singleton-Peel Merge

## Summary

Merge commit `50455316` brought the singleton-peel traversal fix and the
sibling inflation cap into `dev`. Rerunning family-diffusion TBS on the
unchanged 48,879-cell inputs gives 70 clusters instead of one. The tree is
byte-identical; only the gates and traversal changed. Cell-line recovery stays
poor because one cluster holds 58.5% of the cells. On the rerun count nulls,
the production baseline no longer meets the frozen 0.01 target.

## Key Points

- A singleton peel (one leaf child, one internal child) no longer ends the
  descent. Tree-BH runs one BH on a peel's family merged with its internal
  child's family, and the edge gate and sibling frontier pass through it.
  The inflation factor is capped at the family baseline beyond the largest
  calibration parent.
- The cancer root is still `L43204` versus 48,878 cells. Both root edges stay
  non-significant (BH p 0.916734 and 0.434736), but traversal now visits 145
  nodes, with 45 splits, 27 pass-throughs and 148 significant edges.
- ARI against the 190 cell lines rises from 0 to 0.027593 and NMI from 0 to
  0.548219; published `CHOIR_1` labels score ARI 0.664014. The largest cluster
  holds 28,611 cells from 144 lines. 27 lines are recovered at the 90%/90%
  threshold, and 46 of 70 clusters are at least 90% one line.
- On the merged code, the `cancer_pc_hub_disconnected` surrogate also passes
  through its 1/495 singleton root and returns six clusters (ARI 0.9982,
  previously one cluster and ARI 0). `tests/core/test_cancer_pc_surrogate.py`
  now pins six clusters and one pass-through for the hub, none for its control.
- Peel leaves whose edges are not significant join the nearest cluster below
  the peel, so the output is no longer an exact tree cut.
- The 600 frozen count nulls were rerun with all four variants on the same
  hashed datasets. Symmetric/average false partitions go from 0 to 2 of 300
  low-depth nulls (upper one-sided 95% bound 0.020836) and from 262 to 264 of
  300 high-depth nulls. Both new low-depth splits pass through a 1/399
  singleton root. Markov/average rises from 34 to 41 low-depth splits; the
  other five variant-model results are unchanged.
- By cell line: 27 recovered, 30 separated but sharing a cluster with lines of
  the same cancer type (16 skin lines in one cluster, 9 breast lines in
  another), 13 split (mostly colorectal) and 120 dominated by the giant cluster,
  including every head-and-neck line. `CHOIR_1` recovers 79 lines under the same
  definition; 9 are recovered by both.
- The p-values are internally consistent: all lie in [0, 1], BH values never
  fall below raw values, and BH-small edges left non-significant follow Tree-BH's
  scaled family levels. 81 edge and 7 sibling p-values underflow to zero.
- Of the 45 sibling splits, 36 have BH p at most 1e-10 and 5 lie between 0.001
  and 0.01. 34 use the capped `family_baseline` inflation estimator. Recording
  the factors and rerunning without the cap turns 23 of them into
  pass-throughs, but the 70 clusters are identical, so the partition does not
  depend on the cap ([[sibling-inflation-debugging]]).
- 63 of 70 clusters stop at a non-significant child edge. Clusters stopping at
  a split whose smaller child has at most 10 cells hold 65% of the cells. The
  giant cluster's core splits 4 versus 28,606 cells, with child-edge BH p
  0.561712 and 0.088406; pass-through applies only to one-leaf peels.

## Evidence

`raw/inbox/choir-cell-lines/tbs-merged-rerun-20261106.md` records held-fixed
inputs, input and output hashes, the before/after table, recovery definitions,
the post-output process stall and scope limits. The summary, composition and
fragmentation tables are copied beside it. The full ordered suite passed 1,038
tests on the merged working tree after the surrogate test update. The
singleton-peel change, its 122-case benchmark effect and its 120 reseeded null
runs are described in `CHANGELOG.md`.

`raw/inbox/choir-cell-lines/tbs-merged-null-rerun-20261106.md` records the
count-null rerun: unchanged protocol and scripts, dataset hash matching, the
before/after table and the singleton-root mechanism of the new splits. With
300 nulls any single false partition exceeds the 0.01 bound, so 2/300 shows
the target is no longer established, not that the rate clearly exceeds 1%.
Signal mixtures were not rerun. Neither result establishes selected-tree
calibration or cancer recovery.

`raw/inbox/choir-cell-lines/tbs-merged-cell-lines-pvalues-20261106.md` records
the cell-line breakdown, the p-value checks and the stopping analysis. Its
script regenerates the per-line, per-sibling-test and per-cluster tables from
the local run outputs without rerunning TBS.

## Links

- [[choir-cell-lines-data]]
- [[choir-cell-lines-recovery-plan]]
- [[sibling-inflation-debugging]]
- [[cancer-pc-surrogate-evidence]]
- [[root-conditional-sampler]]
