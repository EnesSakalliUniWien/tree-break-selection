# Cell lines and p-values of the merged cancer rerun

Analysis on 2026-11-06 of
`results/analyses/choir_cell_lines_tbs_merged_20261106/` (see
`tbs-merged-rerun-20261106.md`). Tables were produced by
`tbs-merged-cell-lines-pvalues-script.py.txt` from `cluster_assignments.csv`,
`annotations.csv` and `linkage.npy`; no TBS computation was rerun.

## Cell lines

A line is recovered when at least 90% of its cells lie in its dominant cluster
and that cluster is at least 90% that line.

| Outcome | Cell lines | Cells |
| --- | ---: | ---: |
| Recovered | 27 | 10,220 |
| Separate from the giant cluster, sharing a cluster with other lines | 30 | 7,151 |
| Split (dominant cluster holds under 90% of the line) | 13 | 2,952 |
| Dominant cluster is the giant 28,611-cell cluster | 120 | 28,556 |

- Recovered lines include six lung lines (NCIH2110, NCIH460, NCIH1944,
  NCIH446, NCIH522, NCIH727) and three ovarian lines (OAW28, ONCODG1, OVSAHO).
- Shared clusters group lines by cancer type: cluster 68 holds 16 skin lines
  and cluster 33 holds 9 breast lines, including MCF7, T47D and BT474. Smaller
  shared clusters hold 2–3 liver, lung or ovarian lines.
- Split lines are mostly colorectal (8) and gastric (2); their dominant clusters
  hold 28–42% of their cells.
- The giant cluster contains every head-and-neck line (11), 11 of 12 brain,
  31 of 40 lung and 9 of 11 pancreatic lines.
- Published `CHOIR_1` labels recover 79 of 190 lines under the same definition;
  9 lines are recovered by both.

`tbs-merged-cell-line-breakdown-20261106.csv` lists every line.

## P-value checks

- All stored p-values lie in [0, 1]. Edge BH values are never below the raw
  values. 81 edge BH p-values and 7 sibling p-values are exactly zero
  (floating-point underflow); decisions are unaffected but their ordering is lost.
- 328 edges were tested and 148 are significant; 97,428 were ancestor-blocked.
  Some non-significant edges have BH p below the 0.001 edge alpha. This is by
  design: Tree-BH compares each family's BH values with the base alpha scaled by
  ancestor rejection shares (`multiple_testing/tree_bh/correction.py`).
- 78 sibling tests were run: 45 different (the 45 splits) and 33 same. Of the
  45 splits, 36 have BH p at most 1e-10, 4 between 1e-10 and 0.001, and 5
  between 0.001 and 0.01 (maximum 0.00923). Same decisions have BH p from
  0.0112 upward, median 0.7153.
- 34 of the 45 splits, and 39 of the 78 tests, use the capped
  `family_baseline` inflation estimator added by the merge; the rest use
  `local_kernel`. The fitted inflation factors are not saved:
  `Sibling_Test_Statistic` is already divided by reference scale and inflation
  factor (`inflation_adjusted_sibling_tests.py`). Whether those splits depend
  on the cap is therefore not determined.

## Why clusters stop

For each cluster, the largest subtree wholly inside it is taken as its core;
passed-through peel leaves can join a cluster from outside that subtree.

- 63 of 70 clusters stop because a child edge of the core is not significant,
  3 stop at a same sibling decision and 4 are single leaves.
- Among the 63 edge-gate stops, the smaller child has 1 cell in 7, 2–10 cells in
  25, 11–100 in 26 and over 100 in 5. Clusters stopping at a split with a
  smaller child of at most 10 cells hold 31,982 cells (65%).
- The giant cluster's core `N97667` (28,610 cells plus one peel leaf) splits 4
  versus 28,606 cells. Its child edges have BH p 0.561712 and 0.088406, so the
  descent stops. The singleton-peel change passes only through one-leaf peels.

`tbs-merged-sibling-tests-20261106.csv` and
`tbs-merged-cluster-stops-20261106.csv` hold the per-test and per-cluster rows.

## Limits

These are selected-tree p-values without a calibrated reference law. On the
merged code the production baseline no longer meets the frozen count-null
target (`tbs-merged-null-rerun-20261106.md`), so individual splits and stops
are descriptive.
