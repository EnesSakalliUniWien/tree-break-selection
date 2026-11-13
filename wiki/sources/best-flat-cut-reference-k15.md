---
title: Best-Flat-Cut Reference Extended to 11 to 15 Invariant Directions
type: source
status: draft
updated: 2026-11-12
sources:
  - raw/inbox/choir-cell-lines/best-flat-cut-reference-k15-20261111.md
  - raw/inbox/choir-cell-lines/best-flat-cut-reference-k15-scripts-20261111.py.txt
  - tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/best_flat_cut_reference.py
  - scripts/build_best_flat_cut_reference.py
tags:
  - calibration
  - sibling-gate
  - real-data
  - simulation
---

# Best-Flat-Cut Reference Extended to 11 to 15 Invariant Directions

## Summary

The opt-in `sibling_best_flat_cut_reference` now covers 1 to 15 invariant
directions (`MAX_DIRECTIONS` = 15; for more than 10 directions at least 50
cells). The cells for 1 to 10 directions are bit-identical to the earlier
table. Under one Gaussian population the split of a Ward tree never exceeded
the 99% point at k = 11, 13, 15 (0 of 12,000 datasets). On the CHOIR cells the
Markov/Ward fit goes from ARI 0.2975 to 0.8219 on the unchanged tree.

## Key Points

- New cells: n from 50 to 800, k = 11 to 15, B = 1,999, same search and
  seeds; 99% points from 11.089 to 17.339 at n = 50 and from 2.531 to 2.671
  at n = 800; they rise with k and fall with n in every cell.
- For k above 10 the stored null is the null of the searched value. Stronger
  searches beat the table's search in 30 to 70% of datasets by a median of
  0.6 to 2.6%; the best of them exceeds the stored 99% point in 25 of 2,700
  datasets (0.93%), 6 of 300 at n = 800, k = 15.
- The tree's split is far below the reference under one population: median
  ratio to the 99% point 0.23 to 0.30, maximum 0.50; 0 rejections also with
  k chosen by `invariant_directions` in 30 features.
- Real data, unchanged tree: all 190 splits with 11 to 15 directions are
  tested, 92 at p <= 0.01 and 88 accepted by BH; 230 clusters, ARI 0.8219,
  116 lines recovered, largest cluster 1,990 (before: 138, 0.2975, 76,
  7,678; ceiling 3: 217, 0.8329, 119, 1,990).
- Tree after the sprs with the extended reference: 316 clusters, ARI 0.7376,
  merged-only lines 19 against 38, split-only lines 58 against 21.
- Tests: 13 pass in the three reference test files, 534 in
  `tests/statistics`; ruff passes.
- Limits: one Gaussian population; Ward root split only; the strongest
  search is still a search; one real dataset and tree. The 13 benchmark
  cases were rerun on 12 November and are unchanged because none has more
  than 5 invariant directions ([[best-flat-cut-reference-benchmark]]).

## Evidence

`raw/inbox/choir-cell-lines/best-flat-cut-reference-k15-20261111.md` holds the
change, the table of 99% points, the validity counts, the real-data table,
the reading and the limits; protocol, deviations, check scripts and outputs
are in
`raw/inbox/choir-cell-lines/best-flat-cut-reference-k15-scripts-20261111.py.txt`.

## Links

- [[gaussian-best-flat-cut-reference]]
- [[best-flat-cut-reference-real]]
- [[best-flat-cut-reference-benchmark]]
- [[spr-refit]]
- [[best-flat-cut-divisive]]
- [[ceiling-rule-real-and-benchmark]]
