---
title: Root Cause of the Accepted Extra Splits - the Tree Is Built from the Tested Cells
type: source
status: draft
updated: 2026-11-13
sources:
  - raw/inbox/choir-cell-lines/selected-tree-null-over-splitting-20261112.md
  - raw/inbox/choir-cell-lines/selected-tree-null-over-splitting-scripts-20261112.py.txt
  - raw/inbox/choir-cell-lines/peel-pass-default-test-20261113.md
  - raw/inbox/choir-cell-lines/peel-pass-default-test-scripts-20261113.py.txt
  - raw/inbox/choir-cell-lines/peel-pass-opt-in-20261113.md
  - tests/statistics/58_test_edge_null_calibration.py
  - benchmarks/shared/cases/gaussian.py
tags:
  - calibration
  - edge-gate
  - sibling-gate
  - selection
  - null
---

# Root Cause of the Accepted Extra Splits - the Tree Is Built from the Tested Cells

## Summary

The edge test and the sibling test are valid for a split that does not depend
on the tested cells, and the repository's null tests check only that. The
method builds the tree from the tested cells. On one Gaussian population in 30
features today's complete method finds more than 1 cluster in 46 of 60 runs,
in all four operator and linkage cells; with the opt-in
`sibling_best_flat_cut_reference` it finds 1 cluster in 80 of 80.

## Key Points

- Working reference: the two null tests in
  `tests/statistics/58_test_edge_null_calibration.py` build the tree from
  independent noise; no test builds it from the tested data.
- E1, one variable (the distance the tree is built from; 240 cells, 40
  datasets, raw edge p < 0.05 on all edges): 5 features, ward 1.4% against
  17.3%; 5, average 1.2% against 16.9%; 30, ward 2.8% against 26.1%; 30,
  average 2.6% against 15.2%.
- Same arms, parents of at least 40 cells, p <= 0.001: 0.6% against 94.3%,
  0% against 51.5%, 0.9% against 55.4%, 0.3% against 7.5%.
- E2, tie extension off (`extend_rank_through_ties`): k is 2 everywhere and
  the rejections stay (identical counts at 5 features; 59.3% to 56.6% at 30).
  The k rule decides where selected splits pass, not whether the test holds.
- E3, today's complete method on one Gaussian population, seeds 401 to 405,
  more than 1 cluster: 3 of 20 runs at 200 cells in 80 features (the
  registered null case is this size, seed 401, 1 cluster); at 30 features 46
  of 60 (symmetric/average 9 of 15, symmetric/ward 15 of 15, markov/average
  7 of 15, markov/ward 15 of 15), up to 15 clusters.
- E4, same grid with `sibling_best_flat_cut_reference=True`: 1 cluster in 80
  of 80 runs.
- Benchmark null case added the same day: `gauss_null_256n_30f_continuous`
  (256 cells, 30 features, seed 401; continuous suite 10 cases). Today 10
  clusters (default tree) and 4 (Markov/Ward); with the opt-in 1 and 1.
- 13 November, seeds 1 to 50, 30 features, today's complete method: more
  than 1 cluster in 151 of 200 fits at 256 cells and 115 of 200 at 496 cells
  (identity covariance), 197 and 194 of 200 with the unequal variances of
  the CHOIR root. The default cell (symmetric/average) at 496 cells with
  identity covariance: 0 of 50. With peels of fewer than 15 cells passing
  the edge gate ([[k-split-direction-draft]]): 825 of 860 fits against 703,
  and the default cell at 496 cells 28 of 50. With the opt-in on the same
  880 fits: 1 fit split, with and without the peel pass
  (`raw/inbox/choir-cell-lines/peel-pass-opt-in-20261113.md`).
- History: `extend_rank_through_ties` entered with `bb96e2b5` (2026-10-14),
  the basis half and its within-child whitening with `9c5cd805`
  (2026-10-18).
- Limits: 5 seeds per size; identity covariance; the sibling test alone on an
  independent tree was not run; no default was changed, so the new case
  fails for today's test; the opt-in loses real groups on the benchmark
  ([[best-flat-cut-reference-benchmark]]).

## Evidence

`raw/inbox/choir-cell-lines/selected-tree-null-over-splitting-20261112.md`
holds the four phases, the tables of E1 to E4, the root cause and the limits;
the protocol with its amendment, the three scripts and their outputs are in
`raw/inbox/choir-cell-lines/selected-tree-null-over-splitting-scripts-20261112.py.txt`.
The runs of 13 November with and without the peel pass are in
`raw/inbox/choir-cell-lines/peel-pass-default-test-20261113.md` with
`raw/inbox/choir-cell-lines/peel-pass-default-test-scripts-20261113.py.txt`.

## Links

- [[cancer-pc-over-splits]]
- [[cancer-pc-edge-test]]
- [[gaussian-best-flat-cut-reference]]
- [[best-flat-cut-reference-benchmark]]
- [[edge-selection-null-audit-20260601]]
