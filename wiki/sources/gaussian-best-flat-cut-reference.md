---
title: Gaussian Reference Table for the Best Flat Cut
type: source
status: draft
updated: 2026-11-13
sources:
  - raw/inbox/choir-cell-lines/gaussian-best-flat-cut-reference-20261111.md
  - raw/inbox/choir-cell-lines/gaussian-best-flat-cut-reference-evidence-20261111.json
  - raw/inbox/choir-cell-lines/gaussian-best-flat-cut-reference-scripts-20261111.py.txt
  - tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/best_flat_cut_reference.py
  - scripts/build_best_flat_cut_reference.py
  - raw/inbox/choir-cell-lines/peel-pass-production-20261113.md
  - raw/inbox/choir-cell-lines/table-rebuild-20261113.md
tags:
  - calibration
  - sibling-gate
  - simulation
---

# Gaussian Reference Table for the Best Flat Cut

## Summary

Under one Gaussian population the best flat cut's `rho` inside k directions
depends only on (n, k), so a table of its 99% point is a candidate reference
for the sibling test. Simulations map where that table holds: it is finite and
ordered from n = 15, settled up to k = 10 with a strong enough search, and
shows no excess rejections when k is chosen from the same data at n = 400. It
is unusable at n = 10 for k of 3 and above, only a lower bound for k = 15 and
20, and shows a small excess at n = 100 with 50 or 80 features. The recorded
checks changed no production code; the opt-in implementation is described in the
last Key Point.

## Key Points

- The best flat cut equals the best of all splits: on 100 datasets (n = 10 to
  14, k = 2 to 6) the maximum over all 2^(n-1) - 1 splits matches the exact
  flat-cut maximum to 1.1e-15. One reference therefore covers every other
  split in the same k directions.
- 99% points for k = 1 to 6 run from 2.987 to 6.173 at n = 50 and from 1.989
  to 2.302 at n = 800. They rise with k and fall with n in every cell, so
  rounding n down errs on the strict side.
- n below 50: from 5.359 (k = 1) to 62.040 (k = 6) at n = 15 and from 3.207 to
  7.685 at n = 40. At n = 10 the 99% point is 36 to 3,482 for k = 3 to 6 with a
  standard error up to 437. At n = k + 1 `rho` is infinite in some datasets.
- The search that built the first table (fixed grid, re-sorts from the single
  best direction) equals exact enumeration for n up to 40 but misses the best
  flat cut for larger n and k of 4 and above: 1,075 of 2,000 datasets at
  n = 200, k = 6. Its 99% point moves by at most 0.006, but its found `rho` is
  not invariant (changed in 138 of 200 datasets under a linear map) and no
  longer bounds other splits.
- Re-sorts from the top 50 directions plus single-point moves agree with
  300,000 directions in 60 of 60 datasets per cell up to k = 6 and change under
  a linear map in at most 2 of 200.
- k above 6 (n = 100, 200, 400): 99% points 4.566, 3.415, 2.733 for k = 8 and
  5.276, 3.692, 2.882 for k = 10 are settled to about 0.014. The values for
  k = 15 (6.879, 4.361, 3.232) and k = 20 (8.511, 4.919, 3.492) are lower
  bounds; at k = 20 the strongest search tried is still beaten in 8 to 16 of 60
  datasets.
- k chosen by `invariant_directions` from the same data, rejections at the 99%
  point: 159/20,000 with 10 or 30 features; 163/18,000 at n = 400 with 50 to
  200 features; 152/12,000 (0.0107 to 0.0148) at n = 100 with 50 or 80
  features. The grouping by n was chosen after the cells were seen.
- `invariant_directions` returns no directions when n <= p + 1, so the test is
  skipped with n - 1 features or more.
- Scope: one Gaussian population, five covariance shapes. Non-Gaussian
  populations, complete tree fits and selected descendants were not tested.
  The first run of two parts was stopped; 13 of its 42 cells were reproduced
  from the same seeds and the rest rely on its printed summary.
- Implementation (opt-in `sibling_best_flat_cut_reference`, 2026-11-11): the
  module `best_flat_cut_reference.py` compares a split's `rho` with a stored
  null of the best flat cut at (n, k) built by
  `scripts/build_best_flat_cut_reference.py` (n from 15 to 800, k = 1 to 10);
  it fails closed with p = 1 for k = 0, k above 10 and n below 15. See the
  2026-11-11 entry in [[log]] for the table check and the calibration counts.
- Table at k = 1 against the exact best cut of one direction (12 November,
  20,000 Gaussian datasets per n, [[k-split-direction-draft]]): reference
  p <= 0.01 in 1.65% at n = 15, 1.03% at 30, 1.12% at 75, 1.89% at 100, 0.66%
  at 150, 1.40% at 400. The 99% points of a cell with 1,999 draws carry
  this error. Exact 99% points from 99,999 draws: 5.343 at n = 15 (table
  4.888), 2.531 at 100 (2.439), 2.356 at 150 (2.389), 2.091 at 400 (2.073).
- The table's smallest p-value, 1/2,000, limits what Benjamini-Hochberg can
  accept per depth level; with twice the median it cost 54 accepted splits
  on the CHOIR cells.
- Since 13 November the table has 19,999 draws per cell (smallest p-value
  1/20,000), built on lisc in 2 h 45 min. At k = 1 its 99% points match the
  exact ones (5.402 at n = 15, 2.529 at 100). With the peel pass the CHOIR
  Markov/Ward fit goes from 110 recovered lines and ARI 0.8214 to 124 and
  0.8648; one Gaussian population splits in 2 of 880 fits; the benchmark is
  unchanged (`raw/inbox/choir-cell-lines/table-rebuild-20261113.md`).
- Extension on the same day ([[best-flat-cut-reference-k15]]): the table and the
  module cover k = 1 to 15 (above 10 directions n from 50); the cells for k up
  to 10 are unchanged. For k above 10 the null is the searched value's; a Ward
  tree's split exceeded the 99% point in 0 of 12,000 one-population datasets.
- Since 13 November the opt-in also lets a small peel, a node whose smaller
  child has fewer than `MIN_CELLS` = 15 leaves, pass the edge gate like a
  singleton peel ([[k-split-direction-draft]],
  `raw/inbox/choir-cell-lines/peel-pass-production-20261113.md`).

## Evidence

`raw/inbox/choir-cell-lines/gaussian-best-flat-cut-reference-20261111.md` holds
the tables, counts, intervals, search definitions and limits. The evidence JSON
holds the complete output of the n = 50 to 800 table, the enumerator check, the
same-seed re-run and the capped run, plus the first run's printed summaries.
The scripts file archives the protocol, deviations and all eight scripts.

## Links

- [[best-flat-cut-reference-benchmark]]
- [[best-flat-cut-reference-real]]
- [[best-flat-cut-reference-k15]]
- [[flat-cut-sup-count]]
- [[single-population-split-bound]]
- [[three-tests-construction-audit]]
- [[split-count-exclusion]]
- [[gaussian-root-calibration-and-cut-geometry]]
