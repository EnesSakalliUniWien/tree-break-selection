---
title: Gaussian Best-Flat-Cut Reference on the CHOIR Cell Lines
type: source
status: draft
updated: 2026-11-11
sources:
  - raw/inbox/choir-cell-lines/best-flat-cut-reference-real-20261111.md
  - raw/inbox/choir-cell-lines/best-flat-cut-reference-real-scripts-20261111.py.txt
tags:
  - calibration
  - real-data
  - sibling-gate
---

# Gaussian Best-Flat-Cut Reference on the CHOIR Cell Lines

## Summary

On all 48,879 CHOIR cells the opt-in `sibling_best_flat_cut_reference`, as
implemented, is worse than the split-count rule-out and the ceiling-3 rule:
Markov/Ward gives 138 clusters with ARI 0.297 (ceiling 3: 217 and 0.833). The
cause is two implemented limits, not the statistic: every split with more than
10 invariant directions is refused, and the lowest possible p-value 0.0005 is
too large for the BH step.

## Key Points

- Replay: the rule-out arm reproduces 69 and 664 clusters exactly; the same
  nodes, k and `rho` as in the ceiling-3 run.
- Markov/Ward, reference: 138 clusters, ARI 0.297, NMI 0.851, 76 lines
  recovered, 91 lines merged only, largest cluster 7,678 cells. Ceiling 3: 217,
  0.833, 0.940, 119, 33, 1,990.
- Default tree: 53 clusters, ARI 0.018; unchanged in kind.
- 190 of 1,838 Markov/Ward splits have 11 to 15 invariant directions and are
  refused, including all 103 splits with 1,000 cells or more. Ceiling 3
  accepted 77 of them (`rho` median 7.75, minimum 3.7; the table's 99% point at
  n = 800, k = 10 is 2.48).
- 17 of the 48 splits at the lowest possible p-value 0.0005 are not accepted by
  the BH step; 13 of them were accepted under ceiling 3 (`rho` 6.9 to 55.6).
- The limit k up to 10 comes from the search: the measured 99% points for
  k = 15 and 20 are lower bounds only.
- Not run: a larger reference sample, a table beyond k = 10, a restriction to
  the 10 leading invariant directions. No known no-group parts, so error
  control cannot be read from these data.

## Evidence

`raw/inbox/choir-cell-lines/best-flat-cut-reference-real-20261111.md` holds the
setup, the metric table, the split table and the breakdown; protocol, runner,
summary script and output are in
`raw/inbox/choir-cell-lines/best-flat-cut-reference-real-scripts-20261111.py.txt`.

## Links

- [[gaussian-best-flat-cut-reference]]
- [[best-flat-cut-reference-benchmark]]
- [[ceiling-rule-real-and-benchmark]]
- [[split-count-exclusion-real]]
- [[best-flat-cut-reference-k15]]
