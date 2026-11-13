---
title: Why the Edge Test Is Significant inside One Gaussian Group on the cancer_pc Cases
type: source
status: draft
updated: 2026-11-12
sources:
  - raw/inbox/choir-cell-lines/cancer-pc-edge-test-20261112.md
  - raw/inbox/choir-cell-lines/cancer-pc-edge-test-scripts-20261112.py.txt
  - raw/inbox/choir-cell-lines/cancer-pc-edge-k-20261112.md
  - raw/inbox/choir-cell-lines/cancer-pc-edge-k-scripts-20261112.py.txt
tags:
  - benchmark
  - calibration
  - edge-gate
  - selection
---

# Why the Edge Test Is Significant inside One Gaussian Group on the cancer_pc Cases

## Summary

Inside one true group the edge test rejects because the tree chooses the split
from the same cells and the test treats the split as fixed. For the tree's
split it rejects at raw p <= 0.001 in 27% to 79% of the edges, for random
splits of the same sizes in 0.05% to 0.5%. The branch time does not separate
the four operator and linkage cells; the first split inside a group and
Tree-BH do.

## Key Points

- 292 edges at 146 nodes of one true group (at least 40 cells), 4 cells of
  operator and linkage x 4 cases; the recorded test is reproduced at 292 of
  292 edges.
- Tree's split without branch time, raw p <= 0.001: symmetric/average 27.1%,
  symmetric/ward 37.9%, markov/average 42.2%, markov/ward 79.2%. Statistic
  per k, median: 1.16, 1.82, 2.02, 2.60; random splits 0.99 to 1.01.
- Random splits without branch time, raw p <= 0.001: 0.05%, 0.5%, 0.2%, 0.5%.
- Branch time: the multiplier 1 + L / mean L has median 1.00 to 1.03 in all
  four cells; with it the rejections go from 26, 22, 38, 38 to 17, 18, 35, 36.
- Tree-BH tests 12, 12, 24 and 42 of the edges; significant 0, 0, 2 and 34.
  Of the edges it does not test, 17 of 84, 18 of 46 and 33 of 66 have raw
  p <= 0.001.
- The largest node of each group: smaller side 0.09, 0.18, 0.08 and 0.44 at
  the median; both edges raw p <= 0.001 at 0, 0, 0 and 6 of 7 nodes.
- So three cells are protected by one weak, unbalanced split at the top of
  each group, not by a valid level of the edge test.
- Fixed split, form of the covariance: variances only (186 edges) 0.7% to
  1.0% at 0.01; full matrix in the two ward cells 4.2% and 7.0% (20 and 12
  edges).
- What decides k (follow-up, 3,241 nodes recorded): the signal count is 0 or
  1 at all 146 nodes; k comes from the whitening covariance of the basis
  rows and `extend_rank_through_ties`.
- Basis rows whitened by their own pooled within-child covariance (both
  basis children at least 8 rows, basis rows minus 2 at least 30): 29 tied
  eigenvalues, k = 30 at 46 of 51 nodes (28 or 29 at 5).
- Whitened by the node's total covariance: k = max(2, 30 minus test rows) at
  46 of 46 nodes with fewer than 30 test rows; k = 2 at 49 of 49 otherwise.
- markov/ward has the largest share of nodes with k from 28 to 30 (13 of 24;
  14 of 48, 12 of 29, 12 of 45 in the other cells), hence its median 28.
- Tree's split, raw p <= 0.001, by k: 84 of 102 edges at k 28 to 30; 2 of 80
  at k = 30 minus test rows; 38 of 110 at k = 2. Random splits at 0.01:
  2.89%, 0.06%, 0.76%.
- Edges at k 28 to 30 per cell: raw p <= 0.001 at 16 of 28, 18 of 24, 24 of
  24, 26 of 26; reached by Tree-BH 0, 4, 0, 26; significant 0, 0, 0, 26.
- Limits: random splits keep the parent's basis and k; one seed per case;
  identity noise covariance, 30 features; no real data; changes after first
  outputs are listed in the captures.

## Evidence

`raw/inbox/choir-cell-lines/cancer-pc-edge-test-20261112.md` holds the
description of the edge test, the setup, the table of the four variants, the
table of the largest nodes, the reading of the three hypotheses, the
deviations and the limits; protocol, scripts and the complete output are in
`raw/inbox/choir-cell-lines/cancer-pc-edge-test-scripts-20261112.py.txt`.
The rule for k, the eigenvalue ties and the rejections by k are in
`raw/inbox/choir-cell-lines/cancer-pc-edge-k-20261112.md` with protocol,
scripts and output in
`raw/inbox/choir-cell-lines/cancer-pc-edge-k-scripts-20261112.py.txt`.

## Links

- [[selected-tree-null-over-splitting]]
- [[cancer-pc-over-splits]]
- [[cancer-pc-statistics-diagnosis]]
- [[edge-selection-null-audit-20260601]]
- [[best-flat-cut-reference-benchmark]]
