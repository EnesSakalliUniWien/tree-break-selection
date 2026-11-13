---
title: Accepted Extra Splits of Today's Test on the cancer_pc Cases
type: source
status: draft
updated: 2026-11-12
sources:
  - raw/inbox/choir-cell-lines/cancer-pc-over-splits-20261112.md
  - raw/inbox/choir-cell-lines/cancer-pc-over-splits-scripts-20261112.py.txt
  - raw/inbox/choir-cell-lines/cancer-pc-operator-linkage-20261112.md
  - raw/inbox/choir-cell-lines/cancer-pc-operator-linkage-scripts-20261112.py.txt
tags:
  - benchmark
  - calibration
  - sibling-gate
  - edge-gate
---

# Accepted Extra Splits of Today's Test on the cancer_pc Cases

## Summary

On the 4 `cancer_pc` cases today's test finds 10 to 16 clusters on the
Markov/Ward tree because it accepts splits inside one true group: 33 of its 50
split nodes. On the default tree no edge inside a true group is significant, so
the sibling test is never run there; the default tree's 5 to 6 clusters come
from the edge gate, not from the sibling test. Splits inside one pure true
group are accepted only with operator "markov" and linkage "ward" together.

## Key Points

- Split nodes over the 4 cases: Markov/Ward 50 (17 between groups, 33 inside a
  group); default 22 (16 between, 3 inside, 3 pass-through).
- Disconnected cases, Markov/Ward: group 0 (256 cells) is cut into 5 and 8
  pieces; without the hub the groups of 64 cells are halved (27|38, 29|35,
  28|36). Connected cases: groups 4 and 5 are separated only in part, then the
  remaining mixture is cut 7 times.
- 32 of the 33 inside splits: 30 degrees of freedom, inflation factor 1.76 to
  2.70, statistic 51 to 183, raw p 0.0095 to 8.0e-24.
- Nodes with at least 40 cells, all from one true group: default 48,
  Markov/Ward 24. Smaller side share, median: 0.124 against 0.382. `rho` in 30
  features, median: 1.214 against 3.872.
- Same 64 cells of group 4: default 10|54 with `rho` 0.334 and edge p 0.99;
  Markov/Ward 29|35 with `rho` 4.51 and edge p 0.00021.
- Edge gate on those nodes: significant for at least one child in 0 of 48
  (default) and 17 of 24 (Markov/Ward, both children in all 17).
- Sibling test on those nodes: p-value present in 0 of 48 (default, all 48
  `Sibling_Divergence_Skipped`) and 17 of 24 (Markov/Ward); raw p <= 0.01 and
  `Sibling_BH_Different` in 12 of the 17, all 12 split.
- The reference (`sibling_best_flat_cut_reference`) refuses all 24
  Markov/Ward nodes and all 33 accepted inside splits; it also refuses 3 of 17
  between splits on Markov/Ward and 6 of 16 on the default tree.
- Operator against linkage (4 cells, same day): found clusters
  symmetric/average 6, 5, 6, 6 (mean ARI 0.8218); symmetric/ward 6, 12, 6, 9
  (0.6754); markov/average 6, 8, 6, 15 (0.7334); markov/ward 10, 13, 16, 15
  (0.3033).
- Accepted inside splits with all cells from one true group: 0, 0, 0, 13;
  with cells from several true groups (connected cases only): 3, 13, 14, 20.
  The hypothesis "the linkage decides" is not supported.
- Nodes of one true group, the 4 cells in the same order: smaller side share
  0.124, 0.228, 0.206, 0.382; `rho` in 30 features 1.214, 1.369, 2.947,
  3.872; edge significant for both children 0, 0, 0, 17.
- Why the edge test is significant there is in [[cancer-pc-edge-test]]:
  the tree chooses the split from the same cells.
- Limits: one seed per case; the classification rule
  and one self-check were changed after first outputs.

## Evidence

`raw/inbox/choir-cell-lines/cancer-pc-over-splits-20261112.md` holds the
setup, the per-case table, the table of nodes inside one true group, the
reading, the deviations and the limits; protocol, scripts and the complete
output with every split of every run are in
`raw/inbox/choir-cell-lines/cancer-pc-over-splits-scripts-20261112.py.txt`.
The 4 cells of operator and linkage are in
`raw/inbox/choir-cell-lines/cancer-pc-operator-linkage-20261112.md` with
protocol, scripts and output in
`raw/inbox/choir-cell-lines/cancer-pc-operator-linkage-scripts-20261112.py.txt`.

## Links

- [[selected-tree-null-over-splitting]]
- [[cancer-pc-edge-test]]
- [[best-flat-cut-reference-benchmark]]
- [[gaussian-best-flat-cut-reference]]
- [[cancer-pc-surrogate-evidence]]
- [[cancer-pc-statistics-diagnosis]]
