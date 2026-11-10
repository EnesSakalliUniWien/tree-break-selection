---
title: Rule-Out Test on the Real CHOIR Cell-Line Data
type: source
status: draft
updated: 2026-11-08
sources:
  - raw/inbox/choir-cell-lines/split-count-exclusion-real-20261108.md
  - raw/inbox/choir-cell-lines/split-count-exclusion-real-scripts-20261108.py.txt
tags:
  - calibration
  - real-data
  - choir
  - tree-choice
---

# Rule-Out Test on the Real CHOIR Cell-Line Data

## Summary

End to end on all 48,879 cells, the rule-out test (count all possible splits; a split
is real only if it lies in its invariant directions and leaves at most 5% of cells
equidistant) roughly doubles the cell lines Markov/Ward recovers (41 to 80) and cuts
its clusters from 954 to 664. The default tree is unchanged because its giant cluster
comes from the topology.

## Key Points

- Replay of today's test reproduces the 2026-11-06 runs exactly (75 / 954 clusters).
- Markov/Ward: ARI 0.549 to 0.765, NMI 0.863 to 0.901, lines 90/90 41 to 80, lines with
  >= 90% of cells in one cluster 56 to 91; still 664 clusters for 190 lines.
- Default tree: 69-82 clusters and 27-29 lines in every arm; giant cluster 28,611-32,680
  cells.
- Counting all splits without the rule-out is worse than today (1,605 clusters, 24
  lines): small splits pass.
- No no-group parts exist in these data, so the remaining within-line splits cannot be
  called false or real here.

## Evidence

`raw/inbox/choir-cell-lines/split-count-exclusion-real-20261108.md` holds the setup,
replay checks, tables and reading; code and metric lines are in
`raw/inbox/choir-cell-lines/split-count-exclusion-real-scripts-20261108.py.txt`.

## Links

- [[split-count-exclusion]]
- [[choir-cell-lines-recovery-plan]]
- [[subspace-distance-ties]]
