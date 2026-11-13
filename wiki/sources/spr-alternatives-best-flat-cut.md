---
title: Alternative sprs on the Tree against the Best Flat Cut
type: source
status: draft
updated: 2026-11-11
sources:
  - raw/inbox/choir-cell-lines/spr-alternatives-best-flat-cut-20261111.md
  - raw/inbox/choir-cell-lines/spr-alternatives-best-flat-cut-scripts-20261111.py.txt
  - tree_break_selection/tree/spr.py
  - tree_break_selection/tree/elements/partition.py
tags:
  - real-data
  - sibling-gate
  - subtree
  - topology
---

# Alternative sprs on the Tree against the Best Flat Cut

## Summary

An spr inside a tested node moves a clade from one child to the other and
leaves the node's cells unchanged, so the node's invariant directions, its
best flat cut and its reference do not react; the best flat cut's `rho` is
the upper bound of every spr alternative. On the Markov/Ward tree of the
48,879 CHOIR cells no split whose best flat cut is refused passes after sprs
(0 of 1,007), and greedy sprs make the tree's split pass in 122 of the 221
splits where only the best flat cut passes. Nodes are identified by their
leaf sets (`Partition`, taken over from BranchArchitect).

## Key Points

- Measured: 1,285 tested splits with 32 cells or more and 1 to 10 invariant
  directions; 190 splits with 11 to 15, including every split with 1,000
  cells or more, are outside the stored reference.
- Median `rho` by size (32-159 / 160-999 cells): tree split 0.97 / 0.92,
  after the first spr 1.60 / 1.49, after greedy sprs 3.19 / 2.61, best flat
  cut 4.19 / 3.22. Passing at 0.01: tree 15 / 42, greedy sprs 77 / 102, best
  flat cut 116 / 162.
- In the 122 rescued splits the first pass needs a median of 3 sprs
  (quartiles 1 to 8, largest 195; one spr in 38) and moves a median of 0.267
  of the cells. The median line cut share falls from 0.144 (tree split) to
  0.038 at the first pass and 0.021 at the greedy end.
- Greedy sprs end exactly at the best flat cut in 260 of 1,285 splits; in a
  simulated check with two clear groups they stop at a nearby flat cut
  (`rho` 1.78 against 9.2 to 9.9).
- Targeted sprs (fewest clades that must change sides to reach the best flat
  cut): median 16 to 20 where the tree split fails, with single cells as half
  or more of the moves; 2 where tree split and best flat cut both pass.
- Moved clades of 10 cells or more have a majority line of at least 90% in
  883 of 953 greedy and 466 of 500 targeted moves.
- Below the split the best flat cut reacts: a child's best flat cut passes in
  92 of 170 testable children before the sprs and in 25 of 180 after.
- In 4 of 1,285 splits (k = 9 or 10) the greedy sprs reached a higher `rho`
  than the table's search, so the search missed the best flat cut there; no
  verdict changed.
- Limits: one dataset and tree; per-split p-values without the BH step; sprs
  not applied to the fit; greedy counts are upper bounds; the rescued splits
  are selected by the same data that gives their p-value.

## Evidence

`raw/inbox/choir-cell-lines/spr-alternatives-best-flat-cut-20261111.md` holds
the setup, tables, reading and limits; protocol, amendment, scripts,
self-check output and summary output are in
`raw/inbox/choir-cell-lines/spr-alternatives-best-flat-cut-scripts-20261111.py.txt`.
The spr on a `PosetTree` is `tree_break_selection/tree/spr.py`; node identity
is `tree_break_selection/tree/elements/partition.py`.

## Links

- [[split-vs-best-flat-cut]]
- [[best-flat-cut-divisive]]
- [[subtree-descriptions-and-statistical-tests]]
- [[gaussian-best-flat-cut-reference]]
- [[poset-tree]]
- [[spr-subspaces-equidistances]]
- [[spr-eigenvalue-blocks]]
- [[spr-refit]]
