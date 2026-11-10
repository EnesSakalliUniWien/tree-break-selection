---
title: Alternative Separations at Tied Eigenvalues (Prototype)
type: source
status: draft
updated: 2026-11-08
sources:
  - raw/inbox/choir-cell-lines/alternative-separations-prototype-20261108.md
  - raw/inbox/choir-cell-lines/alternative-separations-prototype-scripts-20261108.py.txt
tags:
  - topology
  - spectral
  - alternative-clusters
  - prototype
---

# Alternative Separations at Tied Eigenvalues (Prototype)

## Summary

A prototype that flags subtrees whose top eigenvalues tie (multiplicity m >= 2),
cuts them along the tree into m + 1 blocks and reports every split of the blocks
the tree cannot hold. It finds the designed case (balanced groups) but also flags
ties inside single populations in subtrees the tree chose, and finds no clean tie
between cell lines in the real data.

## Key Points

- Alternatives need m >= 3: with 3 blocks a tree holds all 3 splits.
- Self-check: equal-distance groups give m = 3, 4 pure blocks, 5 of 7 splits held,
  dropped 2+2 splits at 0.998 and 0.990 of the tree's split.
- Simulated balanced data: root flagged in 10/10 per tree with pure blocks; dropped
  splits at 1.03-1.05 of the tree's split.
- Flags inside one population: 0.2-3.9 per no-group dataset, 6.9-27.4 per group
  dataset; for m >= 3 they report alternatives too.
- Real CHOIR data: at >= 500 cells only m = 2 flags (no alternatives possible); at
  >= 100 cells (exploratory) 206 and 308 flags, none with pure blocks from
  different lines.
- Needed before use below the root: a noise band for selected subtrees that
  includes the tree's choice.

## Evidence

`raw/inbox/choir-cell-lines/alternative-separations-prototype-20261108.md` holds
the method, self-check and tables; code and output are in
`raw/inbox/choir-cell-lines/alternative-separations-prototype-scripts-20261108.py.txt`.

## Links

- [[subtree-spectrum-topology]]
- [[split-count-correction]]
