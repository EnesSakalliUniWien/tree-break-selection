---
title: Sprs, Subspaces and Equidistances
type: source
status: draft
updated: 2026-11-11
sources:
  - raw/inbox/choir-cell-lines/spr-subspaces-equidistances-20261111.md
  - raw/inbox/choir-cell-lines/spr-subspaces-equidistances-scripts-20261111.py.txt
tags:
  - real-data
  - subtree
  - topology
  - distances
---

# Sprs, Subspaces and Equidistances

## Summary

An spr's change of `rho` equals the moved clade's transfer cost in the node's
invariant directions, and an equidistant clade is an spr that changes nothing.
On the 1,285 measured splits of the Markov/Ward CHOIR tree, the sprs that
raise `rho` are related to equidistance in the tree's space (tied five times
as often) but come mostly from the change of space: the tree's space and the
invariant directions put 10.7% of the cells on different sides.

## Key Points

- Identity: with transfer cost dW = c b/(b+c) |mu_C - mu_B|^2 - c a/(a-c)
  |mu_C - mu_A|^2, the between/total ratio changes by exactly -dW/n in the
  whitened invariant directions (deviation 5.2e-17 in the self-check). The
  cell margin of 8 November is the one-cell case.
- Cells tied / on the wrong side: tree's space 3.2% / 2.4%, PCs 11.8% /
  12.0%, invariant directions 12.3% / 10.3%, fuzzy space 25.7% / 32.2%.
- Cells whose spr raises `rho`: 10.5% tied and 8.7% on the wrong side in the
  tree's space, against 2.1% and 1.4% for other cells. At least 81% of them
  are neither.
- Cells moved by greedy sprs are tied / wrong in the PCs at 23.0% / 30.8%,
  against 7.7% / 5.1% for cells that stay.
- Ward near-ties: a grandchild spr raises `rho` in 43.4% of merges with a tie
  ratio below 1.1 and in 12.8% at 3 and more; one grandchild spr makes a
  failing split pass in 23 of 4,828 cases.
- Sprs derived from equidistance in the tree's space alone (no best flat cut,
  no search) move 6,218 cells and give 61 passing splits against 57 for the
  tree split, 179 for greedy sprs and 278 for the best flat cut.
- Where tree split and best flat cut both pass, the two spaces disagree on
  1.5% of the cells.
- Limits: one dataset and tree; tie cuts 0.1 and 1.1 are conventions;
  centroid distances only; eigenvalue blocks not separated; nothing applied
  to the fit.

## Evidence

`raw/inbox/choir-cell-lines/spr-subspaces-equidistances-20261111.md` holds the
trace, tables, reading, deviation and limits; protocol, scripts, self-check
output and summary output are in
`raw/inbox/choir-cell-lines/spr-subspaces-equidistances-scripts-20261111.py.txt`.

## Links

- [[spr-alternatives-best-flat-cut]]
- [[subspace-distance-ties]]
- [[subtree-distance-topology-families]]
- [[subtree-descriptions-and-statistical-tests]]
- [[split-vs-best-flat-cut]]
- [[spr-eigenvalue-blocks]]
