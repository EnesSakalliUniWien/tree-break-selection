---
title: Distances in Different Subspaces and Equidistant Cells
type: source
status: draft
updated: 2026-11-08
sources:
  - raw/inbox/choir-cell-lines/subspace-distance-ties-20261108.md
  - raw/inbox/choir-cell-lines/subspace-distance-ties-scripts-20261108.py.txt
tags:
  - topology
  - spectral
  - distances
---

# Distances in Different Subspaces and Equidistant Cells

## Summary

For each split, every cell's distance to the two children's centroids was compared
in the tree's own space, the test PCs, the test metric, and the invariant and fuzzy
subspaces. Real separations leave no cell equidistant; splits inside one population
leave 7-21% of cells equidistant, and the tree's space and the test space put 7-12%
of cells on different sides. In the real data, cell lines are cut where their cells
are equidistant in the test space.

## Key Points

- Between-group splits: 0 tied and 0 wrong-side cells in the tree, test and invariant
  spaces; fuzzy directions agree with the invariant side only about half the time.
- Within-group splits: tied cells 7-11% (tree space), 11-21% (test space); the side
  differs between the tree and test spaces for 7-12% of cells; separation 0.5-1.4 in
  the tree's space against 0.2-0.4 in the test space.
- Real CHOIR data: in splits that cut lines, cut-line cells are tied in the test space
  at 23% (Markov/Ward) and 40% (default) against 7% for other cells.
- The tie cut |m| < 0.1 is a stated convention, not a test.

## Evidence

`raw/inbox/choir-cell-lines/subspace-distance-ties-20261108.md` holds the measure,
sanity check and tables; protocol, code and output are in
`raw/inbox/choir-cell-lines/subspace-distance-ties-scripts-20261108.py.txt`.

## Links

- [[subtree-spectrum-topology]]
- [[alternative-separations-prototype]]
- [[topology-correction-chain-comparison]]
- [[spr-subspaces-equidistances]]
