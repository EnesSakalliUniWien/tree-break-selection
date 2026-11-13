---
title: Tree Splits Against the Best Flat Cut in the Rebuilt Directions
type: source
status: draft
updated: 2026-11-11
sources:
  - raw/inbox/choir-cell-lines/split-vs-best-flat-cut-20261111.md
  - raw/inbox/choir-cell-lines/split-vs-best-flat-cut-scripts-20261111.py.txt
tags:
  - real-data
  - sibling-gate
  - subtree
  - topology
---

# Tree Splits Against the Best Flat Cut in the Rebuilt Directions

## Summary

On the 48,879 CHOIR cells, every tested split was compared with the best flat
cut in the invariant directions of its own subtree. The split the tree offers
(inherited geometry) is in most subtrees unrelated to the best flat cut
(rebuilt geometry). In subtrees with several cell lines the best flat cut is
significant in 192 of 520 and the tree split in 54, and most significant best
flat cuts respect the lines. Inside single cell lines the best flat cut is
significant far more often than a Gaussian population allows.

## Key Points

- Markov/Ward, 1,475 splits with invariant directions: median ARI between the
  tree split and the best flat cut 0.036, 0.020 and 0.077 (32-159, 160-999,
  1,000 cells and more); the tree split is itself a flat cut in 47 cases.
- Median `rho` below 1,000 cells: tree split about 1, nearest flat cut 2.6 to
  3.2, best flat cut 3.3 to 4.2. The nearest flat cut moves about 21% of the
  cells.
- The 37 splits accepted by BH have a median ARI of 0.968 with the best flat
  cut.
- Reference at 0.01 (n >= 15, k <= 10, Markov/Ward): both pass in 57
  subtrees, only the best flat cut in 221, neither in 1,007.
- Two or more cell lines in the subtree (520): best flat cut significant in
  192, tree split in 54; of the 192, 148 cut at most 5% of the cells from
  their line and 80 cut no line.
- One cell line in the subtree (765): best flat cut significant in 38 of 570
  (32-159 cells) and 48 of 195 (160-999 cells); the tree split in 3.
- Default tree at 1,000 cells and more: tree split median `rho` 0.11, ARI
  -0.009 with the best flat cut.
- Limits: dependent nested subtrees of one dataset; cell lines are the only
  truth; k above 10 gives a lower bound; no fit used the best flat cut as the
  split.

## Evidence

`raw/inbox/choir-cell-lines/split-vs-best-flat-cut-20261111.md` holds the
setup, the tables and the limits; protocol, scripts and summary output are in
`raw/inbox/choir-cell-lines/split-vs-best-flat-cut-scripts-20261111.py.txt`.

## Links

- [[subtree-descriptions-and-statistical-tests]]
- [[best-flat-cut-reference-real]]
- [[gaussian-best-flat-cut-reference]]
- [[ceiling-rule-real-and-benchmark]]
- [[subtree-spectral-alignment-and-crossing-clusters]]
- [[spr-alternatives-best-flat-cut]]
