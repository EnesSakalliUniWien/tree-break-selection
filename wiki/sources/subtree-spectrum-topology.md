---
title: Subtree Spectrum, Multiplicity and Separations Outside the Topology
type: source
status: draft
updated: 2026-11-08
sources:
  - raw/inbox/choir-cell-lines/subtree-spectrum-topology-20261108.md
  - raw/inbox/choir-cell-lines/subtree-spectrum-topology-scripts-20261108.py.txt
tags:
  - topology
  - spectral
  - alternative-clusters
---

# Subtree Spectrum, Multiplicity and Separations Outside the Topology

## Summary

Per subtree: leaves, number of splits, eigenvectors above the noise band,
eigenvalues and their geometric multiplicity, and which group separations the
tree cannot hold. A rooted binary tree holds 2g - 3 of the 2^(g-1) - 1
separations of g groups; when the top eigenvalues form one eigenspace
(balanced groups, multiplicity 3) the dropped separations are almost as strong
as the kept one, so they are alternative clusterings the tree cannot show.

## Key Points

- Multiplicity rule (block equality test at 0.01) recovers multiplicity 3 in
  99% of runs; of three distinct eigenvalues (30, 20, 10) the top two merge into
  one eigenspace in 78% of runs at 100 cells; all are resolved at 1,000 (99.8%).
- No-group roots: no direction above the noise band at log-depth sd 0.25; the
  depth direction appears at sd 0.42 (60%) and 1.0 (100%).
- Balanced 4 groups: k = 3, multiplicity 3; the tree holds every 1-vs-3 split
  and one of three 2+2 splits; the other two are 0.95-0.97 as strong.
- Imbalanced 4 groups: multiplicity 1+1+1; dropped separations are near copies
  of kept ones, 8% weaker.
- Density-variation data: the strongest separation (a 600-cell group alone) is
  not a clade in either tree.
- Markov/Ward subtrees of no-group data show 3-7 directions above their noise
  band: subtrees the tree chose are already shaped by that choice.
- Real CHOIR data: default-tree nodes cut lines in 49% of analysed nodes (423
  lines), Markov/Ward in 27% (70); the link to multiplicity is weak.

## Evidence

`raw/inbox/choir-cell-lines/subtree-spectrum-topology-20261108.md` holds the
definitions, calibration and tables; protocol, scripts and output are in
`raw/inbox/choir-cell-lines/subtree-spectrum-topology-scripts-20261108.py.txt`.
Split compatibility and split networks: Huson and Bryant 2006
(doi:10.1093/molbev/msj030, via PubMed).

## Links

- [[split-count-correction]]
- [[topology-correction-chain-comparison]]
- [[choir-cell-lines-recovery-plan]]
