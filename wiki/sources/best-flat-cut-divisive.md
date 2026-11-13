---
title: Divisive Fit at the Best Flat Cut on the CHOIR Cell Lines
type: source
status: draft
updated: 2026-11-11
sources:
  - raw/inbox/choir-cell-lines/best-flat-cut-divisive-20261111.md
  - raw/inbox/choir-cell-lines/best-flat-cut-divisive-scripts-20261111.py.txt
tags:
  - real-data
  - sibling-gate
  - subtree
  - topology
---

# Divisive Fit at the Best Flat Cut on the CHOIR Cell Lines

## Summary

A trial fit that starts with all 48,879 cells and splits each set at its best
flat cut in the set's own invariant directions, while the Gaussian reference
rejects one population at 0.01, gives 323 clusters with ARI 0.874 against the
cell lines. It separates the lines from each other almost completely and
splits 88 lines internally. It is not Tree-Break Selection: it uses no
diffusion distances, no tree, no edge gate, no pass-through and no BH step.

## Key Points

- ARI 0.874, NMI 0.941, ARI against CHOIR_1 0.617; ceiling 3 on Markov/Ward
  has 0.833, 0.940 and 0.538.
- Lines merged only 8 and both 8 (ceiling 3: 33 and 22); lines split only 81
  (ceiling 3: 16); lines recovered at 90/90: 93 (ceiling 3: 119).
- Every line is the majority of at least one cluster; under ceiling 3, 22 are
  the majority of none. Joining clusters with the same majority line gives ARI
  0.957 against 0.841.
- 322 accepted cuts, 169 with more than 10 invariant directions; their
  p-values come from run-time nulls searched with the same search, which is a
  valid reference for the searched value at fixed (n, k).
- 44 clusters have fewer than 20 cells; 9 accepted cuts isolate fewer than 5
  cells.
- Limits: per-test level 0.01 with no family-wise or FDR step; k chosen from
  the same cells; the sides of an accepted cut are not Gaussian samples; one
  dataset; no benchmark run; no complete-fit null calibration.

## Evidence

`raw/inbox/choir-cell-lines/best-flat-cut-divisive-20261111.md` holds the
procedure, the tables and the limits; protocol, scripts, the run-time null
values and the summary output are in
`raw/inbox/choir-cell-lines/best-flat-cut-divisive-scripts-20261111.py.txt`.

## Links

- [[split-vs-best-flat-cut]]
- [[subtree-descriptions-and-statistical-tests]]
- [[best-flat-cut-reference-real]]
- [[gaussian-best-flat-cut-reference]]
- [[ceiling-rule-real-and-benchmark]]
- [[spr-alternatives-best-flat-cut]]
