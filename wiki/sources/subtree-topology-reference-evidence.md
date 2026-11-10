---
title: Subtree Distances and Gaussian Topology Reference Evidence
type: source
status: reviewed
updated: 2026-11-09
sources:
  - raw/inbox/choir-cell-lines/subtree-topology-reference-20261109.md
  - raw/inbox/choir-cell-lines/subtree-topology-reference-20261109-evidence.json
  - raw/inbox/choir-cell-lines/subtree-topology-reference-20261109-sources.py.txt
  - raw/inbox/choir-cell-lines/subtree-topology-mathematics-20261109.md
  - raw/inbox/choir-cell-lines/subtree-topology-mathematics-20261109-wolfram.json
  - raw/inbox/choir-cell-lines/subtree-topology-geometry-20261109.md
  - raw/inbox/choir-cell-lines/subtree-topology-calibration-20261109.md
  - raw/inbox/choir-cell-lines/subtree-topology-bootstrap-20261109.md
tags:
  - subtree
  - topology
  - distances
  - selective-inference
---

# Subtree Distances and Gaussian Topology Reference Evidence

## Summary

Distances constrain a family of subtree topologies after specifying the
reconstruction rule and treatment of equidistances. The diagnostic enumerates
Ward alternatives, fits distances conditional on candidate topologies, and
implements common-Gaussian conditional references for individual cuts, paired
topology scores and searches. The mathematical reference is conditional on its
assumptions; broad native-score power and final-partition control remain open.

## Key Points

- Rectangle axis views give crossing trees, with union two and intersection zero.
  A square gives two Ward trees; an equidistant four-leaf simplex gives 15 binary
  refinements. All simplex additive fits have zero internal lengths, so the
  distance model is an unresolved star.
- Full rectangle additive distance RMSE is .1927 for AB|CD versus 1.0716 for
  AC|BD. Equal rooted topology can still have unequal branch lengths. Additive
  pair distances do not identify the root.
- Fixed-label, Gram-only rank has the stated Beta law. The conditional rotation
  rank allows arbitrary deterministic recomputed scores and selected-subtree
  events under the stronger common-mean Gaussian union null.
- Selected-winner nulls reject 219/1,200 with the naive chain and 60/1,200 after
  conditioning at nominal .05; the exact selected reference rejects 59/1,200.
  The all-label maximum rejects 13/320 and median-selected descendant 8/160.
- The 138 integrated attempts include 0/50 cut-search, 1/40 topology-search and
  2/40 fixed-pair null rejections. All 616 candidate-adjusted tails were verified.
  Intervals are broad and do not certify calibration at .05.
- Both native signal roots recover their planted split but fail to reject at
  99 and 999 transitions. Doubling the shift in an exploratory paired follow-up
  gives p=.034 and .017. This is limited power evidence.
- In 12,000 bootstrap/permutation replicates, the square's two fixed-geometry
  Ward alternatives expand to seven under feature resampling and 14 under
  larger measurement noise. These are stability frequencies, not null p-values.

## Evidence

The combined report and JSON retain protocols, source hashes, Wolfram requests,
geometry fits, finite-panel counts, exact intervals and integration verification.
The source archive includes the local runnable scripts. Full chains and figures
remain in `results/analyses/subtree_topology_reference_20261109/`.
Production source and defaults were preserved; no installations or commits.

## Links

- [[subtree-distance-topology-families]]
- [[gaussian-root-calibration-and-cut-geometry]]
- [[subtree-spectrum-topology]]
- [[alternative-separations-prototype]]
- [[root-validation-generalisation-evidence]]
