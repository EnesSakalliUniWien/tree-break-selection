---
title: Sibling Test Conditioned on the Edge Gate
type: source
status: draft
updated: 2026-11-08
sources:
  - raw/inbox/choir-cell-lines/gate-conditioned-sibling-test-20261108.md
  - raw/inbox/choir-cell-lines/gate-conditioned-sibling-test-scripts-20261108.py.txt
tags:
  - calibration
  - selective-inference
  - sibling-gate
  - edge-gate
---

# Sibling Test Conditioned on the Edge Gate

## Summary

The analysis treated the edge gate of a split as a cutoff on that split's own
sibling score and conditioned the sibling p-value on the gate having opened. A
correction found after the run: the edge test uses only a test half of the cells,
so that cutoff was placed too low (0.47 instead of 0.72 of the observed score,
medians) and the arms are approximations. With today's empirical-null inflation
the approximation lowers false splits but does not fix large splits; most of the
excess comes from the tree's choice of the split.

## Key Points

- Correction: `split_edge_test_rows` gives the edge test only a test half of
  each child's cells (the other half fits the edge basis), so edge statistics are
  not multiples of the all-cell sibling statistic. Holding the basis half fixed,
  the gate cuts the sibling score at a median 0.72 of its observed value; the
  arms used s* = 0.47. An exact version must condition on the basis half and on
  the test half's direction.
- The arms used p = sf(T/c) / sf(s* T/c) on chi2(df), with s* the factor on the
  edge statistics at which the Tree-BH family flips. Copy of Tree-BH matched
  production in 5,300/5,300 fits.
- Markov/Ward, splits >= 40 cells on no-group data: 72-73% pass today, 63%
  with the gate condition and today's c, 88-90% with the gate condition and
  c = 1. Default tree at sd 0.42: 55% to 33%; false partitions 2 to 0 (sd 0.25)
  and 50 to 34 (sd 0.42).
- Default-tree ARI on four real-group models is unchanged; on real splits the
  gate sits far below the score (s* about 1e-4).
- A truncated-likelihood re-estimate of c moves it by under 0.5%; the
  references that drive c are mostly edge blocked and carry no own-gate cutoff.
- On no-group splits the gate opens at about half the observed score for large
  splits (median s* 0.41-0.64) while those splits score 5 to 14 times the
  one-population value.
- Decision rule fixed before results: no arm fixes large splits; all keep real
  groups.

## Evidence

`raw/inbox/choir-cell-lines/gate-conditioned-sibling-test-20261108.md` holds the
derivation checks, all endpoint tables and the reading; the protocol and scripts
are in
`raw/inbox/choir-cell-lines/gate-conditioned-sibling-test-scripts-20261108.py.txt`.
Data and needed-factor definitions follow
[[topology-correction-chain-comparison]].

## Links

- [[path-depth-weight-correction]]
- [[topology-correction-chain-comparison]]
- [[descendant-chain-family-test]]
- [[sibling-inflation-debugging]]
