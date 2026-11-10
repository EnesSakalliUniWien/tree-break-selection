---
title: Exact Sibling Test Given the Edge Gate
type: source
status: draft
updated: 2026-11-08
sources:
  - raw/inbox/choir-cell-lines/exact-gate-sibling-test-20261108.md
  - raw/inbox/choir-cell-lines/exact-gate-sibling-test-scripts-20261108.py.txt
tags:
  - calibration
  - selective-inference
  - sibling-gate
  - edge-gate
---

# Exact Sibling Test Given the Edge Gate

## Summary

Conditioning on the edge test's basis half and on the direction of its test-half
contrast makes the edge gate an exact cutoff on the test-half length. The exact
conditional sibling test is less calibrated than today's: the tree's choice
inflates the test-half length 2 to 7 times, and the reference splits put the
correction for it at about 1.

## Key Points

- Exactness checked on all 935 datasets: production edge statistics scale as f^2
  (22,976/22,976), T(1) equals production (11,488/11,488), test rows match.
- Markov/Ward splits >= 40 cells on no-group data pass 75-79% (today 72%); 10-39
  cells 13-28% (today 1%); default tree sd 0.42 false partitions 50 to 57.
- Median test-half R^2/p of tested no-group splits: 2.0-2.2 (10-39 cells), 3.1-3.7
  (40-159), 4.4-7.3 (160-400); references give c_t 1.0-1.3.
- Default-tree ARI on real groups unchanged or higher; frozen rule: no arm fixes
  large splits, all keep real groups.
- The improvement of the earlier approximation came from pairing a too-low gate
  cut with today's c, estimated for a different statistic.

## Evidence

`raw/inbox/choir-cell-lines/exact-gate-sibling-test-20261108.md` holds the
conditioning, self-checks, tables and reading; protocol, deviations, scripts and
summary output are in
`raw/inbox/choir-cell-lines/exact-gate-sibling-test-scripts-20261108.py.txt`.

## Links

- [[gate-conditioned-sibling-test]]
- [[path-depth-weight-correction]]
- [[topology-correction-chain-comparison]]
