---
title: Topology, Estimated Correction and Chain Test on No-Group Data
type: source
status: draft
updated: 2026-11-07
sources:
  - raw/inbox/choir-cell-lines/topology-correction-chain-20261107.md
  - raw/inbox/choir-cell-lines/topology-correction-chain-scripts-20261107.py.txt
tags:
  - topology
  - calibration
  - selective-inference
  - sequencing-depth
---

# Topology, Estimated Correction and Chain Test on No-Group Data

## Summary

On 900 count datasets with one true group, the default tree splits off a few
cells at the top when depth spread is small, and Ward divides pure noise into
balanced halves. Today's empirical-null inflation corrects small splits about
right but large splits about three times too little. The root chain test never
rejected (0/900), but on 400-cell data with four true groups it rejected only
4/30 tops.

## Key Points

- Default tree, log-depth sd 0.25 / 0.42 / 1.0: median root smaller side 2 /
  5.5 / 125.5 cells; false partitions 2 / 50 / 259 of 300. Markov/Ward: 138.5 to
  161 cells, 296 to 300 of 300, median 8 clusters for one population.
- Default-tree protection at low depth is the child-parent edge gate (opens at
  the root in 1% of datasets), not the sibling correction.
- Markov/Ward root score is 8.7 to 17.6 times the one-population value; needed
  correction 7.7 to 22.1, applied 2.2 to 2.6; 95% to 99% of roots split.
- By parent size (sd 0.25): applied 9.2 vs needed 9.4 for 10-39 cells; 2.5 to
  2.6 vs 7.0 to 7.6 for 40-400 cells.
- Chain test (root of the default tree only, plug-in sigma): 0/300 rejections at
  each depth where production splits 0, 20 and 254 roots; median 2 to 7 accepted
  moves of 199; 9 to 17 s at 400 cells; a 9-step run at 1,000 cells did not
  finish in 10 minutes.
- Added power check: on 400-cell subsets with four groups, production splits
  every root (ARI 0.99 to 1.00) and the chain test rejects 4/30 (median p 0.02).

## Evidence

`raw/inbox/choir-cell-lines/topology-correction-chain-20261107.md` holds the
tables, root-gate shares, the deviations and the chain settings; the protocol
and all scripts are in
`raw/inbox/choir-cell-lines/topology-correction-chain-scripts-20261107.py.txt`.
Real-data topologies are in [[choir-cell-lines-recovery-plan]] (2026-11-06
comparison: 42 vs 134 clean lines).

## Links

- [[descendant-chain-family-test]]
- [[sibling-inflation-debugging]]
- [[sequencing-depth-normalization-literature]]
- [[top-down-traversal]]
