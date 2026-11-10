---
title: Large-Split Correction with Tree Depth, Weight and Whole Paths
type: source
status: draft
updated: 2026-11-08
sources:
  - raw/inbox/choir-cell-lines/path-depth-weight-correction-20261108.md
  - raw/inbox/choir-cell-lines/path-depth-weight-correction-scripts-20261108.py.txt
tags:
  - calibration
  - selective-inference
  - sibling-gate
---

# Large-Split Correction with Tree Depth, Weight and Whole Paths

## Summary

Adding tree depth, edge-blocked-only reference weight and whole paths (ancestors
and descendants) to today's empirical-null inflation does not fix large splits.
With all three, Markov/Ward tested splits of 40 or more cells on no-group data
still pass at 50% instead of about 1%, and false partitions stay at 291 to 300
of 300. Real groups on the default tree are unchanged.

## Key Points

- Candidates swap only the factor of each tested split; production code runs end
  to end. Replay of today's factor matched 11,488/11,488 tested splits.
- Markov/Ward, sd 0.25, splits >= 40 cells: pass share 71.8% (today), 71.9%
  (depth), 66.4% (weight), 66.7% (paths), 49.9% (all); median factor 2.5 to 3.4.
  For 40-159 cells "all" reaches 6.1 against 7.0 needed; 160-400 cells stay at
  2.9 against 7.6.
- Default tree, sd 0.42: depth and all lower false partitions from 50 to 37 of
  300; sd 0.25 stays at 1 to 4.
- Default-tree mean ARI on four signal models is unchanged (within 0.001);
  Markov/Ward ARI rises with "all" (balanced 0.315 to 0.445) because it splits
  less noise.
- Decision rule fixed before results: no candidate fixes large splits; all keep
  real groups.
- Mechanism (descriptive endpoint): at the same size, tested no-difference splits
  score 1.3 to 7 times more than reference splits. The edge gate sorts splits by
  nearly the same difference the sibling score measures, so every reference split
  sits on the low-score side; reweighting them cannot reach the needed factor.

## Evidence

`raw/inbox/choir-cell-lines/path-depth-weight-correction-20261108.md` holds the
protocol summary, all endpoint tables and the deviation; the protocol, deviations
and scripts are in
`raw/inbox/choir-cell-lines/path-depth-weight-correction-scripts-20261108.py.txt`.
The needed-factor definition and the no-group data follow
[[topology-correction-chain-comparison]].

## Links

- [[topology-correction-chain-comparison]]
- [[sibling-inflation-debugging]]
- [[descendant-chain-family-test]]
- [[top-down-traversal]]
