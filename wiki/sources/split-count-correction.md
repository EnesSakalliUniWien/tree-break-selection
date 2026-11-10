---
title: Correction for the Tree's Choice from Split Counts
type: source
status: draft
updated: 2026-11-08
sources:
  - raw/inbox/choir-cell-lines/split-count-correction-20261108.md
  - raw/inbox/choir-cell-lines/split-count-correction-scripts-20261108.py.txt
tags:
  - calibration
  - selective-inference
  - sibling-gate
  - tree-choice
---

# Correction for the Tree's Choice from Split Counts

## Summary

A union bound over every split of a subtree's known leaf set (2^(n-1) - 1 splits)
removes almost all false passes of the largest no-group splits and cuts default-tree
false partitions, without losing real groups, but it lets small splits through.
Restricting the count to the subtree's invariant directions changes little because
noise splits do not follow those directions.

## Key Points

- Markov/Ward no-group splits of 160-400 cells pass 94% today and 0.2-0.4% with the
  all-splits count; 40-159 cells 61-63% to 16-21%; 10-39 cells 1-3% to 85-93%.
- Default tree at sd 0.42: large splits 55% to 1.2%; false partitions 50 to 12 of
  300; real-group ARI unchanged.
- Production's per-subtree Marchenko-Pastur count is always 0 here (each subtree is
  whitened by its own covariance first); invariant directions were read from the
  subtree covariance's eigenvalue bulk instead (no-group root 0, four-group root 3).
- Noise splits are separable in the invariant directions in only 0-9% of cases.
- Large noise splits score about (2/pi) n + df (ratio 1.1-1.3 at 160-400 cells);
  small ones about four times that, so their fixed-split chi-square law already
  fails and no count repairs it.
- Frozen rule: no arm fixes large splits (11-16% at >= 40 cells); both keep real
  groups. A rule combining the count (large) with today's factor (small) was not run.

## Evidence

`raw/inbox/choir-cell-lines/split-count-correction-20261108.md` holds the arms, the
amendment, tables and reading; protocol, deviations, scripts and output are in
`raw/inbox/choir-cell-lines/split-count-correction-scripts-20261108.py.txt`.
Literature context: Kimes et al. 2017 (doi:10.1111/biom.12647) and Huang et al. 2015
(doi:10.1080/10618600.2014.948179), retrieved via PubMed; Cover 1965 from secondary
sources only.

## Links

- [[exact-gate-sibling-test]]
- [[path-depth-weight-correction]]
- [[topology-correction-chain-comparison]]
