---
title: Cancer PC Statistics Diagnosis
type: source
status: reviewed
updated: 2026-11-04
sources:
  - reports/cancer_pc_statistics_debug_20261104/README.md
  - reports/cancer_pc_statistics_debug_20261104/run_diagnostic.py
  - reports/cancer_pc_statistics_debug_20261104/summary.csv
  - reports/cancer_pc_statistics_debug_20261104/root_factors.csv
  - reports/cancer_pc_statistics_debug_20261104/calibration_factors.csv
  - reports/cancer_pc_statistics_debug_20261104/frozen_tree_swaps.csv
  - reports/cancer_pc_statistics_debug_20261104/fixed_pair_reference_check.json
  - reports/cancer_pc_statistics_debug_20261104/independent_verification.json
  - reports/cancer_pc_statistics_debug_20261104/manifest.json
  - reports/cancer_pc_statistics_debug_20261104/statistics_diagnosis.png
tags:
  - statistics
  - calibration
  - diffusion
  - benchmarks
---

# Cancer PC Statistics Diagnosis

## Summary

Eight native runs over the four cancer PC cases and two previously examined
seeds isolate excessive raw same-mean sibling rejection, covariance-related
singleton stopping, empirical-null contamination, and a mean-contrast target
that can miss heterogeneous descendants. No production behavior was changed.

## Key Points

- In the hub/disconnected seed 20261106, 422 of 483 zero-generating-mean
  sibling records have raw p below 0.01. At 183 two-leaf nodes, fitted diagonal
  covariance and fitted null-whitened PCA yield Q times branch multiplier
  approximately 30, but the reference uses one degree of freedom. The raw
  threshold is crossed by 182 pairs; none becomes a native sibling rejection.
- A separate predetermined Gaussian-pair control with known noise isolates
  the estimated-direction problem: projecting onto the tested contrast retains
  a chi-square(30) norm, despite nominal rank one. Using chi-square(1) rejects
  10,000/10,000 simulated null pairs at alpha 0.01, versus 102/10,000 with a
  fixed coordinate and 112/10,000 using the full-space reference. This is not
  the full selected-tree pipeline's error rate.
- The singleton root at seed 20261104 stays closed when branch scaling is
  removed, full projection is used, or common identity covariance is supplied.
  Exact known heterogeneous row-noise covariance opens both root edges. That
  oracle information is unavailable from one observation alone and does not
  provide a calibrated selected test or complete clustering repair.
- At the 16/480 failed root, 11 mean-different records among 494 supported
  records supply 0.708302 of the weighted-statistic numerator. Removing them
  while retaining remaining family weights changes c from 14.843856 to 4.867973
  and root p from 0.876715 to 0.000213350. There are only two stopping-event
  dependency groups; all eight root calibrations fail internal support
  thresholds if enforcement is enabled.
- Frozen topology-plus-branch exchanges preserve one versus six clusters on
  either matched dataset in both disconnected seeds. This establishes the
  fitted-tree contribution without separating topology from branch lengths.
- A connected-case 384-leaf stopping node contains three populations, but its
  1/383 sibling split has exactly equal generating means. Testing that mean
  contrast cannot detect the deeper symmetric mixture. A missing split need
  not be an error under the hypothesis actually tested.

## Evidence

The report states the contrast-weight variance identity, two-leaf algebra,
reference-law control, one-factor interventions and native traversal results.
Source and input/metric parity, 3,960 raw sibling replays, 16 root edge replays,
eight full-trace partition replays and four matched frozen-tree replays passed.
Independent checks recalculate 120 factor p-values, verify 241 source hashes,
and reconstruct the mixed stopping node's population counts. Generator and
calibration tests pass 59 cases. The four-panel figure was visually reviewed.

Across the eight runs, 50 of 116 eligible edges and 24 of 25 sibling frontier
tests reject; none of these rejected contrasts has zero generating mean.
These are realized diagnostic counts, not a selected Type-I-error/FDR estimate.
The selected cases, within-tree dependence, original cancer covariance,
count preprocessing, descendant selective inference and full-scale NNLS remain
outside the established evidence. Production sources retain their pre-existing
working-tree modifications, including the empirical-null family fallback.

## Links

- [[cancer-pc-surrogate-evidence]]
- [[analytic-selected-root-laws]]
- [[conditional-diffusion-monte-carlo-evidence]]
- [[choir-cell-lines-recovery-plan]]
