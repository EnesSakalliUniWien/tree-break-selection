---
title: Branch Distribution and Variance Benchmark
type: source
status: reviewed
updated: 2026-09-28
sources:
  - reports/branch_distribution_variance_20260928/README.md
  - reports/branch_distribution_variance_20260928/protocol.md
  - reports/branch_distribution_variance_20260928/compatibility.md
  - reports/branch_distribution_variance_20260928/evaluation/verification.json
  - reports/branch_distribution_variance_20260928/evaluation/summary.csv
  - tree_break_selection/hierarchy_analysis/statistics/projection/spectral/tree_estimator.py
  - benchmarks/shared/runners/tbs_diffusion_runner.py
tags:
  - calibration
  - diffusion
  - nnls
---

# Branch Distribution and Variance Benchmark

## Summary

The missing internal-distribution and branch-variance comparison completed
48,600 runs on the same 900 Gaussian dataset pairs as the earlier independent
diffusion benchmark. It tested current settings and a report-local discovery-only
internal-state adapter without changing production code.

## Key Points

- Current NNLS at sibling alpha .01 and edge alpha .001: empirical barycenters
  with edge variance produce 25/250 high-dimensional null false partitions;
  leaf-only and branch-length states produce 0/250 there. Branch-length states
  produce 4/250 identity2 false partitions, versus 6/250 leaf-only.
- Independent NNLS branch-length states change signal ARI from .25120 to
  .25707 in ar1_8 and from .15531 to .14621 in identity40, with 1/250 identity2
  false partitions. This is not a uniform recovery improvement.
- NNLS now affects independent tests through branch-length-dependent discovery
  PCA: 6,887/27,900 node p-values and 1,988 ranks change, but only 3/900 final
  partitions differ from diffusion with the same internal-state mode.
- Internal rows do not add independent observations. Discovery-only internal
  PCA preserves the independent Gaussian contract; same-data MP and estimated
  covariance are a different construction.
- The edge-variance flag does not disable sibling branch variance. Toggling it
  preserves raw sibling statistics but changes support, inflation and gate
  availability. A geometric multiplier cannot be added to the independent
  known covariance without correcting the reference scale or specifying a
  different stochastic covariance model.
- The adaptive diffusion dispatcher does not forward internal-state flags.
  The benchmark uses the existing distance runner and verifies actual settings;
  the dispatcher limitation and lack of a public independent internal-state
  option remain unresolved.

## Evidence

All 502,200 node records and 5,400 independent contexts were checked; 7,200
previous baseline rows reproduce. All 49 focused tests pass. No captured
warnings or evaluation errors occurred. The initial pilot export error is
preserved and superseded by a verified 486-run pilot. At root-tail cutoffs
.001/.005/.01/.05, none of the 72 independent cells rejects excess tails after
Holm correction; the finite sample does not prove calibration. Matched datasets,
alpha sensitivity, unavailable tests and dependent node counts remain explicit.

## Links

- [[independent-diffusion-nnls-benchmark]]
- [[independent-gaussian-production-benchmark]]
- [[diffusion-nnls-calibration-open-points]]
