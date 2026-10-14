---
title: Distribution Equation Fitting Comparison
type: source
status: reviewed
updated: 2026-09-28
sources:
  - raw/inbox/wolfram-distribution-fitting-20260928.md
  - reports/distribution_equation_comparison_20260928/README.md
  - reports/distribution_equation_comparison_20260928/compare_distributions.py
  - reports/distribution_equation_comparison_20260928/heldout_comparison.csv
  - reports/distribution_equation_comparison_20260928/rank_mixture_summary.csv
  - reports/distribution_equation_comparison_20260928/manifest.json
  - benchmarks/diagnostics/calibration/selected/hierarchy/selected_hierarchy_geometry_covariates.py
  - benchmarks/diagnostics/calibration/statistics/statistic_distribution_shape_panel.py
  - tree_break_selection/hierarchy_analysis/statistics/projection/selected_gaussian_radial.py
tags:
  - calibration
  - degrees-of-freedom
  - projection
---

# Distribution Equation Fitting Comparison

## Summary

Recovered methods include six log-ratio regression equations, empirical selected
ratio tails, Satterthwaite scale/effective-df fitting and conditional chi-radius
integration. A new independently trained/evaluated Gaussian diagnostic compares
nominal, fitted and exact example laws. Wolfram Alpha and Wolfram Language check
the moment equation and the controlled adaptive-energy conditional law.

## Key Points

- Moment fitting of a times chi-square nu gives a = variance/(2 mean) and
  nu = 2 mean squared/variance. Matching moments does not establish tail validity.
- For two independent standard normal coordinates with fixed projection basis,
  energy fraction q selects rank one when U = Z1 squared / total energy >= q.
  Conditional on that event, the exact survival is the average of
  exp[-t/(2 cos-squared(theta))] over theta from zero to arccos(sqrt(q)).
  Rank two retains chi-square two because its selection is purely directional.
- At alpha .01, energy .9 over both ranks gives nominal rejection 1.4850%,
  fitted scale/df 1.0300%, and exact conditional 1.0135% on 200,000 independent
  evaluation vectors. Conditional rank-one results use only 40,963 vectors.
- Conditioning chi-square two on exceeding four changes support without changing
  rank. Nominal rejection is 7.4818%, moment-fit rejection 2.4512%, and exact
  conditional rejection 1.0369% among 27,293 selected evaluation vectors.
- Fitting can approximate the controlled energy law closely while missing other
  tails. At alpha .001, the radial example's frozen moment fit has population
  rejection 7.364 times alpha; maximum likelihood is worse at 12.365 times alpha.
- The historical source-family holdout gives edge_spectral_modes AUC .966050
  but log-ratio R-squared .062485. A ratio predictor is not a calibrated CDF.
  Seven of 79 historical tail contexts were marked admissible under the old
  support contract; no current production promotion follows.
- Same-data learned PCA, covariance estimation, inferred trees and empirical
  sibling inflation are outside the exact-law derivation. Full pipeline
  selected-tail and clustering calibration remain open.

## Evidence

Training uses 100,000 vectors and evaluation uses 200,000 separate vectors,
with matched draws across examples and four alpha values. The runner checks
256 production kernel results, four existing radial-tail results, ten Wolfram
checkpoints and independent quadrature. Source hashes are unchanged. Baseline
reference/kernel and radial tests passed 8 and 15 tests respectively. The
report includes exact equations, counts, pointwise intervals, fitted parameters,
curves and limitations. No production behavior changed.

## Links

- [[lower-dimension-degrees-of-freedom-audit]]
- [[selected-tail-equation-cloud-run-20260604]]
- [[selected-ratio-tail-law-diagnostic-20260602]]
- [[statistic-distribution-shape-panel-20260613]]
- [[independent-projection-vandermonde-control]]
- [[open-mathematical-questions]]
