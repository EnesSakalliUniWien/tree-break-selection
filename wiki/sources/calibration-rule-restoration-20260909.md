---
title: Calibration Rule Restoration 2026-09-09
type: source
status: reviewed
updated: 2026-09-09
sources:
  - reports/calibration_restoration_20260909/README.md
  - reports/calibration_restoration_20260909/summary.json
  - reports/calibration_restoration_20260909/verification.json
  - reports/calibration_restoration_20260909/validation_protocol.md
  - reports/calibration_restoration_20260909/external_contract_check.json
  - raw/inbox/calibration-restoration-reference-notes-20260909.md
  - tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflation_correction/empirical_null_inflation_estimation.py
  - tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflated_projected_wald_annotation/pipeline.py
  - tests/statistics/35_test_empirical_null_inflation_estimation.py
  - tests/statistics/50_test_sibling_skip_annotation.py
tags:
  - calibration
  - diffusion
  - nnls
  - validation
---

# Calibration Rule Restoration 2026-09-09

## Summary

The user requested restoration of the historical empirical calibration rule
and identification of the statistical validation needed next. The restoration
returns the context-weighted plug-in chi-square sibling p-values and uses them
in traversal-aligned BH when internal support exists. Selected-tail validity
remains unproved.

## Key Points

- The verified 122-case pydiffmap NNLS rerun returns 97 results, 21 skips and
  four no-internal-support outcomes, with zero execution errors. All 97 jointly
  successful August 1/current cases have identical ARI and K; their mean ARI is
  0.797281 in both captures. The restored aggregate is not an improvement over
  August 1's 0.773363 aggregate, which includes three additional zero-ARI rows.
- Restoration recovers 92 previously unsupported results. All 97 saved label
  files were independently rescored; all 399 source hashes match. The 88
  available earlier-run data/truth signatures match; 34 were not captured.

- The local scale estimator, feature-family matching, Gaussian context weights
  and nondeflation floor follow the pre-August-11 rule. The production bypass
  that withheld empirical p-values is removed.
- Dependency-group support diagnostics, optional strict support thresholds and
  the dedicated exact independent common-scale F API remain. Tuple-valued group
  IDs now work in scale-sensitivity diagnostics. The empirical entry point
  rejects exact-F models so their ownership requirements cannot be bypassed.
- Calibration sweeps have new artifact versions, preventing reuse of cached
  withheld-p-value outcomes. Missing internal support still closes sibling
  gates and yields an unsupported benchmark result; no neutral scale or
  fixed-coordinate replacement is inserted.
- The next necessary experiment uses independent held-out null and signal
  datasets and reruns diffusion, tree construction, NNLS, covariance/PCA/rank,
  edge selection, calibration fitting, sibling BH and traversal on each one.
- P-value validity requires superuniformity under the declared null/selection
  context, not exact uniformity. A rejected uniformity test alone does not
  establish inflated false-split error or justify disabling the empirical rule.
- The existing 500-replicate external diagnostic has 27 strata, none admissible,
  and at most 470 matching independent simulations (499 required). Its 19 strata
  with at least 100 matching simulations report no scalar-tail rejection at
  alpha 0.01. This is a different diagnostic, not validation of restored diffusion
  NNLS, and its conservative nonuniform tails must be interpreted accordingly.

## Evidence

The formula and support/ownership boundary have direct numerical regressions.
The protocol records the estimands, rerun requirements, independent sampling
unit, selected-null contamination check, held-out tuning/evaluation separation,
alpha sensitivity and tail-resolution limitations. Existing raw external-study
counts were recomputed in `external_contract_check.json`.

## Links

- [[empirical-null-calibration-reference-law-contract]]
- [[diffusion-nnls-versions-and-support-review-20260909]]
- [[selected-hierarchy-external-calibration-contract-20260602]]
- [[open-mathematical-questions]]
