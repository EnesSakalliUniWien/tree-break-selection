---
title: Lower-Dimension Degrees-of-Freedom Audit
type: source
status: reviewed
updated: 2026-10-17
sources:
  - raw/inbox/lower-dimension-reference-law-notes-20260928.md
  - reports/lower_dimension_df_audit_20260928/README.md
  - reports/lower_dimension_df_audit_20260928/check_reference_law.py
  - reports/lower_dimension_df_audit_20260928/reference_law_summary.csv
  - reports/lower_dimension_df_audit_20260928/reference_law_verification.json
  - reports/lower_dimension_df_audit_20260928/implementation_audit.md
  - reports/lower_dimension_df_audit_20260928/prior_research_audit.md
  - tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_kernel.py
  - tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_reference_distribution.py
  - benchmarks/diagnostics/calibration/statistics/statistic_distribution_shape_panel.py
tags:
  - calibration
  - projection
  - degrees-of-freedom
---

# Lower-Dimension Degrees-of-Freedom Audit

## Summary

The current kernel records degrees of freedom as the actual retained projection
dimension, including complete numerical eigenvalue groups. This bookkeeping
does not establish a chi-square null law when the projection or rank is selected
using the tested data. Two separate audits cover implementation and
previous research; new Gaussian controls isolate these distinct claims.

## Key Points

- For a fixed orthonormal rank-k projection of an isotropic Gaussian null,
  squared projected norm has chi-square k law. Rank counts directions, not
  distinct eigenvalues; a complete repeated eigenvalue block contributes its
  entire retained dimension.
- A one-dimensional direction chosen to align with an eight-dimensional tested
  vector retains its entire squared norm, whose law is chi-square eight. This
  is an illustrative counterexample, not the production PCA algorithm or a
  recommendation to replace every selected rank with ambient dimension.
- In 100,000 independent Gaussian contrasts at alpha .01, the production
  contrast-energy rule at fraction .8 gives 2,148 rejections despite correctly
  recording selected k. Fixed rank one gives 957, an independent random
  direction gives 1,012, and full rank eight gives 985. These controls fix
  covariance and omit tree inference and fitted inflation.
- The June statistic-shape and null-law panels already separated projection
  df from covariance-derived scale/effective df and adaptive selection. The
  September independent-data PCA comparison also shows distortion with rank
  fixed at one; rank variation is not required for a selection effect.
- The earlier held-out energy-.90 experiment found excess across distinct and
  partially tied spectra, while complete-group retention fixed the fully tied
  example. Its full-pipeline comparisons disabled energy selection. The new
  .80 control confirms an existing finding under another bounded design.
- Energy selection defaults to disabled in orchestration but is explicitly
  enabled at .90 in the pancreas benchmark. That application was not rerun.
- With a fixed projection and nonidentity null covariance, the relevant
  quadratic-form weights come from projected null covariance. Satterthwaite
  matches moments using both scale and effective df; it is not an exact
  selected-projection correction.
- The previous fixed-rank result and the new variable-rank counterexample are
  compatible. Neither supplies the full selected-tree calibration law.

## Evidence

The root diagnostic uses seed 2026092803 and checks its vectorized formulas
against 1,024 production-kernel calls. A separate tied-eigenvalue case expands
target rank one to actual df two. Source hashes match before/after execution,
and output includes pointwise Wilson intervals. The root baseline passes 30
tests; these tests establish behavior rather than an adaptive null theorem.
New external-source notes record the fixed-quadratic-form and selective
inference references and their applicability limits. No production behavior
or defaults were changed.
The implementation audit's overlapping baseline passed 41 tests. Earlier June
locked/rerun CSVs are available; some June smoke figures survive only as wiki
summaries and are identified as such in the research audit.

## Links

- [[projected-wald-statistic]]
- [[statistic-distribution-shape-panel-20260613]]
- [[null-law-decomposition-panel-20260613]]
- [[selected-pca-projected-wald-validation]]
- [[eigenspace-benchmark-calibration-reassessment]]
- [[independent-projection-vandermonde-control]]
- [[passthrough-pvalue-normalisation-audit]]
- [[open-mathematical-questions]]
