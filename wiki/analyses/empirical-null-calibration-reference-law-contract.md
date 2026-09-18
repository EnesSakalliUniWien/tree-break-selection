---
title: Empirical-Null Calibration Reference-Law Contract
type: analysis
status: reviewed
updated: 2026-09-09
sources:
  - reports/calibration_restoration_20260909/validation_protocol.md
  - reports/diffusion_nnls_versions_20260909/README.md
  - reports/latest_changes_benchmark_review_20260909.md
  - tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflation_correction/empirical_null_inflation_estimation.py
  - tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflation_correction/types/inflation_model.py
  - tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/pair_testing/collection/child_parent_edge_metadata.py
  - tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/pair_testing/types/sibling_pair_record.py
  - tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflated_projected_wald_annotation/pipeline.py
  - tests/statistics/35_test_empirical_null_inflation_estimation.py
  - tests/statistics/40_test_child_parent_edge_metadata.py
  - tests/statistics/50_test_sibling_skip_annotation.py
  - wiki/sources/selected-hierarchy-external-calibration-contract-20260602.md
tags:
  - analysis
  - calibration
  - statistics
  - selection
---

# Empirical-Null Calibration Reference-Law Contract

## Summary

An estimated inflation factor does not by itself define a calibrated p-value.
For independent focal and calibration chi-square statistics with fixed unit
weights and one common scale, finite calibration uncertainty gives an exact F
reference law. The production same-selected-hierarchy records do not satisfy
that contract: they can overlap in observations, share stopping events, and use
data-derived hierarchy and weights. The 2026-09-09 restoration returns the historical plug-in chi-square p-values
and runs traversal-aligned sibling BH when internal support exists. The
selected-tail law remains unvalidated; `internal_admissible` means availability
under the support policy, not proven selective error control.

## Details

### Implementation review boundary

The initial 2026-09-09 review identified a production branch that bypassed
support-threshold decisions, plus a tuple dependency-group indexing failure.
The restoration removes that bypass and fixes tuple-group masks. The empirical
entry point now rejects exact-F models rather than mislabeling a plug-in
chi-square p-value as exact F. New regressions and the restored integration
checks cover these behaviors. See [[calibration-rule-restoration-20260909]].

### Clustering utility and calibration status

The diffusion NNLS review demonstrates that zero admissible internal support
can coexist with an exact clustering result under an explicitly diagnostic
gate. On the dense four-cluster signal case, pydiffmap NNLS plus
`fixed_coordinate_bh` returns ARI=1 at alpha 0.01, with the same tree and
branch lengths as the unsupported production result. This does not supply a
selected-tail law: fixed-coordinate gates also split null cases.

Removing only the benchmark's outer unsupported check leaves the production
sibling gates closed and returns one cluster on every compatible diagnostic
case. Calibration status should therefore remain distinct from empirical
clustering utility; a diagnostic assignment requires an explicit sibling
decision rule. That diagnostic review preceded the empirical-rule restoration. Zero internal
support still closes production sibling gates and remains an unsupported
benchmark outcome; the restoration does not insert a replacement gate. See
[[diffusion-nnls-versions-and-support-review-20260909]].

### Restricted exact law

Let the focal statistic and (m) calibration statistics satisfy

\[
T_u/(ac)\sim\chi^2_{\nu_u},
\qquad
T_q/(ac)\sim\chi^2_{\nu_q},
\]

mutually independently, with common positive reference scale (a), common
unknown inflation (c), and fixed unit calibration weights. With
(D=\sum_q\nu_q), independence and chi-square additivity give

\[
Y=\sum_q T_q/(ac)\sim\chi^2_D,
\qquad
\widetilde c=\frac{\sum_qT_q}{aD}.
\]

Therefore

\[
\frac{T_u/(ac\nu_u)}{Y/D}
=\frac{T_u}{a\nu_u\widetilde c}
\sim F_{\nu_u,D}.
\]

The exact mode requires a calibration observation-ID set for every calibration
parent. Those sets must be non-empty and pairwise disjoint, and the focal set
must be disjoint from their union. Calibration dependency groups and parent IDs
must also be distinct, all degrees of freedom positive, all weights exactly one,
and all reference scales equal. Its reported p-value is the larger of the exact
F upper tail and the unadjusted chi-square upper tail, which enforces the
one-sided no-deflation rule. Unequal deterministic weights are rejected because
neither Kish effective sample size nor Satterthwaite yields this exact law.

### Selected-hierarchy law

The production calibration sample is extracted from the hierarchy that also
selects the focal sibling pair. Nested blocked descendants can describe one
stopping event, observations overlap across nodes, and `sibling_null_weight` is
computed from selected child-edge evidence. An ordinary chi-square tail after
substituting an estimated scale ignores finite calibration uncertainty. The
exact F derivation also does not apply because its independence and fixed-weight
assumptions fail.

The restored production path uses the fitted scale in the historical plug-in
chi-square tail and enters sibling BH. The decision API retains the calibration
sample and support diagnostics, including the unresolved-law reason. Production
annotations contain adjusted statistics and empirical p-values. The necessary
next validation is held-out, complete-pipeline null and signal simulation.
Superuniformity in the declared null/selection context is the validity target;
nonuniformity alone does not demonstrate anti-conservative error.
Zero-dimensional records remain the deterministic (p=1) case.

### Ownership and support

`CalibrationSample` is the authoritative aligned representation of statistics,
reference scales, degrees of freedom, contexts, feature families, weights,
parent IDs, roles, and dependency groups. Its arrays are copied and read-only;
derived counts and summaries are not stored independently.

Each tested null-like record owns its own dependency group. A blocked descendant
walks upward to the nearest tested, non-significant ancestor edge path and
inherits that stopping event's group and group weight. Supported-group counts,
effective sample size, maximum group weight share, and leave-one-group scale
sensitivity operate on these groups. This prevents nested descendants from
manufacturing independent support, but it does not prove groups independent.

Support thresholds reject invalid values at construction. The immutable
support-policy snapshot pairs the thresholds with support-contract version 2;
serialization converts this typed object only at the output boundary. The old
stopped-or-null count is absent because every retained calibration record is
already strict null-like or edge-blocked.

## Evidence

- `tests/statistics/35_test_empirical_null_inflation_estimation.py` verifies the
  F reference, the unadjusted-p-value floor, observation ownership, invalid
  thresholds, immutable sample alignment, and the restored empirical p-values, strict support decisions and separation
  of the empirical and exact-F APIs.
- `tests/statistics/40_test_child_parent_edge_metadata.py` verifies that nested
  blocked descendants resolve to one stopping-event dependency group.
- `tests/statistics/50_test_sibling_skip_annotation.py` verifies production
  empirical p-values and requested support-threshold enforcement.
- [[selected-hierarchy-external-calibration-contract-20260602]] independently
  shows that scalar mean rescaling does not validate the selected-ratio tail.

## Links

- [[calibration-rule-restoration-20260909]]
- [[internal-calibration-q9-q10-q11-debug-20260605]]
- [[selected-hierarchy-null-support-contract]]
- [[selected-hierarchy-selection-geometry]]
- [[open-mathematical-questions]]

## Open Questions

- What conditional selected-tail law is valid after hierarchy construction,
  edge selection, stopping-event reuse, and focal sibling selection?
- What predeclared dependency-group thresholds are validated across feature
  families and parent-size regimes?
- Can unequal fixed weights be handled by an exact generalized chi-square-ratio
  computation, or should a bounded approximation be developed and labeled?
