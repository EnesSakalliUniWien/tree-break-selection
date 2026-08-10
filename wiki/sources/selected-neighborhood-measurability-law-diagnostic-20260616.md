---
title: Selected Neighborhood Measurability Law Diagnostic 2026-06-16
type: source
status: draft
updated: 2026-08-10
sources:
  - benchmarks/diagnostics/calibration/selected/neighborhood/selected_neighborhood_measurability_law.py
tags:
  - source
  - diagnostics
  - measurability
  - traversal
  - neighborhood
---

# Selected Neighborhood Measurability Law Diagnostic 2026-06-16

## Summary

`selected_neighborhood_measurability_law.py` implements the selected-
neighborhood measurability law as a diagnostic evaluator. It applies direct
selected-family sibling p-values first, permits diagnostic rescue only when
direct measurement is blocked and supported interpolated priors plus topology
coherence are present, and otherwise fails closed with an explicit bottleneck
label.

## Key Points

- The implementation includes an unclipped
  `compute_child_interpolated_null_prior` helper. Probability-domain violations
  raise or mark the row invalid instead of being made valid by clipping.
- The row evaluator produces `split`, `diagnostic_rescue`, or `fail_closed`
  actions and reports bottlenecks such as
  `direct_measurable_not_significant`, `interpolated_prior_unavailable`,
  `support_bottleneck`, `selection_bottleneck_selected_nonnull_only`,
  `bandwidth_bottleneck`, and topology coherence failures.
- The refined candidate-audit implementation joins optional hold-out
  interpolation rows. Its output exposes
  traversal state, effective interpolation support, support weight, bandwidth
  scales, stable/signal neighborhood distances, topology coherence variables,
  selected-tree structural fallback fields, and interpolation behavior labels
  in one candidate-level table.
- The reduced output contract is versioned as
  `selected_neighborhood_measurability_law/v2`. No active focused test or
  non-empty retained result capture currently validates this benchmark-only
  diagnostic, so empirical result claims remain intentionally absent.

## Evidence

- `benchmarks/diagnostics/calibration/selected/neighborhood/selected_neighborhood_measurability_law.py`
  implements the evaluator, summary, CLI, and output manifest.
- The previously cited result directories are empty. Regenerate a result
  capture and add focused contract coverage before promoting this page back to
  `reviewed` or recording empirical conclusions.

## Links

- [[selected-neighborhood-measurability-law]]
- [[selected-neighborhood-bottleneck-law]]
