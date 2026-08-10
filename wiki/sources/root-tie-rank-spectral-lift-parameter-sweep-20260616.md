---
title: Root Tie Rank Spectral Lift Parameter Sweep 2026-06-16
type: source
status: draft
updated: 2026-08-10
sources:
  - benchmarks/diagnostics/calibration/root/tie_rank/root_tie_rank_spectral_lift_parameter_sweep.py
tags:
  - source
  - diagnostics
  - root
  - spectral
  - generator
---

# Root Tie Rank Spectral Lift Parameter Sweep 2026-06-16

## Summary

`root_tie_rank_spectral_lift_parameter_sweep.py` tests whether generated roots
can reach an observed selected spectral excess while satisfying selected
tie-rank and action-edge conditions. Passing rows are standalone
proposal-generator reach evidence; they do not feed a retained tail evaluator
or define production p-values.

## Key Points

- The sweep supports dense-plus-spike, two-block, coherent rank-one, and
  target-conditioned coherent proposal families.
- Selected spectral excess remains the reach target rather than a conditioning
  coordinate.
- Target-conditioned settings prevent generated rows from supporting unrelated
  targets.
- Version 3 replaces tail-oriented statuses and follow-up steps with standalone
  root-metric reach diagnostics.
- All outputs remain diagnostic-only.
- The sweep requires an explicit observed mixed-law table and does not select
  a historical result capture.

## Evidence

- `build_spectral_lift_sweep_target_rows` evaluates action-edge/tie eligibility
  and the remaining spectral lift for each setting and target.
- `summarize_spectral_lift_sweep_targets` reports target coverage and spectral
  reach without claiming calibration.
- No tracked result capture or focused current test supports the former smoke
  run numbers, so this page remains draft.

## Links

- [[selected-neighborhood-measurability-law]]
