---
title: Root Tie Rank Target Conditioned Importance Frontier 2026-06-17
type: source
status: draft
updated: 2026-08-10
sources:
  - benchmarks/diagnostics/calibration/root/tie_rank/root_tie_rank_target_conditioned_importance_frontier.py
tags:
  - source
  - diagnostics
  - root
  - spectral
  - importance-sampling
---

# Root Tie Rank Target Conditioned Importance Frontier 2026-06-17

## Summary

`root_tie_rank_target_conditioned_importance_frontier.py` searches proposal
settings for likelihood-ratio-weighted external-null rows that match an
unsupported target's root component, tie band, selected-ratio action band, and
edge-action band. Matching candidates proceed directly to the selected-root
tail evaluation.

## Key Points

- Generated rows carry `conditioning_target_case_id` to prevent cross-target
  borrowing.
- Proposal rows preserve target-null probability, proposal probability, and
  importance log-weight metadata.
- The optional rejection mode retains only candidates matching the target
  conditioning stratum; its target-specific normalizer remains a diagnostic
  concern rather than a production p-value rule.
- Version 2 renames the former `pre_topology_*` schema and statuses to the
  final conditioning contract and removes replay-oriented outcomes.

## Evidence

- `_conditioning_stratum_key` uses root component, selected tie rank,
  selected-ratio action, and edge action while leaving spectral excess as the
  tail variable.
- `build_target_rows` reports conditioning-stratum hits and marks matched rows
  ready for selected-root tail evaluation.
- No tracked result capture or focused current test supports the former smoke
  run numbers, so this page remains draft.

## Links

- [[root-selected-importance-tail-support-20260617]]
- [[root-selected-spectral-tail-law-with-legacy-overlay-20260617]]
- [[selected-neighborhood-measurability-law]]
