---
title: Root Tie Rank Selected Null Simulation Pilot 2026-06-16
type: source
status: draft
updated: 2026-08-10
sources:
  - benchmarks/diagnostics/calibration/root/tie_rank/root_tie_rank_selected_null_simulation_pilot.py
tags:
  - source
  - diagnostics
  - root
  - calibration
  - selected-null
---

# Root Tie Rank Selected Null Simulation Pilot 2026-06-16

## Summary

`root_tie_rank_selected_null_simulation_pilot.py` is the first executable
selected-null generator for the discrete root tie-rank calibration path. It
builds iid Bernoulli null matrices matched to each binary overlap case's
sample count, feature count, and expected sparse-template marginal feature
rate, sends those matrices through the existing root selected-region replay,
and recomputes target-stratum support against the observed mixed root-law
rows. Version 2 propagates the reduced component/tie/edge/spectral stratum
contract and no longer emits a topology-derived bandwidth count.

## Key Points

- Generated null matrices pass through the canonical TBS context and root
  selected-region extractor.
- Successful roots are converted to mixed-law and calibration-feasibility
  rows before observed-target support is summarized.
- Failures are retained as explicit diagnostic rows rather than silently
  reducing the attempted replicate count.
- Target support output now follows
  `root_tie_rank_selected_null_simulation_pilot/v2`.

## Evidence

- The implementation uses `preloaded` generated null matrices so the canonical
  TBS context builder and root selected-region extractor remain unchanged.
- The retained code computes target-stratum support only after recomputing the
  reduced calibration key.
- No non-empty retained result capture or focused current test supports the
  former numerical pilot claims, so this page remains draft.

## Links

- [[root-tie-rank-calibration-feasibility-20260616]]
- [[root-selected-mixed-region-law-20260616]]
- [[selected-neighborhood-measurability-law]]
