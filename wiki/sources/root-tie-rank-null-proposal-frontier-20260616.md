---
title: Root Tie Rank Null Proposal Frontier 2026-06-16
type: source
status: draft
updated: 2026-08-10
sources:
  - benchmarks/diagnostics/calibration/root/tie_rank/root_tie_rank_null_proposal_frontier.py
tags:
  - source
  - diagnostics
  - root
  - calibration
  - proposals
---

# Root Tie Rank Null Proposal Frontier 2026-06-16

## Summary

`root_tie_rank_null_proposal_frontier.py` extends the root tie-rank
calibration path with explicitly labeled proposal families. Ordinary iid
Bernoulli rows may be calibration candidates, while enriched proposals remain
diagnostic unless their calibration role and likelihood-ratio contract make
them admissible external-null support.

Version 2 builds generated mixed-law rows directly from margin and tie-cell
evidence. It no longer accepts or propagates the retired neighborhood-topology
overlay.

## Key Points

- The diagnostic writes proposal roots, merge margins, tie rows, mixed-law
  rows, feasibility rows, target support, frontier summaries, and failures.
- Proposal reachability and calibration admissibility remain distinct.
- Importance proposal families record target-null probability, proposal
  probability, and log-weight metadata.
- Conditioning keys follow the reduced component/tie/edge/spectral contract.
- The observed mixed-law rows are an explicit predecessor input.

## Evidence

- The implementation attaches proposal metadata only after the canonical root
  margin and tie-cell builders return their rows.
- Feasibility support counts depend on explicit calibration roles rather than
  proposal-family labels alone.
- No non-empty retained result capture or focused current test supports the
  former numerical smoke-run claims, so this page remains draft.

## Links

- [[root-tie-rank-selected-null-simulation-pilot-20260616]]
- [[root-tie-rank-calibration-feasibility-20260616]]
- [[root-selected-mixed-region-law-20260616]]
- [[selected-neighborhood-measurability-law]]
