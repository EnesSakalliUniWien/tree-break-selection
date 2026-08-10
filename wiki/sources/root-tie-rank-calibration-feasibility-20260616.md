---
title: Root Tie Rank Calibration Feasibility 2026-06-16
type: source
status: draft
updated: 2026-08-10
sources:
  - benchmarks/diagnostics/calibration/root/tie_rank/root_tie_rank_calibration_feasibility.py
tags:
  - source
  - diagnostics
  - root
  - calibration
  - selected-region
---

# Root Tie Rank Calibration Feasibility 2026-06-16

## Summary

`root_tie_rank_calibration_feasibility.py` converts the mixed root selected
region
\[
\mathcal E_{\mathrm{root}}
=
\mathcal E_{\mathrm{margin}}
\cap
\mathcal E_{\mathrm{tie}}
\cap
\mathcal E_{\mathrm{rank}}
\]
into explicit calibration strata. The strata use the mixed-region component,
selected tie-rank band, root edge-margin band, and root spectral-ratio band.
The panel reports whether each stratum has enough selected-null support to
estimate a conditional root tail law. Version 2 removes the retired
topology-derived bandwidth coordinate.

## Key Points

- With target alpha `0.01`, plus-one empirical p-value resolution requires at
  least `99` selected-null rows per stratum.
- With target alpha `0.01` and relative tail standard error target `0.25`, the
  tail-precision requirement is `1584` selected-null rows per stratum.
- Only rows with an explicit admissible calibration role count as selected-null
  support.
- Conditioning keys now contain component, tie, edge, and spectral bands only.
- The mixed-region predecessor table is an explicit input rather than a dated
  capture default.

## Evidence

- The implementation calculates plus-one alpha resolution, binomial tail
  precision, support counts, and conservative empirical tail estimates.
- `root_tie_rank_calibration_feasibility/v2` records the reduced conditioning
  schema.
- No non-empty retained result capture or focused current test supports the
  former numerical case study, so this page remains draft.

## Links

- [[root-selected-mixed-region-law-20260616]]
- [[root-selected-tie-cell-burden-20260616]]
- [[selected-neighborhood-measurability-law]]
