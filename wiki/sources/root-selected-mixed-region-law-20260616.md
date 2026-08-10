---
title: Root Selected Mixed Region Law 2026-06-16
type: source
status: draft
updated: 2026-08-10
sources:
  - benchmarks/diagnostics/calibration/root/selected/root_selected_mixed_region_law.py
tags:
  - source
  - diagnostics
  - root
  - selected-region
  - tie-cell
---

# Root Selected Mixed Region Law 2026-06-16

## Summary

`root_selected_mixed_region_law.py` joins root selected-region margins and
rank-aware tie-cell burden into a single diagnostic table. It formalizes the
root event as
\[
\mathcal E_{\mathrm{root}}
=
\mathcal E_{\mathrm{margin}}
\cap
\mathcal E_{\mathrm{tie}}
\cap
\mathcal E_{\mathrm{rank}},
\]
where the rank term records the deterministic lexicographic position of the
selected pair inside each tied minimum merge set. Version 2 removes the retired
selected-neighborhood topology and bandwidth overlay. The panel does not
compute a calibrated p-value and does not promote traversal.

## Key Points

- Rows distinguish smooth, discrete tie-rank, mixed, and unsupported selected
  root regions from margin and tie-cell evidence.
- Discrete regions remain blocked until their tie-rank null law is calibrated;
  smooth regions remain diagnostic-only.
- Relationship output is descriptive and does not turn a correlated root
  coordinate into a production rescue rule.
- The output schema no longer exposes the retired neighborhood-locality
  overlay columns.

## Evidence

- The implementation requires only a root selected-region summary and a
  tie-cell burden table.
- `root_selected_mixed_region_law/v2` records the reduced contract explicitly.
- No non-empty retained result capture or focused current test supports the
  historical numerical claims, so this page remains draft.

## Links

- [[root-selected-region-overlap-case-family-20260616]]
- [[root-selected-tie-cell-burden-20260616]]
- [[selected-neighborhood-measurability-law]]
