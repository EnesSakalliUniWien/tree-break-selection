---
title: Root Selected Region Overlap Case Family 2026-06-16
type: source
status: draft
updated: 2026-08-10
sources:
  - benchmarks/diagnostics/calibration/root/selected/root_selected_region_margins.py
tags:
  - source
  - diagnostics
  - root
  - selected-region
---

# Root Selected Region Overlap Case Family 2026-06-16

## Summary

`root_selected_region_margins.py` replays root-child construction inequalities
for selected roots. It separates tie-heavy discrete selected regions from
smooth first-order margin regions and records the root, edge, tie-rank, and
spectral coordinates needed by downstream diagnostic laws.

## Key Points

- The panel is diagnostic-only and does not define a production rescue rule.
- Root-child merge margins and tied minimum choices determine whether the
  selected region is discrete or smooth.
- Deterministic selected tie-rank coordinates are retained separately from raw
  tie multiplicity.
- The emitted case-level rows feed the tie-cell burden and mixed-region panels.

## Evidence

- `collect_observed_root_selected_region_row` runs the canonical tree context
  and returns the selected-root row plus merge-margin records.
- The module versions its row, relationship, and summary table contracts and
  writes them through the shared diagnostic bundle.
- No tracked result capture or focused current test supports the former
  seven-case numerical claims, so this page remains draft.

## Links

- [[root-selected-region-margins-20260603]]
- [[root-selected-tie-cell-burden-20260616]]
- [[selected-neighborhood-measurability-law]]
