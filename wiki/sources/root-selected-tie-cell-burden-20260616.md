---
title: Root Selected Tie Cell Burden 2026-06-16
type: source
status: draft
updated: 2026-08-10
sources:
  - benchmarks/diagnostics/calibration/root/selected/root_selected_tie_cell_burden.py
tags:
  - source
  - diagnostics
  - root
  - selected-region
  - tie-cell
---

# Root Selected Tie Cell Burden 2026-06-16

## Summary

`root_selected_tie_cell_burden.py` describes the discrete part of the selected
root event. It aggregates tied construction choices with

\[
\sum_t \log m_t,
\]

where \(m_t\) is the tied-minimum multiplicity at construction step \(t\).
The burden is a conditioning coordinate, not a calibrated p-value or monotone
penalty.

## Key Points

- Raw tie multiplicity and the selected pair's deterministic tie rank are
  recorded as distinct coordinates.
- Case rows retain tie-step counts, multiplicity summaries, rank fractions,
  and side-asymmetry measures.
- Relationship tables are descriptive diagnostics only.
- The output feeds the mixed selected-region law.

## Evidence

- `build_root_selected_tie_cell_burden_rows` derives case-level burden and rank
  coordinates from root summaries and merge-margin rows.
- The module writes burden, relationship, summary, and manifest contracts via
  the shared diagnostic reporting layer.
- No tracked result capture or focused current test supports the former
  overlap-run numerical claims, so this page remains draft.

## Links

- [[root-selected-region-overlap-case-family-20260616]]
- [[selected-neighborhood-measurability-law]]
