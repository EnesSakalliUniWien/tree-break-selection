---
title: Root Selected Spectral Tail Law With Legacy Overlay 2026-06-17
type: source
status: draft
updated: 2026-08-10
sources:
  - benchmarks/diagnostics/calibration/root/selected/root_selected_spectral_tail_law_panel.py
tags:
  - source
  - diagnostics
  - root
  - spectral
---

# Root Selected Spectral Tail Law With Legacy Overlay 2026-06-17

## Summary

`root_selected_spectral_tail_law_panel.py` expresses the selected-root
spectral-tail inference target directly. For each observed root, it records
selected spectral excess (S_{\mathrm{root}}), tie rank (T), selected-ratio
action (A), edge-margin action (E), and population-law status (H_u). It
looks for admissible selected-null or external-null support in the same coarse
((T,A,E,H_u)) stratum while leaving (S_{\mathrm{root}}) as the tail
variable.

Version 2 removes the retired topology-derived bandwidth coordinate and
requires callers to provide the feasibility table explicitly. Despite the
historical page name, the retained module no longer owns legacy overlay fields.

## Key Points

- The panel remains diagnostic-only and cannot rescue a production root split.
- Direct selected-null support uses a conservative plus-one empirical tail.
- Importance-weighted external support is reduced to an effective sample size
  before the conservative tail estimate is reported.
- Missing support remains an explicit fail-closed outcome.
- `root_selected_spectral_tail_law_panel/v2` records the reduced schema.

## Evidence

- `root_tail_stratum_key` now conditions on component, (T), (A), (E), and
  (H_u), but not the tail variable itself.
- The command-line contract requires an explicit feasibility input instead of
  pointing to the removed conditioned-topology join artifact.
- No non-empty retained result capture or focused current test supports the
  historical numerical overlay claims, so this page remains draft.

## Links

- [[local-marchenko-pastur-rule]]
- [[selected-neighborhood-measurability-law]]
