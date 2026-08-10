---
title: Root Selected Deformed MP Edge 2026-06-17
type: source
status: draft
updated: 2026-08-10
sources:
  - benchmarks/diagnostics/calibration/root/selected/root_selected_deformed_mp_edge_panel.py
  - benchmarks/diagnostics/calibration/root/selected/root_selected_mixed_region_law.py
  - benchmarks/diagnostics/calibration/root/tie_rank/root_tie_rank_null_proposal_frontier.py
tags:
  - source
  - diagnostics
  - root
  - spectral
  - population-law
---

# Root Selected Deformed MP Edge 2026-06-17

## Summary

`root_selected_deformed_mp_edge_panel.py` computes a diagnostic plug-in
deformed Marchenko--Pastur edge from captured selected-root bulk spectra. It
uses the Silverstein--Choi inverse map

\[
z(v)=-\frac1v+\gamma\int\frac{t}{1+t v}\,dH_u(t)
\]

and solves the right-edge condition on \((-1/t_{\max},0)\).

## Key Points

- The empirical \(H_u\) input excludes the leading identity-MP spikes from the
  captured root spectrum.
- Target and support rows retain identity and deformed edges, edge multipliers,
  and the corresponding spectral-excess values.
- The feasibility and \(H_u\)-observability inputs are required explicitly;
  there are no defaults to retired diagnostic artifacts.
- The result remains diagnostic-only and cannot rescue a production root split.

## Evidence

- `deformed_mp_upper_edge` solves the right-support edge numerically from the
  supplied population spectrum and aspect ratio.
- The panel writes observed and support-side deformed-edge tables for the
  selected-root spectral-tail diagnostic.
- No tracked result capture or focused current test supports the former
  accumulated-run numbers, so this page remains draft.

## Links

- [[root-selected-spectral-tail-law-with-legacy-overlay-20260617]]
- [[local-marchenko-pastur-rule]]
- [[selected-geometry-mp-integral-literature-20260602]]
