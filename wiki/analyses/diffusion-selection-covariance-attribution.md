---
title: Diffusion Selection Covariance Attribution
type: analysis
status: reviewed
updated: 2026-11-02
sources:
  - raw/inbox/diffusion-attribution-20261102.md
  - raw/inbox/diffusion-attribution-wolfram-20261102.md
  - benchmarks/diagnostics/calibration/diffusion_attribution.py
  - benchmarks/diagnostics/calibration/selected/hierarchy/diffusion_attribution_restoration.py
  - tests/core/test_diffusion_attribution.py
  - tests/core/test_diffusion_attribution_restoration.py
  - tree_break_selection/tree/distributions.py
tags:
  - diffusion
  - selection
  - covariance
---

# Diffusion Selection Covariance Attribution

## Summary

Selection alone is sufficient for the measured Gaussian root failure, while
estimated covariance adds a separate finite-sample reference mismatch. Native
projection has no direct root effect in this panel, but affects empirical
inflation through lower-node records.

## Details

Holding entire diffusion trees, root groups and child counts fixed, compare
selected data against fresh independent Gaussian data. With full 24-dimensional
contrasts, known identity covariance, zero branches and inflation 1, nominal
.01 root rejection counts are 50/50 versus 0/50. No PCA or covariance estimation
is needed for the selection failure.

Production guarded covariance gives 8/50 independent rejections:7/25 dense
pooled,1/24 diagonal pooled and 0/1 parent-total fallback. For dense pooled
covariance,55Q/(78×24) follows F(24,55) under a fixed independent Gaussian
partition before numerical ridge. Wolfram computes .2342701 rejection
probability from using the nominal .01 chi-square 24 threshold instead.
Diagonal pooling and parent-total fallback have different reference laws;
the F correction must not be generalized to them or selected data.

Restoring native production projection uses the same basis/rank across
covariance arms. Every root spans all 24 dimensions, so raw Q is unchanged.
Empirical inflation changes through nonroot records: selected-data guarded
descriptive rejections fall from 50 to 34 with full features and to 28 with
native projection. Edge-role assignments are frozen within each testing
matrix. Only 34 selected roots and zero independent roots are focal/eligible;
predictions on skipped roots are not executed gate tests.

## Evidence

[[diffusion-attribution-evidence]] links the complete methods, crossed counts,
raw artifacts and Wolfram inputs/outputs. Independent review verified exact
agreement with prior native factorial root p-values and distinguished 28 root
rejections from 29 final fragmented datasets due to downstream traversal.

## Links

- [[diffusion-attribution-evidence]]
- [[diffusion-topology-branch-rank-isolation]]
- [[analytic-selected-root-laws]]
- [[mathematical-glossary]]

## Open Questions

- What selected-region law captures the actual adaptive diffusion event?
- How should that law incorporate the three covariance regimes and fitted
  empirical-null support/scale rather than treating chi-square as exact?
- Which projection effects emerge when native root spans are reduced, and how
  do resulting conditional tests interact with traversal and multiplicity?
