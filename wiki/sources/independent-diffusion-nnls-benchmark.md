---
title: Independent Diffusion NNLS Benchmark
type: source
status: deprecated
updated: 2026-10-18
sources:
  - reports/independent_diffusion_nnls_20260928/README.md
  - reports/independent_diffusion_nnls_20260928/protocol.md
  - reports/independent_diffusion_nnls_20260928/evaluation/manifest.json
  - reports/independent_diffusion_nnls_20260928/evaluation/root_tails.csv
  - reports/independent_diffusion_nnls_20260928/verification.json
tags:
  - calibration
  - diffusion
  - nnls
---

# Independent Diffusion NNLS Benchmark

## Summary

The independent Gaussian API now fits adaptive Euclidean diffusion, linkage
topology, NNLS branch lengths and PCA exclusively on discovery measurements.
Inference uses the frozen context and separate measurements of the same leaves
with externally known common Gaussian covariance and independent rows.

The API was removed on 2026-10-18: nothing outside tests and the report scripts called it. The report evidence stays, but its scripts no longer run against current code.

## Key Points

- All 900 paired datasets and 18,000 runs completed, with 83,700 node records
  and no warnings/errors/skips. Sources and 2,700 fitted contexts remained
  unchanged. Six new integration tests and 168 focused tests pass.
- NNLS changed all 900 branch fits but none of 27,900 matched node tests.
  Known test covariance has no NNLS branch-length component. Fixed-topology
  diffusion and diffusion plus NNLS consequently give identical test results.
- Conditional diffusion/NNLS gives zero final false partitions in each set of
  250 null datasets at sibling alpha .01 and edge alpha .001. Mean signal ARI
  is .6453/.2512/.1553 for identity2/ar1_8/identity40, compared with direct
  .5663/.1773/.0000. Whole-method gains against current single-matrix methods
  also change independent data, covariance treatment and rank selection.
- Raw conditional root rejection is 6/250 at cutoff .01 for identity40
  diffusion/NNLS: pointwise interval [.0110,.0514], binomial p=.0412. No cell
  rejects after the prespecified 36-cell Holm correction. Finite non-rejection
  and zero final false partitions do not prove calibration.
- Conditional tails change 30 primary-alpha geometry/case cluster counts;
  diffusion/NNLS repeat the same outcomes. Signal recovery can decrease versus
  nominal tails. NNLS success can stop on relative cost change: 695/900 do so,
  with maximum optimality measure .1517. Solver success is not exact optimality.
- This Gaussian Euclidean result does not validate binary-Hamming diffusion,
  estimated covariance, dependent measurements or fitted branch covariance.
- The subsequent [[branch-distribution-variance-benchmark]] tests discovery-only
  internal rows. Branch-length states make NNLS affect the discovery PCA and
  node p-values; the unchanged-NNLS result above applies to this earlier
  leaf-only benchmark, not to every internal-distribution setting.

## Evidence

The frozen protocol, row-level outputs, matching/source-hash verification,
solver diagnostics and figure are in the report directory. The integration
rejects nonconverged NNLS and changed branch lengths/diffusion coordinates;
inference does not refit discovery models. Existing single-matrix callers
retain their input contract and calculation.

## Links

- [[independent-gaussian-production-benchmark]]
- [[branch-distribution-variance-benchmark]]
- [[diffusion-nnls-calibration-open-points]]
- [[adaptive-diffusion-nnls-method-library-audit]]
