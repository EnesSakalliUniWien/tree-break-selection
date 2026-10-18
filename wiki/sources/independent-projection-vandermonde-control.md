---
title: Independent Projection and Vandermonde Control
type: source
status: reviewed
updated: 2026-10-17
sources:
  - reports/independent_projection_vandermonde_20260928/README.md
  - reports/independent_projection_vandermonde_20260928/control_protocol.md
  - reports/independent_projection_vandermonde_20260928/vandermonde_protocol.md
  - reports/independent_projection_vandermonde_20260928/control_evaluation/frozen_manifest.json
  - reports/independent_projection_vandermonde_20260928/control_evaluation/summary.csv
  - reports/independent_projection_vandermonde_20260928/control_evaluation/paired_summary.csv
  - reports/independent_projection_vandermonde_20260928/vandermonde/decompositions.csv
  - reports/independent_projection_vandermonde_20260928/vandermonde/structure_summary.csv
  - reports/independent_projection_vandermonde_20260928/verification.json
tags:
  - calibration
  - projection
  - vandermonde
---

# Independent Projection and Vandermonde Control

## Summary

Independent-data controls identify excess rejection from same-data PCA
eigenvector selection with rank fixed at one. Canonical global-null MP ranks
are almost constant, so those cells cannot separately identify rank-selection
distortion. Vandermonde application succeeds on structured covariances but
does not supply a unique lower-dimensional projection for this study.

## Key Points

- The control crossed basis and rank data sources,
  retaining all 10,000 original evaluation datasets and adding separate
  independent basis/rank training draws. All 486,000 rows completed.
- At fixed rank one, same-data PCA rejects global-null roots 33/1000 and
  37/1000 times (identity/AR), versus 11/1000 and 10/1000 with independent
  PCA. Matched differences are .022 and .027; descriptive paired normal
  intervals exclude zero. Each cell uses independent dataset replicates.
- Canonical evaluation MP rank is one on all 6,000 global-null node datasets.
  Independent MP has three rank-two events, with no changed rejections.
  Mixed-null roots and a prespecified, separately reported changed-threshold
  stress have varying ranks. Held-basis rank-source differences do not show a
  consistent direction; the stress is not a production threshold proposal.
- Exact population identity/AR covariances admit nonunique full-rank
  Vandermonde factorizations. Different anchors reconstruct the same
  covariance. Two known rank-two controls recover their atoms. Eight
  supported application cases reconstruct within 6.22e-15 relative error.
- Across 10,000 datasets each, unrestricted pooled and oracle-whitened PCA
  covariances all fail the Toeplitz check; sample-whitened PCA covariances
  all pass because they are scalar identities. No Toeplitz fitting is done.
- Vandermonde frequencies, weights, atoms, covariance eigenvalues and
  eigenvectors remain distinct. AR atom overlaps reach .288; QR is used
  before principal-angle/projector comparisons. Full-rank span agreement
  does not identify a preferred projection.

## Evidence

The independent verifier regenerates 30,000 streams, matches 30,000 original
node tests, checks all 486,000 chi-square calculations and 486 summary cells,
and reproduces 1,944 NumPy projection/MP calculations. All 117 pre-frozen
source hashes match after execution. Vandermonde has separate before/after
hashes and guarded structural/reconstruction checks. A read-only
review also checked nonorthogonal exact controls. No production changes,
dependency installation, commit or push occurred.

Known covariance, Gaussian sampling and fixed tree membership are essential
to this control. Estimated covariance, practical data splitting, categorical
observations, tree selection, and any feature-order assumption needed for a
Toeplitz model remain open. Numerical rank labels are not exact certificates.

## Links

- [[known-branch-skew-diagnostic]]
- [[vandermonde-reference-notes]]
- [[eigenvector-calibration-construction-review]]
- [[diffusion-nnls-calibration-open-points]]
