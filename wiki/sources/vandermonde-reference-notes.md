---
title: Vandermonde Reference Notes
type: source
status: reviewed
updated: 2026-09-28
sources:
  - raw/inbox/vandermonde-reference-notes-20260928.md
  - benchmarks/diagnostics/spectral/vandermonde.py
  - tests/validation/calibration/197_test_vandermonde_diagnostic.py
tags:
  - vandermonde
  - covariance
---

# Vandermonde Reference Notes

## Summary

Primary-source notes for the PSD Toeplitz Vandermonde theorem and its guarded
numerical use in the independent-projection follow-up.

## Key Points

- The source is Yang and Xie, arXiv:1605.02431, Theorem 1.
- Applicability, uniqueness and numerical rank must be checked separately.
- General covariance matrices are not silently replaced by Toeplitz fits.
- This factorization does not itself establish a chi-square reference law.

## Evidence

The raw note records the consulted primary source and its construction.
The implementation and tests provide recoverable examples, unsupported-input
checks and an explicit full-rank nonuniqueness demonstration.

## Links

- [[independent-projection-vandermonde-control]]
- [[eigensystem-api-reference-notes]]
