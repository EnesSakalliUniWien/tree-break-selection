---
title: Eigensystem API Reference Notes
type: source
status: reviewed
updated: 2026-09-28
sources:
  - raw/inbox/eigensystem-reference-notes-20260928.md
  - reports/mathematical_construction_review_20260928/reproduce_findings.py
  - reports/mathematical_construction_review_20260928/findings.json
tags:
  - eigenvectors
  - reference
  - geometry
---

# Eigensystem API Reference Notes

## Summary

Official SciPy API references distinguish symmetric eigenproblems, general
paired left/right eigenvectors, and principal angles between full subspaces.
The reference notes support terminology for [[eigenvector-calibration-construction-review]].

## Key Points

- `eigh` assumes a symmetric or Hermitian operator; it does not verify symmetry.
- `eig` defines left and right eigenvectors of a general square operator.
- `subspace_angles` compares full column spaces, including unequal ranks.
- Browsed current documentation is newer than the installed SciPy 1.16.3
  runtime. Numerical results were reproduced locally and recorded separately.

## Evidence

- `raw/inbox/eigensystem-reference-notes-20260928.md` contains source URLs,
  access context and concise paraphrases.
- `reports/mathematical_construction_review_20260928/findings.json` records
  the installed runtime and principal-angle/eigenvector residual calculations.

## Links

- [[eigenvector-calibration-construction-review]]
- [[eigenvalue-block-tree-comparison]]
