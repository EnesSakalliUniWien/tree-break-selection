---
title: Independent Audit of the Three Sibling-Test Constructions
type: source
status: reviewed
updated: 2026-11-09
sources:
  - raw/inbox/choir-cell-lines/three-tests-audit-20261109.md
  - raw/inbox/choir-cell-lines/three-tests-supplied-construction-20261109.txt
  - raw/inbox/choir-cell-lines/three-tests-wolfram-20261109.json
  - raw/inbox/choir-cell-lines/three-tests-rerun-evidence-20261109.json
  - raw/inbox/choir-cell-lines/three-tests-audit-scripts-20261109.py.txt
tags:
  - calibration
  - sibling-gate
  - mathematical-audit
---

# Independent Audit of the Three Sibling-Test Constructions

## Summary

The supplied comparison overstates the error control of the ceiling-3 and
best-cut-tail proposals. Ceiling 3 applies to a scalar unimodal population cut
at its mode, not arbitrary selected cuts. A continuous unimodal counterexample
has population rho 7.39984. Wolfram confirms C's variance 1/40 and sphere-3
lambda 13, but also confirms B's boundary-variance mismatch and the impossibility
of uniform marginals in every direction above three dimensions. Fresh reruns
reproduce all 54 continuous-benchmark rows and the saved split calculations.

## Key Points

- Current A uses context weights in parent size and projection dimension; its
  90% rule uses observed contrast energy. Fixed membership alone does not give
  the claimed chi-square law with fitted covariance and selected rank.
- For a uniform fixed cut, the asymptotic rho variance is 19.2/n, versus 12/n
  from B's noncentral chi-square. B's counted tail rejects 40/1,000 sphere-3
  boundary samples at nominal 0.01 when n=100.
- A continuous unimodal long-tail example gives B/C rejection counts
  1,996/1,997 out of 2,000 at n=400, nominal 0.01. Both reject every replicate
  at n=1,600. These are fixed-k tail checks, not complete tree fits.
- C's exact R transformation is correct. Its influence-field correlation has
  an absolute cubic angular term, so the cited C2 smooth-process theorem does
  not apply directly. The EC envelope is not a proved survival bound.
- C's fresh uniform k=1, n=100 check rejects 75/1,000 at nominal 0.05. The
  high-dimensional 0-1/300 results concern searched lower bounds for maxima;
  they do not establish conservative control.
- Cover's count is correct up to its harmless included empty partition. A
  samplewise count alone does not justify the proposed union-bound tails.
- Neither fixed-k checks nor BH establish inference for selected descendants
  or error control of the final partition. No production replacement was made.
- All four full-cell A/B fits reproduce their earlier assignments exactly.
  Markov/Ward: A 954 clusters, ARI 0.548774, 41 recovered lines; B 217,
  0.832948, 119; C 173, 0.815201, 111. Default: A/B/C 75/64/57 clusters,
  all near ARI 0.028. C was run end to end as a diagnostic, with all 1,928
  replacement p-values independently verified across the two trees.

## Evidence

The raw report contains derivations, inspected code paths, confidence intervals,
primary literature links, rerun commands and scope limits. The Wolfram capture
contains the actual evaluator request and returned exact expressions. The
rerun evidence capture preserves completed numerical summaries; detailed local
outputs are under `results/analyses/three_tests_audit_20261109/`.

## Links

- [[single-population-split-bound]]
- [[flat-cut-sup-count]]
- [[ceiling-rule-real-and-benchmark]]
- [[ceiling-rule-literature-comparison]]
- [[analytic-selected-root-laws]]
