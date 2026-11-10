---
title: Gaussian Root Calibration and Cut Geometry
type: analysis
status: reviewed
updated: 2026-11-09
sources:
  - raw/inbox/choir-cell-lines/root-validation-mathematics-20261109.md
  - raw/inbox/choir-cell-lines/root-validation-wolfram-20261109.json
  - raw/inbox/choir-cell-lines/root-validation-followup-20261109-evidence.json
  - raw/inbox/choir-cell-lines/root-validation-followup-20261109-sources.py.txt
  - raw/inbox/choir-cell-lines/three-tests-audit-20261109.md
tags:
  - gaussian
  - calibration
  - cut-geometry
  - selective-inference
---

# Gaussian Root Calibration and Cut Geometry

## Summary

A fixed Gaussian root admits an exact Beta reference for a fixed label contrast
in a shared subspace chosen from the feature Gram matrix. A union bound over all
deterministic label splits covers a selected root split. Its fixed-rank critical
ratio tends to 3; this is a conservative Gaussian threshold limit, distinct from
a population unimodality bound or a selected-pair conditional p-value.

## Details

Let Y be centered data from independent Gaussian rows with one common mean and
arbitrary unknown positive-definite common covariance. Choose V from G=Y'Y only,
with rank(YV)=k and 0<k<n-1. Use that same V for every candidate split. For a
normalized fixed split indicator u,

    R² = u' Proj(col(YV)) u | G ~ Beta(k/2, (n-1-k)/2).

The conditional row orientation is Haar. Gram-only rank choices and tie rules
preserve the argument; a rank depending on the observed contrast need not.
The exact finite-sample transformation is

    R² = n rho/(n-2+n rho),
    (n-1-k) R²/[k(1-R²)] ~ F(k,n-1-k).

For a root split selected by any rule using those rows and that shared score,

    p_all = min(1, (2^(n-1)-1) BetaSF(R²))

is superuniform under the global Gaussian null conditional on G, averaging over
which split is selected. Its deterministic family has all unordered nonempty
label splits. The realized Cover family depends on Y and cannot replace that
count using this unconditional per-label tail. Selected descendant unions also
require a separate argument. At k=n-1 the projector fills the centered row space
and R²=1 for every nonempty split. For fixed k as n grows, the all-label critical
R² tends to 3/4, so critical rho tends to 3. This is not a limit for k/n fixed.

For a fixed true-label split of two components with common within covariance,
rho=pi(1-pi)Delta² in retained coordinates. A ceiling-3 criterion therefore
requires Delta above 3.464, 4.330 and 7.947 at shares .50,.20,.05 to exceed the
population boundary. This exposes a balance penalty, but does not on its own
explain all selected four-component tree outcomes.

A scalar multiplier can correct scalar Gaussian covariance inflation under a
fixed-rank model. It cannot generally reproduce selection-induced distribution
shape, fitted-covariance tails and contrast-selected rank simultaneously. The
October 2026 Chen–Gao polar construction supplies a source-backed direction for
unknown shared covariance: condition on a common-mean Gaussian union's mean,
outside rows, feature Gram and the selected-pair event, then use a reversible
random-index rank. Retaining the deterministic recomputed TBS score follows from
the rank argument. The later opt-in diagnostic in
[[subtree-distance-topology-families]] implements that reference and replays the
native raw sibling statistic under a prespecified Ward path/split event. Its
finite panels do not validate full TBS eligibility or traversal. Global-root
whole-procedure Haar replay is a different derived option and answers a global
null rather than the selected mean-contrast question.

## Evidence

The raw derivation contains proofs, current source references and primary
literature links. Four Wolfram requests/responses verify the relevant algebra,
moments and limits. At nominal .01, spectrum-only rank rejects 120/12,000 and
contrast-energy rank with a substituted Beta tail rejects 256/12,000. Exhaustive
127-split maxima at n=8,k=2 with all-label Bonferroni reject 131/12,000; selected
flat-cut geometry and full-span saturation were checked independently.

## Links

- [[root-validation-generalisation-evidence]]
- [[three-tests-construction-audit]]
- [[single-population-split-bound]]
- [[analytic-selected-root-laws]]
- [[root-conditional-sampler]]
- [[subtree-spectrum-topology]]

## Open Questions

- How can the complete TBS root procedure be replayed under the Gaussian row
  orientation law at useful cost and with deterministic failure handling?
- Which null and full selection event should govern selected descendants and
  traversal, including pass-through and the final partition?
- Can a defensible model for fitted count/PCA inputs support this reference?
- Can power be retained for rare, overlapping or high-dimensional groups while
  accounting for selection without an empirical heuristic correction?
