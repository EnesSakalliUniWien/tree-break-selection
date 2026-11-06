---
title: Diffusion Selected Law Open Question Answers
type: analysis
status: reviewed
updated: 2026-11-02
sources:
  - raw/inbox/diffusion-selection-event-answers-20261102.md
  - raw/inbox/diffusion-boundary-vandermonde-answers-20261102.md
  - raw/inbox/diffusion-covariance-traversal-answers-20261102.md
  - raw/assets/diffusion-selection-event-20261102/reentry.json
  - tree_break_selection/space_separation/diffusion.py
  - tree_break_selection/tree/distributions.py
  - benchmarks/diagnostics/spectral/vandermonde.py
tags:
  - diffusion
  - selection
  - roots
  - covariance
---

# Diffusion Selected Law Open Question Answers

## Summary

The conditional reference measure is specified for known identity covariance;
complete recovery of its diffusion selection region is not yet implemented.
Actual full-history selection can re-enter in one feature, Newton alone cannot
establish complete boundary recovery, and the current operator is outside the
PSD Toeplitz Vandermonde theorem. Estimated covariance, fitted inflation and
traversal require additional conditional constructions and error targets.

## Details

### 1. What preserves the observed merge sequence?

Let the observed root groups define unit contrast vector e. Under independent
Gaussian rows with known identity feature covariance and contrast null
e' E[X]=0, set z=X'e, r_obs=||z||, u=z/||z|| and A=(I-ee')X. Conditional on A
and u before selection, r has chi_p density. Refit the entire constructor at
X(r)=A+e r u'. Define S_M as the nonnegative radii where the actual method is
defined and returns the observed complete merge sequence M, including its
tie behavior. For positive mass, the model-level selected tail is

\[
p_M=\frac{\int_{S_M\cap[r_{obs},\infty)}f_{\chi_p}(r)\,dr}
{\int_{S_M}f_{\chi_p}(r)\,dr}.
\]

At each observed step, the selected pair's average diffusion-Euclidean
distance must beat every other active pair. These inequalities are not the
quadratics of the earlier squared-Euclidean diagnostic. The actual map refits
neighbors, NNKDE pilot bandwidths, internal BGH grid choice/rounded dimension,
outer median epsilon, symmetrization, normalization, retained spectrum and
linkage. Specifying outer epsilon=median does not disable internal BGH.

Full history is stronger conditioning than a root partition or unordered
clades. Holding extra discrete states fixed can define a further explicit
selection event. Holding exact continuous bandwidths or eigenvectors fixed
can constrain r to zero-mass sets and cannot silently be used with this
positive-mass radial law.

### 2. Can there be multiple roots, and does dimension determine the count?

The actual continuous preset has a reproducible full-sequence re-entry case:
seed 20261105, n=9, p=1, retained diffusion dimension 1. Membership is true at
r/r_obs=1, false at 1.1, and true at 1.2, with matching neighborhoods checked.
Four localized numerical transitions occur near ratios .9479474, 1.0707273,
1.1195953 and 1.3991668. This establishes that multiple changes do not require
high feature dimension. It is not a complete enumeration of S_M.

The gap changes the order of two disjoint merges while preserving the root
partition. Therefore it does not establish disconnected root-partition-only
selection; the union over compatible histories can bridge the gap. Neighbor
lists, internal KDE discrete choices and retained rank agree at each bracket's
endpoints. No general interval-count formula or monotonic relation to n, p or
rank has been established.

### 3. What can Newton–Raphson do?

Newton can refine a locally isolated smooth boundary. It cannot show that all
boundaries were found. Simple eigenpair derivatives need nonzero spectral gaps;
whole separated projectors are appropriate at exact ties when distances are
rotation-invariant. Numerical tie grouping introduces additional switches.

A certified procedure needs complete branch coverage, validated comparison
enclosures, tangency and jump handling, and treatment of the unbounded tail.
Ordinary bisection of sampled membership changes does not certify unique
zeros. Compactifying r with a trigonometric substitution does not resolve
the endpoint or turn the pipeline into a polynomial.

Unresolved chi mass must be compared with selection probability. If certified
selected numerator/denominator masses are N,D with D>0 and unresolved null
mass is at most delta, coarse safe p bounds are N/(D+delta) and
min(1,(N+delta)/D). A small unconditional omitted tail alone is insufficient.

### 4. Does the trigonometric/Vandermonde construction apply?

The theorem used by the existing diagnostic requires Hermitian PSD Toeplitz
input. The actual nonzero symmetric zero-diagonal diffusion operator has trace
zero, so it cannot be PSD: nonnegative eigenvalues with zero sum would all be
zero. The n80,p24 attribution-seed probe also finds Toeplitz residual 1.2624
and minimum eigenvalue -.31748. Its coordinate Gram is PSD within roundoff
but has Toeplitz residual 1.3974.

Consequently direct application is excluded. A structured circulant model,
Toeplitz fitting or changing the operator would be a different declared model,
requiring fresh approximation/selection analysis. Fourier coordinates alone
do not create the missing structure.

### 5. How do covariance and fitted inflation enter the law?

Dense pooled F, diagonal sum-of-F and parent-total Beta references have
different independent-Gaussian assumptions. They cannot simply replace
chi-square after tree selection. Unknown unrestricted covariance remains a
nuisance problem; the existing unknown-scalar construction requires its
stronger common-mean Gaussian null.

Once a nuisance-free path X(r), base density f0 and selection event S are
justified, a deterministic fitted statistic T can be handled by recomputing
all covariance, projection, rank, edge-role and empirical-scale calculations
along that path:

\[
p_T=\frac{\int_S\mathbf1\{T(X(r))\ge T(X(r_{obs}))\}f_0(r)\,dr}
{\int_Sf_0(r)\,dr}.
\]

This upper-level set may itself be disconnected. It is not the existing
chi-square plug-in formula, and freezing the observed fitted scale changes
the statistic. Undefined support and selection of which test is reported must
be covered by the conditioning/reporting contract. Computing this integral
for the actual pipeline remains unimplemented.

### 6. Does validity extend to projection, multiplicity and traversal?

Not automatically. A root conditional test, all-node FDR, selective levelwise
FDR and global-null final fragmentation are different targets. Native edge
TreeBH and sibling BH on edge-reachable depth families need a verified match
to the relevant p-value/dependence assumptions. Mean-equality edge nulls are
not automatically nested homogeneity hypotheses. This is a missing theorem
mapping, not a proof that every hierarchical procedure is impossible.

If final fragmentation requires rejection of a valid root gate, its global-null
probability is bounded by that root test's level. Current passthrough can cross
a nonrejected root sibling gate, so a root-only argument does not establish
that bound. Preserving passthrough needs additional simultaneous/hierarchical
control and a definition linking rejected hypotheses to reported boundaries.
No traversal or threshold change was made by this investigation.

## Evidence

[[diffusion-open-question-evidence]] links all three agents' primary notes,
source/theorem references and reproducible probes. Parent verification used
uninstrumented diffusion calls to reproduce the three decisive radial points
and unchanged root partition, and matched recorded source hashes. Structural
probe assertions match the native diffusion coordinate Gram.

## Links

- [[diffusion-selected-root-diagnostic]]
- [[diffusion-open-question-evidence]]
- [[diffusion-selection-covariance-attribution]]
- [[analytic-selected-root-laws]]
- [[diffusion-nnls-calibration-open-points]]
- [[vandermonde-reference-notes]]

## Open Questions

- Extend the partial exact-real certificates in [[diffusion-selected-root-diagnostic]]
  to informative conditional mass coverage; reconcile native positive self-distance
  artifacts before claiming validity for its observed merge history.
- Quantify the cost of full-history conditioning relative to root-only events;
  the demonstrated history gap is bridged by an unchanged root partition.
- Establish unknown-covariance nuisance treatment and fitted-statistic
  level-set computation under explicit models.
- Prove an error target for the present complete hierarchy/traversal procedure,
  or explicitly validate a separately specified procedure with such a guarantee.
