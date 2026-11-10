---
title: Subtree Distance Topology Families and Conditional References
type: analysis
status: reviewed
updated: 2026-11-09
sources:
  - raw/inbox/choir-cell-lines/subtree-topology-reference-20261109.md
  - raw/inbox/choir-cell-lines/subtree-topology-mathematics-20261109.md
  - benchmarks/diagnostics/calibration/selected/hierarchy/subtree_topology_geometry.py
  - benchmarks/diagnostics/calibration/selected/hierarchy/subtree_gaussian_reference.py
  - benchmarks/diagnostics/calibration/selected/hierarchy/subtree_topology_tests.py
  - tests/statistics/71_test_subtree_gaussian_reference.py
  - tests/statistics/72_test_subtree_topology_geometry.py
  - tests/statistics/73_test_subtree_topology_tests.py
  - raw/inbox/choir-cell-lines/subtree-clustering-research-20261109.md
  - benchmarks/diagnostics/calibration/selected/hierarchy/subtree_topology_candidate.py
  - benchmarks/diagnostics/calibration/selected/hierarchy/subtree_family_reference.py
tags:
  - subtree
  - topology
  - distances
  - calibration
---

# Subtree Distance Topology Families and Conditional References

## Summary

The relevant count is a family conditional on a fixed subtree, its distances,
the reconstruction rule, and the treatment of ties. Eigenspace rank alone
does not determine that family. Candidate topology and fitted distances must
both be compared; zero-length edges can make several binary refinements one
unresolved distance topology.

## Details

For each complete spectral block or leading union, the geometry module returns
the row-space projector, projected coordinates, distances and strict affine cuts.
Rotations within a repeated eigenvalue block preserve that block's geometry.
The Ward family explores every permitted minimizing merge and removes duplicate
histories. The union across spaces counts alternatives supported by at least
one view; the intersection finds common outputs. Crossing splits can describe
distinct structure that cannot all be edges of one tree.

`fit_additive_distances` fits nonnegative edge lengths for a supplied topology
and reports fitted pair distances, RMSE and largest residual. Complementary
root edges share an incidence column and are represented by their length sum.
Its loss is different from Ward's merge criterion. Rooted RF, edge splits,
cophenetic distances and additive residuals answer different comparison questions.

Cut enumeration is bounded at ten leaves. Rooted topology enumeration, tied-Ward
families and additive fitting are bounded at six. The default Ward family
tolerance is 1e-10 of the maximum original squared distance. Near-tie families
are deterministic diagnostics, not statistical confidence sets. An additive
positive-edge metric identifies an unrooted tree; an ultrametric model supplies
rooted hierarchy information. Arbitrary Euclidean subspaces need not satisfy
either model.

The Gaussian core conditions on a subtree mean, centered feature Gram and
outside rows. A supplied event also conditions on selection. The symmetric
orthogonal kernel and random-position path produce a weak upper rank with
all path positions included. Numerical failures remain unresolved; native
floating-point event membership is not certified.

The fixed-contrast Beta helper requires fixed membership and Gram-only rank;
its all-label bound covers label search, not selecting membership. The
integration provides native raw sibling replay, paired topology score
differences, full-span affine cut maxima, and maxima over regenerated spectral
Ward families. Candidate tails share the maximum reference. This is complete
Gaussian-null control conditional on the specified event, not strong FWER,
FDR or complete TBS selection control. Ward node paths must be prespecified;
additional selection among paths needs a larger event.

The native score adapter refits Ward, NNLS, spectral context, guarded covariance
and .90 adaptive contrast rank. It retains the raw sibling score and does not
use the empirical-inflation p-value. The alternative controls use the existing
full-space sibling score with zero branch lengths and no adaptive rank.

Bootstrap frequencies depend on the resampling law and its units. They can
sample within a fixed geometric family or reach topologies outside it. Simple
row relabeling preserves an equivariant distance procedure except for arbitrary
tie handling; it is not the conditional Gaussian orientation experiment.

The later research implementation adds a separate bounded candidate selector
for up to 256 leaves, with all-pair distance fits but incomplete topology
coverage. Exact exhaustive bounds above remain unchanged. A new common-event
Ward family module extends beyond the earlier complete-null maximum: under
independent homogeneous Gaussian true-null nodes and fixed outside observations,
conditional node ranks followed by Holm support partial-null family control.
It does not establish diffusion/search/gate/traversal inference. Its development
power and the candidate's runtime failed promotion prerequisites; see
[[subtree-clustering-production-gates]].

## Evidence

[[subtree-topology-reference-evidence]] records exact bounded examples, five
Wolfram evaluations, Gaussian calibration panels, 138 integration attempts,
four exploratory power follow-ups, and 12,000 bootstrap/permutation replicates.
The combined source retains the canonical row-order correction and separates
nominal calibration, geometric recovery, score power and numerical limitations.

## Links

- [[subtree-topology-reference-evidence]]
- [[subtree-clustering-research-evidence]]
- [[subtree-clustering-production-gates]]
- [[gaussian-root-calibration-and-cut-geometry]]
- [[subtree-spectrum-topology]]
- [[subspace-distance-ties]]
- [[root-conditional-sampler]]
- [[descendant-chain-family-test]]

## Open Questions

- Which model and distance uncertainty region make a family a confidence set
  for the population topology rather than a set of allowed algorithm outputs?
- How can larger subtrees be covered without presenting sampled candidates
  as exhaustive or replacing the null correction with an observed count?
- What power does the retained native score have across imbalance, overlap,
  high dimension and model violations at practical budgets?
- How should additional TBS gate selection, diffusion topology, traversal and
  heterogeneous subtree nulls enter an end-to-end error-control procedure?
