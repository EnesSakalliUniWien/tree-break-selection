---
title: Subtree Context and Eigenvalue Roles in Native Statistics
type: analysis
status: reviewed
updated: 2026-11-09
sources:
  - raw/inbox/choir-cell-lines/subtree-context-statistics-20261109.md
  - reports/subtree_context_statistics_20261109/structure/report.md
  - reports/subtree_context_statistics_20261109/math/README.md
  - reports/subtree_context_statistics_20261109/native_summary.md
  - reports/subtree_context_statistics_20261109/relationships/README.md
  - benchmarks/diagnostics/calibration/selected/hierarchy/subtree_topology_context.py
  - reports/subtree_context_statistics_20261109/run.py
tags:
  - subtree
  - eigenvalues
  - hypotheses
  - production
---

# Subtree Context and Eigenvalue Roles in Native Statistics

## Summary

A surviving clade need not preserve its incoming edge comparison or outgoing
sibling partition. Compare the labeled node, immediate parent and child partition
before interpreting changes in statistics. Eigenvalue magnitudes have different
roles in diffusion geometry, covariance whitening and native PCA rank selection.

## Details

### Structural cases and hypothesis identity

For focal clade C, case 1 is C absent; case 2 is both original child clades absent;
case 3 is C present with its original parent absent; case 4 is C present with
both original children and the original parent absent. They overlap. Record
case 2 with and without C present. Roots have no parent comparison, and leaf
singletons cannot disappear when the leaf universe is fixed.

Presence is distinct from immediate adjacency. For surviving C, the four
same/different immediate-parent by same/different-child-partition combinations
identify which hypotheses persist. Changing the parent changes the incoming
contrast and complement. Changing children changes the outgoing contrasts.
If parent and children persist, global NNLS branches, their positive mean,
edge reachability, empirical-null support and traversal can still change.

The six-leaf trees `((((0,1),(2,3)),4),5)` and
`((((0,2),(1,3)),5),4)` supply a minimal case-4 witness. C=0123 has different
parents and crossing child partitions. Separate whole covariance eigendirections
can generate these trees; local rebuilding then rotates their inherited axes.

### What the statistic sees

The continuous sibling schematic is V0=(1/nA+1/nB) Sigma_hat,
g=1+(bA+bB)/(2 mean_positive_branch), V=g V0+ridge I,
z=whiten(d,V), T=||Qz||². A fixed projection without ridge gives inverse-g
scaling; native ridge and adaptive rank prevent a universal identity.
Collected raw statistics, displayed empirical-inflation statistics, active BH
decisions and traversal are distinct records. A missing clade has no native
p-value. An existing node can retain a finite raw p-value but lack an active
decision, or be unvisited after an ancestral boundary.

Incoming edges use the parent basis-half context; sibling tests use the all-leaf
context. Their eigenvalues cannot be interchanged. Covariance fallback and
full-node whitening also mean the basis-half split does not establish complete
independence. The same focal rows have the same raw covariance, while guarded
covariance and native projections can depend on parent/child context.

### Eigenvalue magnitudes

- Diffusion eigenvalues weight retained coordinates by max(lambda,0)^t.
  Squared-distance coefficients use exponent 2t, but actual distance energy
  also depends on eigenvectors and degree normalization.
- Covariance eigenvalues affect standardized contrasts when that covariance
  is used for whitening. Raw trace share is variation in declared feature units.
- Native null-whitened PCA eigenvalues order components, identify numerical
  ties and positive availability, and enter MP counts. They do not directly
  divide or weight T once Q and the retained rank are fixed. The .90 rule
  accumulates squared contrast coordinates, not eigenvalue weights.

Whole-block local reconstruction must discard a block that straddles the
numerical zero boundary rather than retain an unstable partial block. A zero
local covariance trace makes a variance share undefined. Rank mismatch is
separate from angular difference between equal-rank spans.

### Connections to calibration and production

The weaker raw covariance direction can retain a crossing factor omitted by a
larger-variance view. Better full-distance fitting can favor a different factor.
Factor partitions and their joint refinement therefore need separate recovery
metrics. Inherited coordinates, induced parent graphs and rebuilt child graphs
must likewise retain separate provenance.

The exact per-label Gaussian Beta reference fixes the union and contrast and
uses a shared Gram-only subspace/rank. The all-label correction covers label
selection on that fixed union. Its fixed-rank large-sample critical ratio 3 is
not a general population ceiling. Native contrast-adaptive rank, selected
descendants and changed hypotheses require further justification.

The conditional Ward diagnostic already replays the raw sibling score and
refits fixed-topology NNLS at accepted states. It does not validate the complete
jointly selected diffusion/subspace/branch/gate/traversal procedure. Observed
topology counts and bootstrap frequencies cannot replace that reference.
Production promotion still requires independent recovery, distance, runtime
and calibration evidence under the declared complete procedure.

## Evidence

[[subtree-context-statistics-evidence]] links 160 native fits on 16 reused
matrices, 16,128 dependent context comparisons, 5,120 saved alternative outcomes,
corrected covariance-span metadata, exact checks and native kernel controls.
The null and previous runtime failures remain unresolved; no preset changed.

## Links

- [[subtree-paths-and-missing-neighborhoods]]
- [[subtree-context-statistics-evidence]]
- [[subtree-spectral-alignment-and-crossing-clusters]]
- [[subtree-clustering-production-gates]]
- [[subtree-distance-topology-families]]
- [[gaussian-root-calibration-and-cut-geometry]]
- [[three-tests-construction-audit]]

## Open Questions

- Which declared alternative partitions should a production output preserve,
  and should incompatible partitions be represented by more than one hierarchy?
- Which common selection event covers changed node membership, parent/children,
  candidate spaces, branch fit, eligibility and traversal with useful power?
- Can candidate pruning preserve weaker grouping directions and whole-block
  invariance while meeting the existing runtime prerequisite?
