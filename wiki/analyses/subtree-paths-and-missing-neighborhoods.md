---
title: Subtree Paths and Missing Neighborhoods in Statistical Comparisons
type: analysis
status: reviewed
updated: 2026-11-09
sources:
  - raw/inbox/choir-cell-lines/subtree-path-statistics-20261109.md
  - reports/subtree_path_statistics_20261109/mathematics.md
  - reports/subtree_path_statistics_20261109/panel/summary.json
  - reports/subtree_path_statistics_20261109/panel/endpoint_summary.json
  - reports/subtree_path_statistics_20261109/panel/same_path_example.json
  - benchmarks/diagnostics/calibration/selected/hierarchy/subtree_path_context.py
  - reports/subtree_path_statistics_20261109/run.py
tags:
  - subtree
  - paths
  - hypotheses
  - production
---

# Subtree Paths and Missing Neighborhoods in Statistical Comparisons

## Summary

Two surviving nodes always have a connecting path in a tree. The original
intermediate clades and their neighborhood can disappear or move off that path.
The fixed full-data endpoint contrast survives, while native intermediate and
endpoint tests can compare different groups with different fitted quantities.

## Details

### Structural identity

Match nodes by exact labeled leaf sets. Nested endpoints retain their ancestry
relation; disjoint endpoints remain on separate branches. Along nested paths,
a surviving intermediate clade must remain on the path by set inclusion.
Across branches, a surviving original LCA can become an ancestor above a new
smaller LCA and leave the connecting path. This is distinct from global absence.

The measured neighborhood is the union of parents, children and siblings of
path vertices, excluding the path itself. Preserve role membership, global
absence and displacement outside the new neighborhood separately. Numeric
linkage IDs locate source records but do not identify clades across trees.

### Mean contrasts and covariance

For fixed normalized clade indicators a_C, oriented increments
w_j=a_Cj-a_C(j-1) satisfy sum(w_j)=a_end-a_start. Hence full-data mean contrasts
telescope. Under independent rows with common covariance Sigma,
Cov(w_i^T X,w_j^T X)=(w_i^T w_j) Sigma. The endpoint variance coefficient is
1/n_U+1/n_V-2|U intersect V|/(n_U n_V).

Full-data increments on a fixed nested chain are uncorrelated, and independent
for Gaussian rows. This still does not make their squared standardized scores
equal the endpoint score. On a cross-branch path, the two increments adjoining
a binary LCA are positively collinear when oriented along the path. Their cross
covariance must be included. The exact four-leaf example gives summed edge
scores 42 but direct endpoint score 4.5. It supplies no new selected-path p-value.

### Native contexts and the disappearing neighborhood

The native incoming edge uses a child's test-half mean and the pooled test
halves of its parent's children. At an intermediate node, its incoming and
outgoing means can therefore differ. Native edge projections, covariance,
ridge, branch factors and ranks also need not match across edges. Their
Wald scores and p-values cannot be replaced by a sum based on mean telescoping.

Missing intermediate parents replace nested contrasts and complements. Missing
children replace outgoing sibling partitions. Missing siblings change parent
context and potential empirical-null records. Changing the LCA replaces the
central branch comparison. Ancestors above the connecting path can alter root
reachability, and global branch normalization or empirical correction can change
numeric evidence even when a local path and its neighborhood survive.

The recorded example preserves a 32-to-16-leaf path while the lower endpoint's
children change from 14|2 to 8|8. Its raw T rises from 4.2517 to 78.7587 and an
unavailable sibling gate/boundary becomes active rejection/split. The same 16
observations form one true simulated joint group, so the new split fragments
that group and reduces AMI. A larger score is not evidence of better recovery.

### Connection to earlier tests and eigenvalue roles

Keep the full endpoint contrast, actual native path hypotheses and root-to-node
reachability as separate records. The earlier eigenvalue distinctions still
apply: diffusion weights distance geometry; covariance used for whitening scales
contrasts; native PCA eigenvalues affect ordering, ties and rank. A disappearing
node does not transfer its spectral context to a replacement edge.

Degree-two graph suppression preserves distance when adjacent lengths are
summed. Alternative binary clades generally also change off-path attachments
and receive new NNLS fits, so their fitted length need not equal that sum.
Neither path length nor the number of missing clades supplies a null reference.
Fixed-union Beta and declared-event Ward references keep their existing scope;
selection among paths, endpoint unions and neighborhoods needs its own coverage.

## Evidence

[[subtree-path-statistics-evidence]] records 569,954 saved-data paths, with the
main interpretation restricted to 1,890 internal-endpoint paths, exact Wolfram
checks and 31 passing focused tests. No production behavior changed. The study
does not establish selected-null or final-partition error control.

## Links

- [[subtree-path-statistics-evidence]]
- [[subtree-context-and-eigenvalue-statistics]]
- [[subtree-spectral-alignment-and-crossing-clusters]]
- [[subtree-clustering-production-gates]]
- [[gaussian-root-calibration-and-cut-geometry]]

## Open Questions

- Which endpoint comparisons should remain fixed when a chosen subtree is
  replaced by an alternative topology, and which require a new hypothesis?
- Which reference covers selection of paths and neighborhoods with useful power,
  including native rank, branch, gate and traversal decisions?
- Can preserving declared endpoint and child partitions prevent harmful
  fragmentation while satisfying existing recovery and calibration requirements?
