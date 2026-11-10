---
title: Subtree Spectral Alignment and Crossing Cluster Partitions
type: analysis
status: reviewed
updated: 2026-11-09
sources:
  - raw/inbox/choir-cell-lines/subtree-spectral-alignment-20261109.md
  - reports/subtree_spectral_alignment_20261109/math/README.md
  - reports/subtree_spectral_alignment_20261109/alternatives/README.md
  - reports/subtree_spectral_alignment_20261109/summary.json
  - reports/subtree_spectral_alignment_20261109/parent_child/production_span_summary.json
  - benchmarks/diagnostics/calibration/selected/hierarchy/subtree_spectral_alignment.py
  - tests/statistics/77_test_subtree_spectral_alignment.py
  - tests/statistics/78_test_subtree_alignment_native_replay.py
tags:
  - subtree
  - eigenspaces
  - topology
  - interpretation
---

# Subtree Spectral Alignment and Crossing Cluster Partitions

## Summary

Different spectral views can expose incompatible groupings of the same leaves.
A single rooted tree can omit a grouping even when its finer joint groups are
recovered. Projection can separately remove a group distinction. These effects
must be distinguished from arbitrary eigenvector signs, basis changes inside
repeated eigenblocks, and changes of coordinates between operators or nodes.

## Details

### Which equality is being assessed

Compare ranks and projectors in shared coordinates. Singular values of Q_A^T Q_B
give principal-angle cosines, but zero angles alone do not distinguish equality
from containment of a lower-dimensional span. Equal covariance eigenvalues do
not identify a common eigenspace: diag(4,1) and diag(1,4) have orthogonal leading
directions. A prescribed reflection answers a separate mapped-comparison
question; arbitrary fitted orthogonal alignment makes equal-rank span mapping
automatic. Numerical eigenvalue groups are not statistical population equality.

For symmetric nonnegative W and positive D, S=D^(-1/2) W D^(-1/2) and
P=D^(-1) W have corresponding eigenvectors U, R=D^(-1/2)U and L=D^(1/2)U.
Then L^T R=I. Known degree transport recovers a shared symmetric representation;
raw Euclidean distances can still differ. Real covariance matrices are symmetric,
so their operator left/right eigenvectors can coincide. Rectangular-data left
and right singular vectors instead live in observation and feature coordinates.

### Parent and child relationships

With common features, sample covariance satisfies
(n-1)C_P=(n_L-1)C_L+(n_R-1)C_R+n_L n_R/n (mu_L-mu_R)(mu_L-mu_R)^T.
Between-child means may create a parent direction absent within either child.
Native Wald/MP contexts use node-specific whitening; raw covariance comparisons
do not silently substitute for comparisons of those native contexts.

Restricted parent diffusion coordinates, an induced parent-weight subgraph and
a newly fitted child kernel are distinct constructions. Restriction can reduce
rank. Induced graphs change normalization; refitting can also change neighbors,
bandwidths, scale and weights. Match child leaf IDs/order and recompute rank
before comparing parent and child spans.

Keep unweighted selected eigenvector windows separate from the production spans
after nonnegative eigenvalue clipping and diffusion-time weights. The saved
same-input check finds zero-weight child modes and confirms that inherited and
refitted spans still differ after applying production rules. These differences
are geometric descriptions, not tests of common population structure.

### Alternatives represented by one tree

For four groups, AB|CD and AC|BD have four nonempty intersections. Both cannot
be splits of one tree. The tree may still represent four joint groups, so the
missing alternative does not imply erased data or failed fine partitioning.
In the measured compatible example, leading-rank truncation independently
reduces four-group clade recovery from 20/20 to 6/20.

An exact repeated eigenblock is invariant to its internal basis. The square's
whole tied rank-two geometry permits both tied Ward outputs, but neither axis
is uniquely chosen by its spectrum. At finite sample size, a population tie can
appear as separate sample eigenvalues; candidate recovery then depends on the
numerical grouping and reconstruction rule. Search maxima evaluated using truth
are descriptive oracle bounds, not fitted methods.

### Consequences for further production work

Retain alternative partitions, their compatibility and their subspace/metric
provenance as diagnostics. Evaluate factor partitions separately from their
joint refinement. Measure loss from rank truncation and distinguish inherited
from rebuilt child geometry. Any production selection among these alternatives
requires a reference covering the complete search, covariance/rank decisions,
branch fits and gates. Neither symmetry nor an alternative count supplies it.

## Evidence

[[subtree-spectral-alignment-evidence]] records the 80 paired matrices, 5,120
candidate rows, 160 native fits, 320 child comparisons and exact checks. Both
native configurations still reject the common-mean Gaussian root in 20/20
datasets in this panel. Geometric recovery therefore does not establish a
calibration repair or a promoted production improvement.

## Links

- [[subtree-context-and-eigenvalue-statistics]]
- [[subtree-spectral-alignment-evidence]]
- [[subtree-clustering-production-gates]]
- [[subtree-distance-topology-families]]
- [[subtree-spectrum-topology]]
- [[subspace-distance-ties]]

## Open Questions

- Which grouping distinctions should be retained when the desired output is a
  single hierarchy versus several compatible or crossing partitions?
- Can a practical candidate policy preserve weaker group distinctions and
  whole-block invariance while meeting the existing runtime requirement?
- Which selective reference covers the chosen subspace unions, topology and
  downstream native gates without destroying useful power?
- How do these mechanisms behave with imbalance, higher dimension, non-Gaussian
  data and native whitening transported into common coordinates?
