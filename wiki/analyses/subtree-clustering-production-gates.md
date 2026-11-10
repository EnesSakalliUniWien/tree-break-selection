---
title: Subtree Clustering Equations and Production Gates
type: analysis
status: reviewed
updated: 2026-11-09
sources:
  - reports/subtree_path_statistics_20261109/README.md
  - reports/subtree_context_statistics_20261109/README.md
  - raw/inbox/choir-cell-lines/subtree-clustering-research-20261109.md
  - reports/subtree_clustering_research_20261109/calibration/factorial_cell_0_audit.md
  - reports/subtree_clustering_research_20261109/calibration/wolfram_queries.md
  - reports/subtree_clustering_research_20261109/recovery_comparison.json
  - reports/subtree_clustering_research_20261109/promotion_gate.json
  - benchmarks/diagnostics/calibration/selected/hierarchy/subtree_family_reference.py
  - benchmarks/diagnostics/calibration/selected/hierarchy/subtree_statistic_attribution.py
  - benchmarks/diagnostics/calibration/selected/hierarchy/subtree_topology_candidate.py
  - benchmarks/validation/sweeps/subtree_research_validation.py
tags:
  - subtree
  - equations
  - calibration
  - production
---

# Subtree Clustering Equations and Production Gates

## Summary

Topology, fitted distances, covariance/rank, eligibility and empirical correction
must be examined together. The implemented diagnostics identify calibration and
runtime failures before production promotion. Better distance fit alone does not
establish improved clustering or a valid selected reference.

## Details

### Distances and alternative topologies

For a candidate topology T and a common declared pair-distance vector d, fit
nonnegative branches by minimizing ||A_T b-d||^2. Comparing this loss across
candidate topologies answers a distance-representation question. Ward minimizes
another objective; the native statistic, empirical reference and traversal can
see different information. Zero-length internal edges can make distinct binary
refinements one unresolved metric. The compressed split-Gram implementation
preserves all pairs and the dense NNLS objective.

The bounded candidate keeps the existing topology, then adds one deterministic
Ward topology from each whole spectral block and leading span. It does not
enumerate every near-tie or every geometric topology for large subtrees. No
eigenspace count, observed candidate count, residual or bootstrap frequency is
used as a replacement null correction.

### Branch factors and the native reference

The native factor is g=1+(bL+bR)/(2 mean_positive_branch). Positive global branch
rescaling leaves g unchanged. With fixed projection and unregularized covariance,
the score scales inversely with g. A fixed ridge added after variance scaling
breaks the exact identity T_unit=g T_native; rank selection can also change.
The implemented unit-branch control recomputes the statistic.

Replacing linkage lengths by NNLS can change edge eligibility before a sibling
p-value exists. Energy truncation can change the collection used for empirical
correction even when a particular root's raw score, degrees of freedom and
edge p-values remain unchanged. The null factorial verifies both mechanisms.
Disabling energy truncation retains the native data-dependent edge/MP projection;
it does not produce the fixed-Gram-rank Beta reference.

### Conditional node families

For each true-null node, require independent Gaussian rows with a common mean
and covariance, independent of outside rows. Condition on its mean, centered
feature Gram, outside observations and one shared labeled Ward subtree/ancestry
event. The reversible orthogonal kernel with independent random anchor supplies
a conservative conditional rank, subject to numerical event evaluation.
Conditioning fixes the hypothesis family. Conditional node validity followed by
Holm controls false rejection of homogeneous nodes under the stated model,
including partial-null families with heterogeneous outside means.

This is stronger family scope than the earlier complete-null search maximum,
but does not establish equal-average-only inference, full diffusion/search/gate
selection, traversal or final-partition error control. With m hypotheses and B
transitions, the grid 1/(B+1) must reach alpha/m; hence B>=ceil(m/alpha)-1 is
necessary for a first Holm rejection. Unresolved calls remain unresolved and
enter Holm conservatively with p=1. Minimum resolution is not demonstrated power.

### Acceptance and observed results

AMI is primary; matched F1, pair split/merge shares, ARI and statistical profiles
are supplementary. There is no separately adopted 0.02 rare-group F1 margin.
Evidence validators enforce exact coverage, valid status, paired seed and data
identity, distinct dataset keys, and explicit phase-2 provenance. Development or
unresolved rows cannot produce a frozen-validation pass.

The development configurations show 47/50 selected-root rejections on a Gaussian
null and a mean Markov/Ward AMI decrease of 0.05550. The alternative topology
candidate also exceeds the fivefold runtime gate on two of six inputs. Neither
configuration comparison nor conditional Ward power justifies a new preset.
Frozen panels remain unrun under the preregistered failed-prerequisite rule.

## Evidence

[[subtree-paths-and-missing-neighborhoods]] extends this comparison to surviving
endpoint pairs. Full-data mean contrasts telescope, while native edge contexts
and squared scores do not compose into an endpoint test. Six unchanged internal
paths change a boundary/split decision with altered neighborhoods in the saved
panel. Path, neighborhood and root reachability need separate provenance.

The later [[subtree-context-and-eigenvalue-statistics]] trace identifies missing
clades and changed incoming/outgoing hypotheses separately from unavailable gates
and unvisited nodes. Even matched parent and child groups can change raw scores
through the refitted branch context. Eigenvalue magnitude has distinct distance,
whitening and rank roles; none replaces the complete selection reference.

The later [[subtree-spectral-alignment-and-crossing-clusters]] study adds
parent/child degree and feature-coordinate comparisons and demonstrates missing
crossing factors separately from loss due to rank truncation. Its requirements
extend candidate evaluation; it does not replace the calibration or runtime gates.

[[subtree-clustering-research-evidence]] links 4,000 native development fits,
300 restricted Gaussian/attribution records, 10,000 perturbations, statistical
profiles, Wolfram checks, isolated timing, source archives and verification.
The primary report distinguishes implemented diagnostics from unestablished
production inference and unrun frozen experiments.

## Links

- [[subtree-clustering-research-evidence]]
- [[subtree-context-and-eigenvalue-statistics]]
- [[subtree-paths-and-missing-neighborhoods]]
- [[subtree-spectral-alignment-and-crossing-clusters]]
- [[subtree-distance-topology-families]]
- [[gaussian-root-calibration-and-cut-geometry]]
- [[diffusion-topology-branch-rank-isolation]]
- [[root-validation-generalisation-evidence]]

## Open Questions

- Which joint reference preserves useful information while including diffusion,
  alternative search, branch fits, empirical eligibility and traversal?
- Can prespecified or proved candidate pruning meet the fivefold runtime cap
  while preserving measured coverage against the small exact oracle?
- Can the conditional family test attain useful power across size, overlap,
  effective rank and partial-null structure at practical budgets?
- Which distributional model and distance uncertainty region justify confidence
  statements about population topology rather than algorithm-output stability?
