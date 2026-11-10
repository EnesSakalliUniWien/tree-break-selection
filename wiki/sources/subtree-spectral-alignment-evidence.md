---
title: Subtree Spectral Alignment and Alternative Cluster Evidence
type: source
status: reviewed
updated: 2026-11-09
sources:
  - raw/inbox/choir-cell-lines/subtree-spectral-alignment-20261109.md
  - raw/inbox/choir-cell-lines/subtree-spectral-alignment-20261109-manifest.json
  - raw/inbox/choir-cell-lines/subtree-spectral-alignment-20261109-sources.py.txt
  - reports/subtree_spectral_alignment_20261109/protocol.md
  - reports/subtree_spectral_alignment_20261109/summary.json
  - reports/subtree_spectral_alignment_20261109/parent_child/production_span_summary.json
  - reports/subtree_spectral_alignment_20261109/math/wolfram_evidence.json
  - reports/subtree_spectral_alignment_20261109/verification.json
tags:
  - subtree
  - eigenspaces
  - alternatives
  - validation
---

# Subtree Spectral Alignment and Alternative Cluster Evidence

## Summary

The confirmed follow-up covers subspace equality despite sign/basis differences,
prescribed reflections and matching eigenvalues, for diffusion and covariance
matrices, including left/right and parent/child comparisons. Implemented
diagnostics and finite experiments show how different spectral combinations
can expose crossing partitions that one topology omits. No production preset
or statistical reference changed.

## Key Points

- Added rank-aware principal-angle, containment and projector comparisons,
  explicit reflections, covariance decomposition, reversible diffusion bases,
  and rank-aware parent-row restriction.
- On the same 80 simulated matrices, the alternative panel records 5,120
  candidate outcomes; native replay records 160 fits and 320 child comparisons.
  Every fit/comparison completed, and regenerated input hashes match both panels.
- Across 20 crossed datasets, the family contains both factor splits in all
  20, while the full Ward tree contains only factor A. The full tree nevertheless
  contains all four joint clades in every dataset. Crossing factors cannot both
  be tree splits; omitting a factor grouping is distinct from losing joint groups.
- For equal population amplitudes, the family contains both splits in 11/20.
  The exact tied square has two Ward tie outputs inside its whole rank-two block.
  Selecting axes inside a repeated block is additional information.
- For compatible nested groups, the full tree contains both contrasts and all
  four groups in 20/20; a leading rank-two view retains all four in only 6/20.
- Common-target NNLS residual improvement does not improve the full tree's
  factor-coverage counts in this shortlist. Null datasets also generate 63
  distinct alternatives each, so the count is not a population or null correction.
- Degree transport reduces left/right block discrepancy to roundoff; maximum
  transported gap is 1.65e-15. Production distance reconstruction error is
  1.39e-14. Restriction of the selected unclipped parent Markov-right span loses
  rank in 239/320 comparisons; induced versus
  refitted child operators differ above 1e-10 in 316/320. These are descriptive
  geometric counts, not independent-data rejection rates.
- A separate same-input check applies production clipping and time-three weights.
  Nonpositive modes disappear in 154/320 induced and 160/320 refitted child
  windows. Inherited and child production-coordinate spans still differ in all
  320 comparisons. This adds no independent datasets or native clustering fits.

## Evidence

The primary report retains protocols, source hashes, all attempts, matrices and
native spectral/branch/traversal traces for representative datasets. Exact
Wolfram queries verify degree transport, a child-eigenproblem counterexample,
isospectral covariance/reflection examples, and the scatter identity. Primary
papers are cited in the mathematical report. The final focused suite passes
48 tests after a 30-test baseline; final style/wiki checks are recorded in the
verification artifact. Production algorithm files remain unchanged.

## Links

- [[subtree-spectral-alignment-and-crossing-clusters]]
- [[subtree-clustering-research-evidence]]
- [[subtree-clustering-production-gates]]
- [[subtree-distance-topology-families]]
