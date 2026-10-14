---
title: Eigenvector Geometry and Calibration Construction Review
type: analysis
status: reviewed
updated: 2026-09-28
sources:
  - reports/mathematical_construction_review_20260928/README.md
  - reports/mathematical_construction_review_20260928/reproduce_findings.py
  - reports/mathematical_construction_review_20260928/findings.json
  - reports/mathematical_construction_review_20260928/rank_invariance_repair.md
  - reports/mathematical_construction_review_20260928/rank_invariance_repair.json
  - reports/eigenspace_reassessment_20260928/README.md
  - tree_break_selection/eigenvalue_groups.py
  - tests/statistics/55_test_eigenspace_invariance.py
  - tree_break_selection/space_separation/diffusion.py
  - tree_break_selection/hierarchy_analysis/decomposition/backends/eigen/projection.py
  - tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_kernel.py
  - tree_break_selection/hierarchy_analysis/statistics/projection/spectral/marchenko_pastur.py
  - benchmarks/diagnostics/path_b/recursive_pvalue_geometry.py
  - README.md
  - raw/inbox/eigensystem-reference-notes-20260928.md
tags:
  - eigenvectors
  - calibration
  - review
---

# Eigenvector Geometry and Calibration Construction Review

## Summary

The 2026-09-28 three-agent review found basis-invariance defects in diffusion
and local PCA truncation, an optional adaptive Wald rank rule with an invalid
fixed-rank reference, an unequal-rank principal-angle diagnostic error, and a
Hamming self-neighbor edge case. Missing separate left eigenvectors of the
current symmetric operators does not explain these findings. The report
separates demonstrated defects from unresolved operator and calibration choices.

The subsequent rank repair now preserves complete numerical groups across
diffusion, local PCA, MP test-dimension metadata and Wald projection. Its 25 new
regressions and full 693-test `make check` passed. The selected reference law,
unequal-rank angle diagnostic and Hamming self-neighbor issue remain open.

The subsequent [[eigenspace-benchmark-calibration-reassessment]] executed the
matched benchmark and independent calibration studies. It confirms remaining
adaptive-rank and complete-pipeline calibration failures; geometric invariance
alone does not resolve the reference law.

## Details

### Review snapshot before the rank repair

- A six-cycle with one retained diffusion component changes matched distances
  by 0.216506 after sample permutation; retaining the whole tied pair changes
  distances only at numerical precision.
- Equal covariance eigenvalues 0.5, 0.5 with the MP floor of one select different
  physical rank-one projectors after feature permutation, with difference
  `sqrt(2)` in Frobenius norm. Larger configured floors do not generally protect
  a boundary inside a larger tied group.
- The optional adaptive contrast-energy rank rule gives p-values 0.002699796
  and 0.011108997 for the same statistic nine under a 45-degree rotation of a
  repeated eigenspace. One million independent Gaussian contrasts produce
  14,720 rejections at nominal alpha 0.01. This isolates a contrast-level issue;
  it is not a full-pipeline false-split estimate. Core/generic defaults disable
  the option, while the pancreas application enables energy fraction 0.90.
- The principal-angle helper truncates both bases to the smaller rank. It can
  label a contained one-dimensional child orthogonal to its two-dimensional
  parent, and change that answer when the parent basis is reordered.
- Hamming `k_neighbors=1` includes self, then diagonal removal leaves all
  diffusion distances zero on a four-row distinct binary example.

At review time the adaptive-cosine tie fix protected only its own module. The
follow-up extracted that rule into `tree_break_selection/eigenvalue_groups.py`
and applied it to the other rank boundaries.

### Verified rank repair

The six-cycle sample-permutation distance discrepancy falls from 0.216506 to
9.71445e-16. The two-dimensional covariance fixture retains both tied axes and
its feature-permutation projector discrepancy becomes zero. Both bases in the
adaptive Wald rotation fixture now return rank two, statistic nine and p-value
0.011108997. MP raw threshold counts are unchanged; the test dimension expands
through the boundary group.

Diffusion removes the leading mode only when its eigenvalue is isolated and
retains a repeated leading group whole. PCA keeps its positive-rank constraints
and returns actual eigenvalues rather than flooring small values into artificial
ties. Missing known tied parent directions are rejected explicitly.

Rank expansion can increase memory and runtime. The later reassessment reran
all 752 benchmark attempts, finding 33 changed successful partitions, and
evaluated 4,200 independent null/signal datasets with the complete pipeline.
At sibling alpha .01, adaptive NNLS produced false partitions on 316/500 iid
Bernoulli null datasets and 263/500 categorical null datasets. Its generic
configuration disables adaptive contrast-energy rank. Separately, held-out
isolated adaptive Wald rates remain above target for distinct and partially
tied eigenvalue groups. The earlier one-million-contrast count remains a
pre-repair snapshot; current results are recorded separately.

### Structural and interpretation questions

Dense full-rank continuous covariance estimated from the same leaf observations
used for whitening produces spectral covariance `(m-1)/m I` at zero ridge.
Positive ridge changes each eigenvalue to `(m-1)/m * c_i/(c_i+ridge)`; its
effect need not be small for ill-conditioned covariance. The live node calculation
gives nearly identical eigenvalues 0.95, 0.95 both before and after a large mean
separation in a 20-by-2 example. This establishes a limitation of that spectral
signal summary, not a complete clustering failure; diagonal/multiple-block
covariance and internal spectral rows require separate analysis.

Cosine, covariance and the implemented normalized diffusion operator are
symmetric. Their algebraic and geometric multiplicities coincide. For the
associated positive-degree Markov operator, similarity to the symmetric operator
also guarantees diagonalizability. Numerical eigenvalue groups are not exact
multiplicity certificates.

For `S=D^(-1/2)WD^(-1/2)` and its orthonormal basis `U`, the Markov pairs are
`R=D^(-1/2)U`, `L=D^(1/2)U`, with `L^T R=I`. Their Euclidean left-left and
right-right angles differ from orthogonality in the appropriate degree-weighted
metrics. The report includes all pairwise angles and block principal angles.
GO sample/feature loadings already implement the correct paired singular-vector
recovery and should not be confused with this square-operator distinction.

The project does not define “shadows”; projected coordinates, projected vectors,
energies, spectral contributions and oblique projections remain distinct.

### Documentation and module reachability

The report records contradictory calibration descriptions within README and
the traversal class, stale annotation ownership, a nonexistent cited angle test,
and an index entry presenting retired legacy code as active. The first PCA axis
is called invariant without a declared transformation group/action. These
conflicts reflect incomplete propagation among documentation surfaces.

The static audit found one diagnostic enrichment helper unused inside the
repository, but documented as opt-in. Selected-Gaussian modules are intentional
tested prototypes; no core module lacked incoming imports in the reviewed scan.
No module was deleted. Historical selected-Gaussian underflow and narrow-interval
errors are now repaired in live source and their regressions passed.

## Evidence

- `reports/mathematical_construction_review_20260928/README.md` contains precise
  source locations, assumptions, derivations, priorities and unresolved questions.
- `reports/mathematical_construction_review_20260928/reproduce_findings.py`
  reproduces the observations, records 17 source hashes and compares the
  vectorized null experiment with 128 actual kernel calls.
- `reports/mathematical_construction_review_20260928/findings.json` contains
  exact results, seeds, versions and all illustrative left/right angle matrices.
- Agents passed 184 existing targeted tests; those tests miss the reproduced
  basis-rotation/permutation failures. Vulture exited zero without findings.
- The original review did not rerun the eigenvalue-block benchmark. The later
  reassessment supplies matched reruns and independent calibration evidence,
  without treating one-cluster benchmark collapse as recovery or calibration.
- `reports/mathematical_construction_review_20260928/rank_invariance_repair.json`
  records matched before/after probes and repaired source hashes. The initial
  regression run had 18 failures and three passing controls; all 25 final
  regression cases pass, and full `make check` passes 693 tests.

## Links

- [[eigensystem-api-reference-notes]]
- [[eigenvalue-block-tree-comparison]]
- [[current-adaptive-diffusion-subspace-tree-pipeline]]
- [[open-mathematical-questions]]
- [[empirical-null-calibration-reference-law-contract]]
- [[diffusion-nnls-calibration-open-points]]
- [[eigenspace-benchmark-calibration-reassessment]]

## Open Questions

- Which diffusion operator, metric, signed spectrum and stationary-space
  treatment should the method promise?
- What should MP identify after same-observation dense covariance whitening?
- How sensitive are retained ranks to the numerical-tie tolerance, and how do
  complete-group ranks affect benchmark results, memory use and runtime?
- What selected reference law includes covariance, PCA, rank, topology,
  branch-time factors and calibration-record selection?
- What is the declared local/global error target, and what accounts for
  per-depth testing, dependent nodes and data-selected calibration weights?
- Which selected-tail correction resolves the independently observed false
  splits while preserving signal recovery on fresh evaluation datasets?
