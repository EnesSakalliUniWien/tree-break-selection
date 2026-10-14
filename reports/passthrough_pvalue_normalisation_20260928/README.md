# Passthrough, p-values, and normalisation

Date: 2026-09-28. Two Luna agents independently audit normalisation and the
p-value/clustering relationship. This report adds diagnostics; it does not
change production methods or defaults.

## What the matched comparison establishes

Across the existing 10,000 matched Gaussian datasets, all 30,000 node-level
raw p-values and independently reconstructed BH-adjusted p-values are identical
between passthrough disabled and enabled. Enabling passthrough changes the
partition count in 3,971 datasets and reaches 7,942 additional internal-node
occurrences. These are repeated analyses of existing datasets, not new trials.

The equality of raw p-values in that benchmark is **by construction**: it injects
the same known-covariance tests into both traversals. The production code also
separates the stages: `apply_traversal_aligned_sibling_bh_results` computes BH
over edge-reachable frontiers before `GateEvaluator.decision` applies
`passthrough`. The BH frontier continues through eligible parents regardless of
their sibling rejection. The traversal boolean is absent from that BH function.

| Pattern | Changed partition counts, identity / AR(0.6), each out of 1,000 |
| --- | --- |
| Global null | 8 / 8 |
| Root signal only | 0 / 0 |
| Signal only below A | 989 / 988 |
| Separate signals below A and B | 991 / 987 |
| Strong root, weak A | 0 / 0 |

For signal only below A, A is reached and split in 11 / 12 datasets without
passthrough and all 1,000 / 1,000 with it. Under the global null, any reached
false split increases from 12 / 15 to 20 / 23. A `pass_through` decision is not
itself a significant sibling split. Descending nonetheless refines the final
partition, and permits descendant rejections to affect clustering.

## Why the displayed p-value distribution can change

A distribution restricted to reached nodes changes its membership. For
example, at null node B in the A-only identity case, the mean reached raw p-value
is 0.4214 among 11 reached nodes without passthrough and 0.4937 among 1,000 with
it. No individual node p-value changed. These selected summaries are
descriptive; overlapping nodes are not independent calibration replicates.

`node_reach_pvalue_summary.csv` keeps pattern, covariance, node and traversal
setting separate. `paired_partition_changes.csv` reports dataset-level paired
changes. `verification.json` records input and script hashes and verified counts.

## Production feature-matrix probe

The separate real-data probe uses `data/feature_matrices/feature_matrix.tsv`
(626 observations, 456 binary features), the endotype app's Hamming/average
linkage path, edge alpha .001, sibling alpha .01, and
`projected_wald_inflation` with support contract version 2. It computes one
gate annotation bundle and reuses it for both traversal settings. This is the
plain Hamming tree path, not an executed diffusion-tree comparison.

All 1,251 annotation rows remain unchanged, including the 345 available
empirical-inflation plug-in sibling p-values and their BH corrections. On this
path, `Sibling_Divergence_P_Value` is the inflation-adjusted pre-BH value,
not the original unadjusted nominal chi-square p-value. Availability of this
empirical rule does not establish selected-tail validity.

| Quantity | Passthrough disabled | Passthrough enabled |
| --- | ---: | ---: |
| Final clusters | 1 | 109 |
| Reached internal nodes | 1 | 182 |
| Reached leaves | 0 | 35 |
| Significant split decisions | 0 | 32 |
| Pass-through decisions | 0 | 76 |

The root's inflation-adjusted pre-BH p-value is 0.02693. The median of those
plug-in p-values among reached nodes changes from 0.02693 (one node) to
0.44612 (182 nodes).
That is a change in the selected set, with no changed node p-values. These
real data have no declared population-null or clustering truth in this probe:
109 clusters cannot be interpreted as improved recovery or as 108 false splits.

`production_probe.json` preserves configuration and source/input hashes;
`production_node_pvalues.csv` records per-node values, equality, reach and
decisions. `pvalue_probe.py` reproduces the production annotation/traversal
comparison. This is one input and one configuration, not a calibration study
across tree inference methods.

## Factors that can actually change p-values

The projected-Wald kernel computes `T = ||R z||²` and nominal
`p = chi2.sf(T, k)` for a supplied orthonormal projection with retained rank k.
PCA eigenvalues identify the selected subspace; they do not weight this
reference law. The empirical-inflation path additionally changes the reference
scale before sibling BH. Thus the following must be held fixed in a passthrough
comparison:

| Factor | Route to a changed p-value and clustering |
| --- | --- |
| Data, tree and child memberships | Change the mean contrast, sample sizes, covariance and selected node hypothesis. |
| Covariance, whitening and ridge | Change the standardized contrast z and therefore T. |
| Projection subspace | Changes the retained contrast energy; basis rotations within a complete retained orthonormal subspace leave T unchanged. |
| Retained rank and complete eigenvalue groups | Change both T and degrees of freedom; their net effect on p need not be monotone. |
| Calibration reference and fitted inflation | Change the survival probability for the same statistic. |
| Edge eligibility and BH family | Can change adjusted p-values for unchanged raw p-values. The passthrough boolean alone does not rebuild these families in the inspected path. |
| Traversal reach | Changes which already annotated rejections affect the partition, and which p-values appear in a reached-only summary. |

The independent-projection control already found, at fixed rank one and a fixed
global-null root, 33 / 37 rejections using the evaluation-data PCA basis versus
11 / 10 using an independent-data basis, per 1,000 datasets at alpha .01.
Passthrough cannot resolve this upstream selection effect. Nor does a small
p-value measure eigenvalue-cluster distance or final-clustering quality.

The cached global-null root comparisons quantify upstream changes separately
from passthrough (raw rejection counts at alpha .01; each cell is identity /
AR(0.6), each out of 1,000 matched datasets):

| Upstream configuration | Rejections | Retained df |
| --- | --- | --- |
| Known covariance, fixed tree, full coordinates | 12 / 15 | 8 |
| Known covariance, fixed tree, same-data PCA/MP | 33 / 37 | 1 |
| Parent pooled sample covariance, full coordinates | 11 / 10 | 8 |
| Parent pooled sample covariance, same-data PCA/MP | 11 / 10 | 8 |
| Known covariance, same-data average-linkage root | 765 / 762 | 8 |

`upstream_null_root_factor_summary.csv` also preserves paired gained/lost
rejections, median p-values and median paired differences. The PCA comparison
above changes both basis and rank; only the separate independent-projection
control isolates basis reuse at fixed rank. Selecting a different root changes
child membership and sample sizes together with selection. These comparisons
do not identify independent additive contributions of every upstream factor.

## Normalisation belongs to specific stages

The GO spectral-block workflow uses feature weighting and row-normalized cosine
geometry to construct diffusion distances, but supplies the original typed
feature matrix to TBS. Changing the tree coordinates can change memberships and
therefore tests; it does not directly rescale the feature observations passed
to the node tests.

- Preserve Bernoulli 0/1 and categorical simplex encoding for their respective
  covariance models. Global z-scoring is not a replacement for these models.
- Continuous node covariance uses scatter divided by N−1. The default
  `parent_total` covariance includes between-child mean variation; it is not
  within-child residual covariance. The PCA covariance divisor N serves a
  different construction and must not be silently replaced just to match it.
- Project the Wald contrast in the same whitening coordinates used to define
  the projection basis. Applying a second unrelated standardization changes
  the statistic and its reference assumptions.
- The custom diffusion routine embeds eigenvectors U of the symmetric operator
  `S = D^(-1/2) W D^(-1/2)` as `U lambda^t`. For positive degrees, the
  corresponding random-walk matrix `P = D^(-1) W` has paired right vectors
  `D^(-1/2) U` and left vectors `D^(1/2) U`, whose pairing is the identity.
  These coordinate choices produce different sample distances. The code uses
  pydiffmap's adaptive kernel followed by its own symmetric normalization;
  it does not apply every step of `pydiffmap.DiffusionMap`.
- Passthrough consumes completed gates. It requires no additional data
  normalization or refitting merely because traversal continues below a node.

These distinctions do not demonstrate a new normalization defect causing the
passthrough effect. Two reference-law sensitivities remain: separate Bernoulli
variance blocks need not remove cross-feature correlation, and nonzero ridge
applied before versus after sample-size covariance scaling can give different
whitening coordinates. Together with same-data tree/basis/rank selection,
these require controlled calibration; they cannot be repaired by rescaling
observed p-values after traversal. The normalization agent's report gives the
exact source locations and applicability limits.

## Evidence and scope

- `audit_observed_pvalues.py`: independent BH reconstruction and paired audit
  of the cached three-node Gaussian benchmark.
- `normalisation_audit.md`: normalisation/data-flow agent's source audit.
- `pvalue_audit.md`: p-value/clustering agent's source audit and production probe.
- Prior evidence: `../passthrough_benchmark_20260928/README.md` and
  `../independent_projection_vandermonde_20260928/README.md`.

The Gaussian study fixes the tree and known covariance and deliberately opens
edge prerequisites. It does not establish general error control for inferred
trees, diffusion geometry, empirical inflation or categorical observations.
The targeted baseline passed 22 tests. The production agent's separate baseline
passed six tests (overlapping coverage, not 28 distinct tests). Root inspection
of the production CSV independently verified all before/after annotation fields
named for p-values, traversal counts, and the recorded input plus seven current
source hashes. These source hashes record the post-run state, not a complete
pre-run source snapshot. Diagnostic Ruff and final wiki checks are recorded in
the wiki log. No production source was modified by this audit.
