# Diagnosis of the one-cluster cancer cell-line result

Read-only investigation on 2026-11-03 after the user asked why TBS performed
poorly. Production code, thresholds, the original result, and its labels were
not changed. Diagnostic outputs are in the run folder's `diagnosis/` directory.

## Established failure mechanism

The diffusion/average-linkage hierarchy ends in a one-versus-48,878 root split.
The singleton is `Pool15_19#AAAGATGCATTGAGCT-1` from `ASPC1_PANCREAS`, a cell line
with 290 observations. All 88 of that cell's kernel neighbors have the same
cell-line label. Across the complete symmetrized kernel, 987,220 of 1,030,508
directed stored edges connect matching cell lines (fraction 0.9579935;
weight-weighted fraction 0.9624884). These are descriptive label-agreement
measurements, not an independent validation or a proof of cell-state identity.

TBS tests child-parent mean contrasts in the original 30-PC representation,
while the tree is inferred from diffusion-coordinate distances. At the root,
the alternating basis/test split leaves test sample sizes 1 and 24,439, with
24,440 in their parent. The singleton contrast therefore has sampling variance
factor `1 - 1/24440`, approximately one. The large overall dataset does not
give that singleton test a large effective sample size.

The `guarded_within_child` covariance policy requires at least eight leaves in
each child. The singleton violates that guard. Object-identity checks verified
that the root resolver returns its parent-total covariance, which includes
variation across the entire mixture. This is the implemented fallback, not a
new covariance choice made for the diagnosis.

Both original root-child edges fail the 0.001 threshold. Tree-BH only schedules
a child's family after rejecting its incoming edge. Consequently all 97,754
descendant nodes are marked ancestor-blocked. Pass-through does not override
these closed edge prerequisites. The single-cluster output follows this
hierarchical decision; it is not a test establishing global homogeneity.

## Frozen-tree variance control

Reconstructed the saved linkage and original NNLS fit with the same PCA data,
pair sample, seed, tolerance, iteration cap and four-thread BLAS limit. The
original residual RMSE was reproduced. The saved optimized edges then allowed
replaying only edge annotation without recomputing topology or fitting labels.

The corrected replay explicitly supplies the original 30-dimensional projection
basis, 0.90 energy fraction, minimum dimension 2, covariance guard and alpha.
It reproduces the original root statistics. An initial replay omitted the
explicit basis dimension; its mismatch was caught by an assertion and that
replay was discarded.

| Root child | Original variance multiplier | Original adjusted p | Adjusted p without branch variance |
| --- | ---: | ---: | ---: |
| Singleton `L43204` | 2.011004 | 0.916734 | 0.120680 |
| Other 48,878 cells `N97755` | 1.229813 | 0.724561 | 0.120680 |

The multiplier is `1 + branch_length / mean_branch_length`. Removing it alone
still produces zero significant edges and only two reached Tree-BH tests.
Thus branch variance worsens these root p-values but is not sufficient to
explain or repair the one-cluster outcome. This control did not replace the
original clustering result.

## Geometry and benchmark limits

A diagnostic global height cut at the known count 190 has ARI 0.076053 and NMI
0.663746. This uses the known count, is not TBS, and is not an oracle upper bound
on arbitrary cuts of the same hierarchy. It flags global tree-cut quality for
further investigation; it does not prove that every adaptive cut must be poor.

Among the top 100 merges, 12 have a singleton smaller child and 26 have a smaller
child of at most five cells; the median smaller-child size is 36.5. Within
ASPC1, diffusion-coordinate norm correlates 0.977394 with square root kernel
degree. This is consistent with density-related coordinate amplitudes affecting
the hierarchy, but does not isolate their causal contribution or justify
changing the coordinate operator.

The ranking benchmark contains 122 cases with 30–1,000 observations, verified
from its input fingerprints. Its mean-ARI ranking does not validate a
48,879-cell, 190-cell-line analysis. Passing software tests does not establish
scientific performance or selected-hierarchy calibration in this regime.

## Evidence and unresolved questions

- `tbs-root-variance-diagnostic.json` captures the exact replay, root edge
  lengths, variance factors, statistics, p-values, and reached-test counts.
- `tbs-structure-diagnostic.json` records kernel label agreement, root-cell
  identity, covariance fallback, tree-cut diagnostics and benchmark sizes.
- Owning code is `tree_break_selection/tree/distributions.py`,
  `hierarchy_analysis/statistics/contrast_covariance.py`,
  `hierarchy_analysis/statistics/multiple_testing/tree_bh/correction.py`, and
  `hierarchy_analysis/decomposition/gates/gate_evaluator.py`, relative to the
  production package root.

The demonstrated failure is an unbalanced upper-tree split coupled to a
top-down edge prerequisite. PCA dimension, diffusion coordinate normalization,
NNLS accuracy, alternative tree construction, and a statistically justified
way to examine deeper subtrees have not been isolated in this diagnosis.
No proposed repair has been validated for error control or cancer-line recovery.
