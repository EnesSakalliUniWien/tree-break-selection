---
title: Mathematical Glossary
type: concept
status: reviewed
updated: 2026-11-04
sources:
  - manuscript/sections/method/notation.tex
  - manuscript/sections/method/representation.tex
  - manuscript/sections/method/edge_test.tex
  - manuscript/sections/method/sibling_test.tex
  - manuscript/sections/method/decomposition.tex
  - tree_break_selection/tree/feature_space.py
  - tree_break_selection/tree/poset_tree.py
  - tree_break_selection/hierarchy_analysis/statistics/contrast_covariance.py
  - tree_break_selection/hierarchy_analysis/decomposition/backends/eigen/operators.py
  - tree_break_selection/hierarchy_analysis/statistics/projection/projection_dimension_estimation/projection_dimension_estimators.py
  - tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_kernel.py
  - tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_reference_distribution.py
  - tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/pair_testing/wald_statistic/sibling_divergence_test.py
  - tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflation_correction/empirical_null_inflation_estimation.py
  - tree_break_selection/hierarchy_analysis/statistics/multiple_testing/tree_bh/correction.py
  - tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflated_projected_wald_annotation/fdr_annotation.py
  - tree_break_selection/hierarchy_analysis/decomposition/gates/gate_evaluator.py
  - raw/inbox/wolfram-mathematical-construction-checks-20261102.md
  - reports/statistical_construction_audit_20261104/mp_audit.md
  - reports/statistical_construction_audit_20261104/reference_law_audit.md
  - reports/statistical_construction_audit_20261104/knowledge_audit.md
tags:
  - mathematics
  - glossary
  - method
---

# Mathematical Glossary

## Summary

The mathematical objects of Tree-Break Selection, ordered from observations to
the final partition. Symbols follow the manuscript; implementation qualifications
follow the current source. Exact identities, covariance assumptions, reference
laws, and algorithmic decisions are distinguished throughout. See
[[mathematical-construction-wolfram-review]] for the checked derivations and gaps.

## Details

### Observations, feature space, and hierarchy

| Term / symbol | Definition and conditions | Owning source |
| --- | --- | --- |
| Working sample-feature matrix, \(X\in\mathbb R^{n\times p}\) | Rows are observed samples; columns are raw working coordinates. A numeric matrix alone does not declare the statistical model. | `notation.tex`; `feature_space.py` |
| Feature-space contract, \(\mathcal F\) / `FeatureSpace` | Explicit partition of raw coordinates into `FeatureBlock` objects, each declaring family, chart, covariance model, and contrast dimension. | `feature_space.py` |
| Bernoulli block | One binary coordinate, identity chart, model variance \(q(1-q)\). Separate Bernoulli blocks do not encode their cross-feature covariance. | `feature_space.py`; `contrast_covariance.py` |
| Categorical block | One-hot probabilities \(q_1,\ldots,q_K\); `simplex_drop_last` uses \(K-1\) contrast coordinates. Model covariance is \(\operatorname{diag}(q_*)-q_*q_*^\top\), where \(q_*=(q_1,\ldots,q_{K-1})\). | Same sources |
| Continuous block | Identity chart with declared `empirical_gaussian` covariance. The contrast code accepts a full covariance matrix or a diagonal vector for each block. Estimated covariance does not establish an exact finite-sample Gaussian reference law. | Same sources |
| Contrast dimension | Sum of block contrast dimensions; can differ from raw dimension \(p\), especially for categorical blocks. It is distinct from the retained projection dimension. | `feature_space.py` |
| Rooted binary hierarchy, \(\mathcal T=(V,E)\) | Candidate nested groups. For a full binary tree with \(n\) leaves: \(n-1\) internal nodes, \(2n-1\) nodes, \(2n-2\) edges. | `representation.tex` |
| `PosetTree` | Directed representation with parent-to-child edges, leaf labels, descendant-set operations, and node distributions. The class wraps `networkx.DiGraph`; the binary-tree condition is also checked by the split gate. | `poset_tree.py`; `gate_evaluator.py` |
| Descendant leaf set, \(L(u)\) | Observed sample indices below node \(u\); \(n_u=|L(u)|\). For binary children, \(L(u)=L(l(u))\sqcup L(r(u))\). | `representation.tex` |
| Empirical node distribution, \(\theta_u\) | \(n_u^{-1}\sum_{i\in L(u)}X_i\): frequencies for discrete coordinates and means for continuous coordinates. Internal nodes are summaries, not additional observations. | `representation.tex`; `poset_tree.py` |
| Barycentric weight, \(\beta_u\) | \(n_{l(u)}/n_u\), so \(\theta_u=\beta_u\theta_{l(u)}+(1-\beta_u)\theta_{r(u)}\). | `representation.tex` |
| Barycentric balance / leverage | \(b_u=\min(\beta_u,1-\beta_u)\); \(\ell_u^{\mathrm{bar}}=\log[\max(\beta_u,1-\beta_u)/\min(\beta_u,1-\beta_u)]\). These describe split imbalance; they are not calibrated penalties. | `notation.tex` |
| Normalized branch time, \(t_{u,c}\); \(\bar t\) | Tree branch-time observation and mean positive branch time. Branch-time variance factors are modeling choices, distinct from topology and leaf counts. Unit root-to-leaf path sums in the manuscript refer to its linkage-derived convention. | `notation.tex`; `sibling_divergence_test.py` |

### Contrasts, covariance, and projection

For this subsection, \(\Sigma\) denotes a common leaf covariance under fixed
membership and independent leaves; \(\Sigma_0(A,B)\) denotes a contrast covariance.
These are different matrices.

Current edge annotation splits rows within each already selected child: one
half supplies the PCA basis and the other supplies child/parent means and
sample counts. The returned sibling spectral context uses all descendant
leaves. Thus the edge formulas below apply to its testing-half summaries;
they do not describe a full-data edge statistic. Membership selection remains
data-dependent, and continuous contrast covariance can still use all leaves.
Row splitting alone does not establish a calibrated selected reference law.

| Term / symbol | Definition and conditions | Owning source |
| --- | --- | --- |
| Child-parent contrast | \(\Delta_{c,u}=\theta_c-\theta_u\), comparing nested groups. | `edge_test.tex`; `contrast_covariance.py` |
| Sibling contrast | \(\Delta_{l,r}=\theta_l-\theta_r\), comparing disjoint children. Disjointness alone does not imply independence for correlated leaves or data-selected membership. | `sibling_test.tex`; `contrast_covariance.py` |
| Nested covariance scale | \(\operatorname{Cov}(\Delta_{c,u})=(1/n_c-1/n_u)\Sigma\) under the fixed-membership common-covariance assumptions above. | `contrast_covariance.py`; Wolfram capture |
| Disjoint two-sample covariance scale | \(\operatorname{Cov}(\Delta_{l,r})=(1/n_l+1/n_r)\Sigma\) under the same assumptions. | Same sources |
| Branch-time variance multiplier | Edge: \(1+t_{u,c}/\bar t\) only under `normalized_branch_length`; default edge policy is `none`. Adaptive projected sibling: \(1+(t_{u,l}+t_{u,r})/(2\bar t)\), regardless of the edge policy. Fixed-subspace sibling gates omit this factor. | `contrast_covariance.py`; `sibling_divergence_test.py`; [[projected-wald-statistic]] |
| Covariance whitening, \(z\) | For model contrast covariance \(V=CC^\top\), the code computes \(z=C^{-1}\Delta\) blockwise using Cholesky factors. This yields identity covariance only when \(V\) matches the actual covariance, including the relevant dependencies. | `ContrastCovariance.whitened_vector`; `contrast_covariance.py` |
| Ridge | Positive diagonal regularization, commonly \(10^{-12}\) here. It makes numerical factorization possible; it does not establish a reference law at probability boundaries. Exact identities below omit ridge. | `contrast_covariance.py` |
| Parent spectral covariance | Covariance of centered, null-whitened descendant rows, used to choose directions. It is distinct from the contrast covariance used to standardize a test. | `operators.py`; [[projected-wald-statistic]] |
| Primal / dual spectral matrices | For centered active tangent rows \(Y\in\mathbb R^{m\times d}\), \(C=Y^\top Y/m\), \(G=YY^\top/m\). Their positive eigenvalues agree. Current code preserves tangent scale; it does not perform an additional Pearson standardization. | `operators.py` |
| Projection basis, \(P_u\) | Orthonormal rows, \(P_uP_u^\top=I_k\), chosen from available parent PCA directions. PCA orders candidate directions. The nominal or contrast-energy rank, extended through numerical eigenvalue groups, determines retained rows; the statistic does not divide by PCA eigenvalues. | `projected_wald_kernel.py`; `projected_wald_reference_distribution.py` |
| MP upper edge | \(\lambda_+=(1+\sqrt{d/m_{\mathrm{MP}}})^2\) on the unit null-covariance scale. A working spectral threshold; its applicability needs local spectral assumptions. | `projection_dimension_estimators.py` |
| Raw MP signal count | Number of eigenvalues strictly above \(\lambda_+\). Can be zero. | Same source |
| Test projection dimension, \(k\) | The spectral field stores a nominal rank after the floor, available positive directions and numerical eigenvalue-group extension. It need not equal the raw MP count. The final kernel rank can differ when contrast-energy selection is enabled. | Same source; `projected_wald_kernel.py` |
| Contrast-energy rank | Optional rank target retaining a declared fraction of the tested projected contrast energy, extended through eigenvalue groups. Its dependence on the tested contrast requires separate null-law justification. | `projected_wald_kernel.py` |
| Projected-Wald statistic | Edge: \(T^{\mathrm{edge}}_{u\to c}=\|P_uz^{\mathrm{edge}}_{u\to c}\|^2\); sibling: \(T_u^{\mathrm{raw}}=\|P_uz_u^{\mathrm{sib}}\|^2\). The kernel reports scale \(a_u=1\) and df \(\nu_u=k\). | `projected_wald_kernel.py` |
| Fixed-subspace chi-square law | If \(z\sim N(0,I_d)\) and \(P\) is fixed, or independent with conditional isotropy, then \(\|Pz\|^2\sim\chi_k^2\). Orthonormal rows alone do not imply this law. | `projected_wald_reference_distribution.py`; Wolfram capture |
| Zero-dimensional context | The kernel returns statistic 0, df 0, p-value 1 only for a numerically zero contrast; a nonzero contrast at requested rank 0 is rejected. This is an explicit convention, not a usual positive-df chi-square test. | `projected_wald_kernel.py` |
| Missing edge basis | The tree-wide edge wrapper emits statistic 0, df 0, p-value 1 before calling the kernel when its parent basis rank is zero. This means no edge test, not zero observed contrast or evidence of equality. | `child_parent_divergence_annotation/tree_testing.py` |

The manuscript uses \(u\to c\) for a parent-to-child edge annotation and the
contrast numerator \(\theta_c-\theta_u\). Some derivation text labels that same
contrast \(c\to u\); read its numerator to disambiguate the notation.

### Empirical null, multiplicity, and final clusters

| Term / symbol | Definition and conditions | Owning source |
| --- | --- | --- |
| Reference scale and df, \(a_u,\nu_u\) | Current raw orthonormal reference uses \(a_u=1\) and df equal to the final retained basis row count, after any adaptive replacement and tie extension. These are reference parameters, not a proof of actual selected-null moments. | `projected_wald_kernel.py` |
| Calibration context, \(x_q\) | \((\log s_q,\log n_q)\), where \(s_q\) is resolved sibling projection dimension. Local axes depend on family and positive bandwidth; Bernoulli prediction omits the parent-size axis. | `empirical_null_inflation_estimation.py` |
| Null / context weights, \(\pi_q,w_q(u)\) | Calibration roles use edge status; stopped-edge records can receive weight one despite a nonzero generating mean contrast. Local Gaussian-kernel weights use the same feature family. Dependency groups describe support; they establish neither null status nor independence. | Same source |
| Empirical inflation, \(\hat c\) | Local estimate \(\max[1,\sum_qw_qT_q/\sum_qw_qa_q\nu_q]\), with family baseline when no axes are active. The current working tree also caps certain larger-parent extrapolations at the family baseline. | Same source; see review for working-tree qualification |
| Effective sample size, \(n_{\mathrm{eff}}\) | Weight concentration measure \((\sum w)^2/\sum w^2\). Support diagnostics use dependency-group weights where applicable; this is not automatically a count of independent observations. | Same source |
| Inflation-adjusted sibling p-value | Active empirical plug-in tail \(\Pr(\chi^2_{\nu_u}\ge T_u^{\mathrm{raw}}/(a_u\hat c_u))\). Matching a scale or checking support does not prove selected-tail validity. | Same source |
| Unsupported calibration context | An explicit undefined/skipped state when required internal/family support is missing; it is not evidence for equality of siblings. Additional sparse-support thresholds apply when enforcement is enabled. | Same source; `gate_evaluator.py` |
| Tree-BH edge family | Outgoing edges of a reached parent; local BH level is base alpha multiplied by ancestor rejection proportions. Only rejected child edges open descendant families. Stored adjusted p-values are within-family adjustments, evaluated against each family's level. | `tree_bh/correction.py` |
| Traversal-aligned sibling BH | BH separately on eligible parents at each edge-reachable depth frontier. Frontier expansion follows edge prerequisites and available sibling results, not sibling rejection alone. Algorithmic adjustment does not establish selected-hierarchy FDR. | `fdr_annotation.py` |
| Edge indicator, \(S^{\mathrm{edge}}(c)\) | Supported child-parent edge rejection after Tree-BH; used as an annotation by traversal. | `notation.tex`; `gate_evaluator.py` |
| Sibling indicator, \(S^{\mathrm{sib}}(u)\) | Available, unskipped sibling rejection after the configured sibling correction. | Same sources |
| Tree break, \(B_{\mathrm{TBS}}(u)\) | Binary parent, at least one significant child edge, and an open sibling gate. | `decomposition.tex`; `gate_evaluator.py` |
| Pass-through | With pass-through enabled, descend despite a closed current sibling gate only when binary/edge prerequisites hold, a deeper full split is reachable, and any configured pass-through support allows it. A significant sibling somewhere below is insufficient. | `gate_evaluator.py` |
| Terminal cluster nodes, \(\mathcal C_{\mathrm{TBS}}\) | Boundary nodes returned by traversal. Their descendant sets are disjoint and cover the input leaves for a valid rooted hierarchy. `TreeDecomposition` converts annotations into these clusters. | `decomposition.tex`; [[tree-decomposition]] |
| Selected region, \(\mathcal S_u\) | Data event selecting hierarchy membership, projection/rank, edge opening and focal testing/traversal context. The needed law is conditional on this event, not merely on an unselected fixed node pair. | [[method-proof-web]]; [[mathematical-construction-wolfram-review]] |

### Exact connection between edge and sibling contrasts

Write \(l=l(u)\), \(r=r(u)\). The barycenter identity gives
\[
\theta_l-\theta_u=(1-\beta_u)(\theta_l-\theta_r),\qquad
\theta_r-\theta_u=-\beta_u(\theta_l-\theta_r).
\]
With a shared covariance baseline, compatible whitening, no branch-time factors
and no ridge, nested standardization cancels the barycentric coefficients:
\[
z^{\mathrm{edge}}_{u\to l}=z^{\mathrm{sib}}_u,\qquad
z^{\mathrm{edge}}_{u\to r}=-z^{\mathrm{sib}}_u.
\]
The corresponding split action is
\[
n_l\|\theta_l-\theta_u\|^2+n_r\|\theta_r-\theta_u\|^2
=\frac{n_ln_r}{n_u}\|\theta_l-\theta_r\|^2.
\]
Wolfram verified the scalar coordinate identities; summing coordinates gives
the Euclidean vector identity. They explain why edge and sibling statistics
can depend on the same contrast. They do not assert that the implemented
statistics coincide when branch factors, covariance estimates, ranks, or bases differ.

## Evidence

- The frontmatter lists inspected manuscript and implementation sources. Table
  basenames identify those owners; links connect existing synthesis.
- `raw/inbox/wolfram-mathematical-construction-checks-20261102.md` preserves exact
  evaluator inputs and returned outputs, including assumptions and a corrected
  initially unevaluated expression.
- On 2026-11-02, 48 focused existing tests passed across wiki lookup, covariance,
  the projected reference helper, and gate evaluation. This checks implementation
  contracts, not the mathematical validity of selected p-values.

## Links

- [[mathematical-construction-wolfram-review]]
- [[wolfram-mathematical-construction-checks]]
- [[projected-wald-statistic]]
- [[tree-break-selection]]
- [[tree-decomposition]]
- [[poset-tree]]
- [[top-down-traversal]]
- [[method-proof-web]]
- [[open-mathematical-questions]]

## Open Questions

- Under which feature dependencies does the declared block covariance actually
  yield an isotropic contrast? How should covariance estimation be accounted for?
- What conditional law handles same-data topology, PCA, contrast-energy rank,
  edge selection, and supported empirical-null estimation together?
- Which guarantees apply to the final partition after both multiplicity stages
  and pass-through? See the review for current manuscript/code discrepancies.
