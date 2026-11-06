# Covariance, fitted inflation, multiplicity and traversal: answers 5–6

Date: 2026-11-02. Scope: source inspection and mathematical synthesis; no production changes or threshold tuning. Dataset-level evidence is the existing 50-replicate attribution experiment, not a new calibration experiment.

## Q5. What conditional law accounts for covariance estimation and fitted inflation?

**Answer:** no single reference law presently covers the production covariance regimes, selected diffusion partition, selected projection, empirical-null weights and fitted scale. An ordinary F substitution repairs one fixed-partition covariance problem; it does not remove the tree-selection event. Fitted inflation adds a further random, selected quantity whose joint distribution is not supplied by the current chi-square plug-in rule.

### Established independent-partition laws

For a partition independent of iid Gaussian test observations with covariance I, write n=n_L+n_R, p=24, nu=n-2, s=1/n_L+1/n_R and delta=mean_L-mean_R. Before the numerical ridge:

| Covariance actually used | Law under the stated independent global Gaussian null | Production root condition in n=80,p=24 panel |
| --- | --- | --- |
| Known I | Q=delta' delta/s has chi-square(p) law | All roots |
| Ordinary dense pooled within-child S | ((nu-p+1)/(p nu)) delta' S^-1 delta/s has F(p,nu-p+1) law | Both children have at least 25 leaves |
| Diagonal pooled within-child variances | Q is the sum of p independent F(1,nu) variables under covariance I | Guard passes but at least one child has fewer than 25 leaves |
| Dense parent-total S_total | Q/(n-1) has Beta(p/2,(n-p-1)/2) law | Guard fails; in this panel one root has a 7-leaf child |

The diagonal statement uses independent features and the same pooled residual df in every coordinate. It does not extend to correlated unknown covariance. These laws describe different random matrices; no common F correction applies to the table. A child stores diagonal covariance when child_count-1<p; the guarded resolver pools variances when either child is diagonal. Production adds 1e-12 to contrast covariance, slightly perturbing the unregularized laws. In the recorded dense roots, the largest observed statistic perturbation was 3.2183e-8.

Evidence: `tree_break_selection/tree/distributions.py`, `_finalize_continuous_covariances`, `_pooled_residual_degrees_of_freedom`, `resolve_node_continuous_covariance_by_block`; `tree_break_selection/hierarchy_analysis/statistics/contrast_covariance.py`, `_block_contrast_covariance`; `raw/assets/diffusion-attribution-20261102/base_runs.csv`, `base_covariance_strata.csv`; `tests/core/test_diffusion_attribution.py`; `raw/inbox/diffusion-attribution-wolfram-20261102.md`.

Known I isolates the selection problem: selected observations reject 50/50 roots at nominal .01; fresh independent observations on the same trees reject 0/50. Production covariance gives 8/50 independent rejections: 7/25 dense pooled, 1/24 diagonal pooled, 0/1 parent-total. Exact ordinary F gives 0/25 independent dense rejections, but 25/25 selected dense rejections. These finite counts illustrate separate causes; they do not establish exact tail probabilities.

### What the selected law must additionally contain

For known fixed isotropic covariance, a fixed root contrast v and its orthogonal residual can be separated. Conditional on the residual, contrast direction and the replayed selection event, the contrast radius has a truncated chi law. That is a conditional selection correction, rather than covariance estimation. The available `gaussian_root_region.py` diagnostic replays average linkage on **squared Euclidean** distances. It does not replay the actual adaptive diffusion event.

For **unknown common scalar variance** under the stronger root null that every observation has the same mean, conditioning additionally on centered total norm and residual direction eliminates sigma. The resulting ratio has a truncated F(p,(n-2)p) law. Our implemented scalar special case uses h=|a|/||w|| and (n-2)h²~F(1,n-2), truncated to its replayed merge-sequence region. It is not ordinary untruncated F applied after clustering. This construction and the stronger null are described in [Yun and Barber, Section 3, Theorem 2 and K=2 discussion](https://arxiv.org/html/2301.12999v3#S3). Merely equal cluster-average means allows within-cluster nuisance means and noncentral residual scatter, so that stronger-null derivation does not follow.

Unknown unrestricted feature covariance requires another nuisance-elimination argument or a justified independent covariance experiment. Conditioning on a plug-in estimate S does not make S the true covariance: conditional on S, the independent contrast still has covariance s Sigma. A Wishart estimate can yield an F pivot only under the required joint Gaussian/Wishart independence and projection conditions. If projection, membership or variance estimate is selected using the same observations, those conditions must be re-established or their selection included in the conditioning event.

### What the production fitted scale actually does

`empirical_null_inflation_estimation.py` fits the weighted ratio

`c_raw = sum_i w_i Q_i / sum_i w_i (reference_scale_i * df_i)`.

This resembles a common-scale gamma likelihood estimator only under a suitable fixed-weight reference model. Here the records, roles and weights are selected from the same hierarchy. The source restricts fit support to positive-df, positive-weight records with a null-like/stopped-edge role. That role is assigned from edge decisions; it is not proof that their scientific nulls are true. `child_parent_edge_metadata.py` uses the product of adjusted child-edge p-values as a weight, or weight 1 for untested/stopped edges. Descendants may share the ancestor stopping event. Dependency groups reduce support-count duplication but do not prove independence or the fitted tail law. The estimator can contain a focal record if it qualifies for calibration; no universal leave-one-focal exclusion is implemented. Focal `is_null_like` records are skipped by the active adjusted-test loop.

Prediction matches feature family and smooths over log projection dimension and log parent count; Bernoulli disables the parent-count axis. Bandwidth is the weighted context standard deviation. Both family baseline and local prediction have lower floor 1. If a target parent lies above every calibration parent and its local factor exceeds the family baseline, the factor is capped at that baseline. The evaluated rule remains `chi2.sf(Q/(reference_scale*c_hat),df)`. Threshold diagnostics for group counts, effective weights and leave-one-group sensitivity establish internal support policy, not a sampling distribution. Enforcing them is optional in the decision function. The current docstrings explicitly mark the selected-tail law unvalidated.

Matching a mean scale does not match covariance anisotropy, tail shape or selected truncation. Since c_hat>=1, the adjusted chi-square p-value is pointwise no smaller than the raw chi-square p-value, but the raw p-value is already invalid in the selected examples. An independent fitted scale would still have estimation uncertainty. In the special auxiliary model Q/c~chi-square(k), and independently V/c~chi-square(m), Q/(V/m) has k F(k,m), not chi-square(k). Arbitrary kernel weights and selected overlapping records do not supply that auxiliary model.

There is nevertheless a precise **conditional-law construction** once a nuisance-free path has been justified. Let X(r) be that path, f0 its known chi or angular/F-induced density, and R the set reproducing the declared selection. Recompute every fitted object along X(r), including covariance regime, projection/rank, edge roles, calibration pool, weights, bandwidth and cap. For T(r)=Q(X(r))/c_hat(X(r)), its conditional upper-tail p-value is

`p = integral_R 1{T(r) >= T(r_observed)} f0(r) dr / integral_R f0(r) dr`.

This is a pushforward of the valid conditional path law, not chi-square(T_observed). The numerator may be a disconnected level set because refitting c_hat and rank need not preserve monotonicity in r. Holding observed c_hat fixed instead computes another statistic and does not establish the fitted statistic's conditional law. Undefined support cases and any event selecting which focal test is reported must be included in the declared conditioning/reporting contract. The formula is a mathematically specified route for known covariance or the restricted unknown-scalar-variance model; an efficiently resolvable actual diffusion event and its fitted-statistic level sets have not been implemented or validated here.

**Minimal next construction for Q5:** first use the established known-I, full-space root conditional law with the *actual* diffusion selection event (owned by the selection-event work). Then add unknown scalar sigma through a declared common-mean root null and the norm/direction conditioning above, replaying that same event. Keep empirical c_hat outside the inferential construction until a joint selected model or independent predictive calibration law is established. A separate independent-discovery/independent-test experiment can freeze tree and basis, project raw test rows into that fixed subspace and use their ordinary pooled Hotelling F when residual df and covariance rank support it. None of these proposals is already a production guarantee.

## Q6. Does validity survive projection, multiple testing and traversal?

**Answer:** a valid root selective p-value establishes one conditional root claim. It does not by itself establish calibrated projected node p-values, hierarchical FDR or a bound on erroneous final fragmentation. Those require declared hypotheses, jointly applicable p-values and an error metric matched to the actual reporting/traversal rule.

### Projection and rank

For a genuinely isotropic null contrast and a fixed independent orthonormal k-dimensional subspace, projected energy has chi-square(k) law. An independent raw projection with unknown covariance can instead use an appropriate projected test-row covariance and finite-sample F construction. Simply calling a same-data PCA matrix fixed does not create independence. The production sibling context uses all leaves; edge splitting after selected membership does not remove membership selection.

`projection_dimension_estimators.py` compares eigenvalues to the MP boundary, applies a minimum dimension and limits positive available directions. `projected_wald_kernel.py` can choose the smallest PCA prefix containing a requested fraction of the **observed contrast's** energy. Both choices are random when learned from tested observations. Extending ranks through complete numerical eigenvalue groups ensures invariance to rotations inside repeated eigenspaces; it does not establish a chi-square reference. A selected radial construction must freeze the basis/rank using appropriate independent or nuisance information, or account for their selection along the path. MP rank is spectral metadata, not an exact finite-sample p-value theorem for these selected, estimated-covariance rows.

### Multiplicity implemented, and the missing theorem match

`multiple_testing/tree_bh/correction.py` applies BH to each reached child-edge family and multiplies subsequent family alpha by the ancestor family's rejection proportion. Only rejected edges lead to descendant families. Its corrected p-values are family-local BH adjustments compared against the family's own alpha.

The primary [TreeBH paper, Theorems 1–3](https://arxiv.org/pdf/1705.07529) controls explicitly defined selective levelwise errors under valid p-values and specified dependence/selection conditions. Theorem 1 gives PRDS conditions involving within-family and certain nonancestor p-values; arbitrary ancestor dependence alone is not disqualifying. Its intersection/Simes specialization has additional structure. No corresponding joint PRDS/selection proof was established here for diffusion-selected tree, learned PCA, estimated covariance and shared fitted inflation. The general-dependence variant in the paper is a different procedure, not an automatic property of unmodified BH. Ordinary BH itself needs valid p-values and independence or appropriate positive dependence; see [Benjamini and Yekutieli, Theorems 1.2–1.3](https://www.math.tau.ac.il/~ybenja/MyPapers/benjamini_yekutieli_ANNSTAT2001.pdf).

Scientific interpretation also matters: child-parent **mean equality** is not a nested homogeneity hypothesis. A child containing equal-sized subgroups with means +a and -a can have the same mean 0 as its parent, while its grandchildren differ. Thus the intersection-null interpretation cannot be assumed. This does not by itself rule out every general selective-family theorem; it means the hypotheses and theorem mapping must be specified rather than inferred from tree topology.

`inflated_projected_wald_annotation/fdr_annotation.py` builds sibling families from the **edge-reachable frontier** and applies BH at the same sibling alpha at each depth. The next frontier includes children of every eligible edge-open candidate parent, irrespective of its sibling rejection. It does not spend alpha across depths or implement the edge TreeBH rejection-proportion recursion. Frontier eligibility and adjusted sibling p-values share edge/tree information. Applying BH separately to selected depth families is not itself a global FDR or partition-error theorem. [Benjamini and Bogomolov](https://arxiv.org/abs/1106.3670) specifically distinguish within-family control from error among selected families.

### Traversal and the required error target

`gate_evaluator.py` requires binary structure and at least one significant child edge. Without passthrough, failure of the local sibling gate stops the subtree. With passthrough, a locally closed sibling gate can be crossed when a reachable descendant has an open split and optional support guards permit it. `SPLIT` and `PASS_THROUGH` both recurse. These rules change which final boundaries are reported; support guards are policy constraints, not established selective probability bounds.

At least four distinct targets must be separated:

1. **Root selective type-I error:** conditional on the declared selection and nuisance quantities, Pr(p_root<=alpha | selection,nuisance,H_root)<=alpha.
2. **Edge/node FDR:** E[V/max(R,1)] for a specified collection of rejected mean-contrast or homogeneity hypotheses; the collection and nulls must be stated.
3. **Selective levelwise/family FDR:** the weighted selected-family criterion of a particular hierarchical theorem; it is not automatically ordinary all-edge FDR.
4. **Global-null fragmentation:** under iid common-mean rows, Pr(number of final clusters>1). It concerns a whole dataset, not one node or two correlated root edges. False final separation under partial alternatives needs another definition, such as splitting a truly homogeneous subtree.

Without passthrough, if any final fragmentation necessarily requires rejecting one valid root gate, root type-I control bounds global-null fragmentation by that gate's alpha, even though it proves nothing about descendant FDR. Passthrough defeats a bound based only on the local root sibling p-value: descendant evidence may cross a nonrejected root sibling gate. The measured difference between 28 root sibling rejections and 29 fragmented datasets in the prior native factorial illustrates why root and partition outcomes must be counted separately.

**Minimal next validated construction for Q6:** declare global-null fragmentation as the first target, use the actual-selection root p-value, and make root nonrejection stop all traversal. This yields a direct implication proof, with no claim of local/partial-alternative guarantees. If traversal including passthrough is required, a separate benchmark construction can freeze the entire tree and raw projection using independent discovery data, derive valid test-data p-values for every predeclared internal node, and apply one all-node Bonferroni/Holm family. A union bound then controls any false rejection without dependence assumptions, despite later choosing a subset of rejections. To interpret every passthrough boundary as a homogeneity claim, use nested homogeneity nulls and prove that a split within a homogeneous region requires a rejected true null; control of mean-contrast rejections alone does not establish this stronger partition claim. This construction is a conservative validation target, not an authorized change to production thresholds or procedures.

## Unresolved boundary and evidence status

The actual diffusion selected event, nuisance-free law for all covariance regimes with shared fitted c_hat, and error control for the present combined edge/sibling/traversal procedure remain unresolved. The independent Gaussian laws and restricted analytic root constructions are established under their stated assumptions; they are not universal clustering validity results. No production files, thresholds, dependencies or services were changed for these answers.
