# What calibration still needs

The immediate necessary step is to measure false splits for the **complete
restored diffusion NNLS method on independent simulated datasets with a known
null**, alongside power on signal datasets. Restoring the empirical rule is a
clustering behavior decision. Its mathematical status has not changed.

## The actual requirement

The old rule estimates a nondeflating local scale and returns

\[
\hat c_u=\max\left(1,
  \frac{\sum_q w_q(u)T_q}{\sum_q w_q(u)a_q\nu_q}\right),
\qquad p_u=\Pr\{\chi^2_{\nu_u}\ge T_u/(a_u\hat c_u)\}.
\]

For valid null p-values, the relevant condition is
`Pr(p_u <= t | the declared selection/null context) <= t`. Exact uniformity
is stronger than necessary: conservative, nonuniform p-values can satisfy this
condition. A rejected uniformity KS test alone is therefore **not a reason to
disable the method**. Conversely, a high ARI and stable mean inflation do not
show that its small p-value tail satisfies this inequality. Traversal-aligned
BH also needs its own error-control argument under the actual selection and
dependence; a per-node tail check alone cannot establish final-method FDR.

Selected clustering changes the hypothesis/test distribution. Gaussian
hierarchical-clustering selective inference illustrates this issue, but does
not supply a ready-made theorem for this method's diffusion and NNLS stages.
[Gao, Bien and Witten](https://arxiv.org/abs/2012.02936).

## Concrete next experiment

1. Freeze code, input generation, alpha, geometry, solver parameters, projection
   rules, calibration weights, bandwidths and traversal. Use the restored
   pydiffmap average-linkage NNLS preset first; treat graphtools and other
   linkages as separate methods. Define the null and the intended error rate
   before generating results. Keep unsupported geometry explicit.
2. Generate independent datasets from declared Bernoulli, categorical and
   supported continuous null models, including feature dependence and sample
   imbalance relevant to the intended use. Native continuous inputs are
   incompatible with the Hamming preset; discretized Gaussian cases test a
   different representation. Use planted signal and mixed null/signal trees to
   examine power and false descendant splits. Assign null truth from the
   generator, not from whether the method failed to reject an edge.
3. On **every replicate**, refit diffusion, construct the tree, refit NNLS branch
   lengths, estimate node covariance/PCA/rank, rerun edge selection, select
   internal calibration records and weights, refit inflation, then run sibling
   BH and traversal. Freezing the observed tree tests a different experiment.
   Do not replace the empirical rule by fixed-coordinate BH for this check.
4. Record raw and adjusted statistics/p-values, selected focal contexts, truth,
   weights, dependency groups, internal support, reached nodes, split decisions,
   final K and ARI. Report tail exceedance at alpha 0.005, 0.01 and 0.05; false
   root/descendant splits; `Pr(K > 1)` under a whole-dataset one-cluster null;
   and replicate-level false-discovery proportion on mixed null/signal data.
   These are different estimands and must not be substituted for each other.
5. Keep dataset/seed replicates independent. Nodes in one hierarchy overlap;
   hundreds of selected nodes from one dataset are not hundreds of independent
   simulations. Report context eligibility, support and unsupported rates
   separately, and provide both conditional-on-eligible and unconditional
   outcomes. Closing every gate can suppress false splits while destroying
   power, so include paired signal recovery and cluster-count errors.
6. Split **datasets/seeds** into development and held-out evaluation pools.
   Choose bandwidth/threshold/tail corrections on development data only, then
   evaluate the frozen choice on held-out seeds with uncertainty estimates
   based on independent datasets. Diagnose selected-null contamination by
   checking generator truth among weighted calibration records. Test sensitivity
   to alpha, group concentration, sample size, dimension and feature family.
7. Size the study for the intended precision and number of contexts. The local
   external-tail diagnostic requires 499 *matching independent simulations* at
   alpha 0.01 for its resolution rule `1/(m+1) <= 0.2*alpha`; that is a resolution
   floor, not a power/precision guarantee. For an independent Bernoulli error
   indicator near probability 0.01, approximately 9,900 replicates give standard
   error 0.001 (`B = p(1-p)/SE^2`). Rare eligible contexts require more attempted
   datasets. Prespecify tolerance, confidence bounds and multiplicity across
   contexts; do not stop simply when a favorable estimate appears.

If error exceeds the declared target, the next repair is a selected-tail
correction validated on held-out datasets, or a narrower formally derived
selection law. Do not promote an external correction from a mean fit alone.
For an exchangeable Monte Carlo reference, use a valid finite-reference rank
calculation such as `(1 + exceedances)/(B + 1)` where its assumptions apply;
estimated-null or mismatched selected contexts do not become exact through
that formula. [Phipson and Smyth](https://gksmyth.github.io/pubs/PermPValuesPreprint.pdf).

## What existing evidence establishes

The existing 500-replicate external selected-hierarchy study rebuilds hierarchy
and edge selection, but is a different diagnostic from the restored diffusion
NNLS pipeline. Its raw CSV has no admissible external-calibration stratum and
at most 470 matching independent simulations per stratum (499 required). Among
strata with at least 100 matching simulations, reported scalar-tail rejection
at alpha 0.01 is zero. This supports conservatism in that diagnostic, not proof
that the restored rule is valid or invalid. Its small uniformity KS p-values
show distributional mismatch, which must be distinguished from anti-conservative
tail error. Recomputed details: `external_contract_check.json`.

The earlier matched support panel demonstrates that fixed-coordinate gates can
recover signal with zero internal support, but also split null datasets. It is
evidence against treating “unsupported” as a clustering-quality score; it is
not validation of those replacement gates. The restored empirical rule still
requires some supported calibration records. The no-support benchmark outcome
is retained, and no neutral scale or replacement gate is inserted.

The exact independent common-scale F API remains valid only under its explicit
mutual-independence, chi-square, common-scale, fixed-unit-weight and observation
ownership assumptions. Independent external seeds alone do not make the focal
data-selected statistic chi-square. The empirical entry point now rejects exact
F models so it cannot silently bypass those assumptions.

## Evidence and scope

- `tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflation_correction/empirical_null_inflation_estimation.py`
- `tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflated_projected_wald_annotation/fdr_annotation.py`
- `raw/assets/benchmark-results/selected_hierarchy_external_contract_20260602_500/{manifest.json,external_calibration_contract.csv}`
- `reports/diffusion_nnls_versions_20260909/{README.md,support_review_results.csv}`
- `raw/inbox/calibration-restoration-reference-notes-20260909.md`

This document specifies the necessary next validation study. It does not claim
that the study has run or that a selected-tail/FDR guarantee has been proved.
