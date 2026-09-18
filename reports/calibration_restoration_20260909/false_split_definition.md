# False splits and required data conditions

CAL-01 definition, 2026-09-10. Simulation truth and coverage are specified here;
calibration itself remains to be validated.

## Definition

A **local false split** is a sibling rejection used by traversal when the
generator-defined sibling null is true:

\[
H_0^{\mathrm{sib}}(u):\theta_{L_u}=\theta_{R_u}.
\]

`theta` means the population feature parameters, not observed child averages:
Bernoulli probabilities, categorical probability vectors, or continuous means
in the declared feature representation. This is the sibling null in
`manuscript/sections/method/sibling_test.tex`. Equal continuous means need not
mean equal distributions; equal binary/categorical marginal probabilities need
not mean equal joint distributions when features are dependent.

For a fixed feature representation and known per-observation generator
parameters `theta_i`, the population parameter of a selected child A is

\[
\theta_A=\frac{1}{|A|}\sum_{i\in A}\theta_i.
\]

The selected sets depend on observed data; `theta_i` remain generator parameters.
Do not replace them with sample averages, majority labels or observed p-values.
This population-versus-empirical distinction follows equations (2)–(3) of
[Gao, Bien and Witten](https://arxiv.org/html/2012.02936v3#S1).
Their Gaussian result does not establish calibration for diffusion NNLS.

### Three outcomes to record separately

| Outcome | Definition |
|---|---|
| Local false sibling rejection | A true sibling null is rejected; record whether the test was reached and whether traversal used its rejection. |
| Final false partition under a global null | One declared homogeneous population produces a returned partition with `K_hat > 1`. Count once per dataset when estimating its probability. |
| Fragmentation on signal data | One true population appears in multiple returned clusters. This is a recovery error, not automatically a false sibling rejection. |

The implementation distinguishes `pass_through` from `sibling_gate_open_split`.
Passing through a node to reach descendants is not itself a sibling rejection;
an unvisited annotated rejection is not a split used by traversal.

Under a global null, every observation has the same population parameter.
Sampling noise can produce visibly separated children without making their
population null false. A successfully returned `K_hat = 1` has no final
false-partition event. Unsupported/skipped/error runs are unavailable outcomes,
not successful one-cluster results; their rates must be reported separately.
CAL-02 still determines conditional and unconditional error estimands.

### Truth in mixed null/signal datasets

A parent containing only one true homogeneous population is a clear null
context. Two pure children from populations with different parameters form an
alternative. Mixed children require their generator-defined parameter averages:

- Both children come from a mean-zero population: true null, regardless of
  their observed means.
- Both children contain equal proportions of populations with means -1 and +1:
  both population means are zero, so the local mean null is true.
- Different proportions of those two populations give unequal means, even
  when both children have the same majority label.

Thus `K_hat > K_true` alone cannot identify false sibling splits: fragmentation
and merging can coexist while the total K stays correct. Use selected-node
truth and the final label contingency table, alongside ARI and K errors.

If a generator/representation lacks an exact population-parameter truth
calculation for mixed children, mark those contexts unclassified for local
Type I error. They may still contribute to partition-recovery metrics.

Sample-median binarization and other data-fitted encodings add selection.
Refit them in each replicate and retain provenance. Start the exact local-truth
panel with native binary/categorical generators or fixed population cut points.
Study sample-fitted encodings separately; raw continuous means do not establish
their mixed-child encoded truth. Pure-population nulls remain available when
the generator and encoding preserve exchangeability within that population.

## Data conditions that must be represented

Establish results separately for each supported feature representation and
method preset. Distinguish baseline model conditions from explicit model-stress
conditions; success under one does not establish validity under the other.

| Condition | Required coverage |
|---|---|
| Population structure | A global one-population null; separated and overlapping signal populations; homogeneous descendants inside signal datasets. |
| Feature representation | Native Bernoulli; one-hot categorical with multinomial blocks; supported continuous inputs; separately labeled discretized Gaussian variants. Native continuous data are incompatible with the Hamming preset. |
| Marginal probabilities | Rare and frequent binary events; near-boundary probabilities; balanced and rare categorical levels; low/high category counts; constant coordinates. |
| Sample size and balance | Small and large datasets; tiny selected children; balanced and strongly imbalanced populations and child sizes. |
| Dimension and signal | Feature dimension below, near and above sample size; sparse/dense differences; weak/strong separation; many irrelevant features. |
| Feature dependence and covariance | Independent features as a baseline; correlated blocks; low-rank/anisotropic covariance; unequal covariance with equal means where applicable. Covariance-only differences remain a mean null, not automatically mean signal. |
| Observation dependence | Independent observations first; separately labeled repeated, batch-correlated or phylogenetic observations when such use is claimed. Independence across simulation datasets does not make observations inside each one independent. |
| Contamination and geometry | Known outliers/nuisance contamination; duplicates/ties; nearly singular spectra and degenerate distances. Declare whether outliers are truth populations or contamination before scoring. |
| Calibration support | Zero, weak and substantial support; concentrated weights; records sharing dependency groups. Measure support as an outcome of full replay rather than manufacturing it by deleting calibration records. |

Include interactions such as rare categories with small children, correlated
high-dimensional noise, and strong signal with little internal null support.
Report unsupported geometry and execution failures as coverage limitations.
Exact grid values, replicate counts, precision and acceptance limits remain
CAL-08, CAL-09, CAL-11 and CAL-12. Each simulated cell must retain its population
truth, feature/observation model, representation, seed, parameters and status.

## Evidence and implementation boundary

- `manuscript/sections/method/sibling_test.tex` specifies the sibling null.
- `tree_break_selection/tree/feature_space.py` specifies the three feature-block
  families and their covariance/contrast representations.
- `_true_context_labels` in
  `benchmarks/diagnostics/calibration/statistics/mixed_internal_calibration_sweep.py`
  identifies pure-parent nulls but uses differing majority labels as a signal
  diagnostic. That signal flag is not exact mixed-child population truth.
- `benchmarks/validation/statistics/selected_edge_type1_geometry.py` records
  the global-null final event as `found_clusters > 1`, with run status.
- `tree_break_selection/hierarchy_analysis/decomposition/gates/gate_evaluator.py`
  distinguishes sibling-open splits from pass-through traversal.
- `benchmarks/shared/generators/gaussian_cases.py` separates continuous and
  sample-median-binarized Gaussian representations.
- `reports/diffusion_nnls_versions_20260909/README.md` records support, geometry
  and alpha limitations motivating the coverage above.

This is a definition and coverage specification. Exact mixed-context truth
calculations and the full held-out simulation panel still need implementation
and evaluation; no runtime behavior is changed here.
