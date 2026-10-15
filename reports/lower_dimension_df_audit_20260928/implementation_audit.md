# Projected-Wald dimension and degrees-of-freedom implementation audit

## Finding

For the audited production paths, a smaller retained projection dimension is
carried through consistently as the projected-Wald degrees of freedom. I found
no path that substitutes the ambient feature count or raw Marchenko-Pastur
signal count for the dimension actually used by the quadratic statistic.
However, this bookkeeping does not justify the chi-square reference when the
projection rank is selected from the same tested contrast. The remaining issue
is calibration of the reference law, not an off-by-one or ambient-dimension df
assignment.

## Dimension-to-statistic trace

The feature space defines the actual contrast coordinate count: categorical
blocks contribute `raw_dimension - 1`, while Bernoulli and continuous blocks
contribute their declared dimensions ([feature_space.py:43](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/tree/feature_space.py:43), [feature_space.py:134](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/tree/feature_space.py:134)). Contrast construction whitens those coordinates with the specified covariance blocks before returning `z` ([contrast_covariance.py:49](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/statistics/contrast_covariance.py:49), [contrast_covariance.py:314](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/statistics/contrast_covariance.py:314)). Thus the spectral basis and test contrast share the same tangent-coordinate width.

Node-local spectral decomposition removes constant coordinates. In the dual
case it caps available directions at `n_samples - 1`; the MP dimension estimate
also caps its test dimension to positive available eigen-directions
([decomposition.py:25](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/decomposition/backends/eigen/decomposition.py:25), [projection.py:56](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/decomposition/backends/eigen/projection.py:56), [projection_dimension_estimators.py:161](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/statistics/projection/projection_dimension_estimation/projection_dimension_estimators.py:161)). PCA projection rows are recovered in the full contrast width and are orthonormal eigendirections ([projection.py:12](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/decomposition/backends/eigen/projection.py:12)).

The kernel applies the selected basis to `z`, sums the squared retained
coordinates, and passes that projected vector to the common reference helper
([projected_wald_kernel.py:135](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_kernel.py:135), [projected_wald_kernel.py:154](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_kernel.py:154)). The helper sets `df = len(projected_vector)` and computes `chi2.sf(T, df)` with scale 1 ([projected_wald_reference_distribution.py:84](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_reference_distribution.py:84)). The rank is therefore the statistic's actual retained rank, including any tie-group extension.

## Caller and edge-case audit

- **Edge tests:** child-parent contrasts are whitened, sent to the shared kernel, and return its df directly ([child_parent_projected_wald_test.py:35](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/statistics/child_parent_divergence/child_parent_projected_wald/child_parent_projected_wald_test.py:35), [child_parent_projected_wald_test.py:53](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/statistics/child_parent_divergence/child_parent_projected_wald/child_parent_projected_wald_test.py:53)). Tree annotation stores that returned value ([tree_testing.py:94](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/statistics/child_parent_divergence/child_parent_divergence_annotation/tree_testing.py:94)).
- **Sibling tests:** sibling contrasts use the same covariance whitening and shared kernel (`sibling_z_scores.py:27`, removed 2026-10-15 and inlined into `sibling_divergence_test.py`; [sibling_divergence_test.py:58](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/pair_testing/wald_statistic/sibling_divergence_test.py:58)). Pair records store the returned df as both `degrees_of_freedom` and `sibling_projection_dimension` ([record_collection.py:182](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/pair_testing/collection/record_collection.py:182), [record_collection.py:259](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/pair_testing/collection/record_collection.py:259)).
- **Inflation:** empirical sibling fitting uses `reference_scale * df` for the expected statistic; plug-in adjusted p-values preserve the record df; the independent exact-F route checks and uses the focal df ([empirical_null_inflation_estimation.py:230](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflation_correction/empirical_null_inflation_estimation.py:230), [empirical_null_inflation_estimation.py:279](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflation_correction/empirical_null_inflation_estimation.py:279), [empirical_null_inflation_estimation.py:567](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflation_correction/empirical_null_inflation_estimation.py:567), [empirical_null_inflation_estimation.py:621](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflation_correction/empirical_null_inflation_estimation.py:621)). This keeps scale/rank arithmetic aligned but inherits the validity limits of the input reference law.
- **Zero rank:** exact zero-dimensional kernel context is accepted only for a zero contrast, returning statistic 0, df 0, scale 1, and p-value 1 ([projected_wald_kernel.py:71](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_kernel.py:71)). An all-constant spectral node returns dimension zero. Adaptive selection with a nonzero candidate basis instead floors selected rank at one ([projected_wald_kernel.py:122](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_kernel.py:122)).
- **Ties:** rank targets extend through the complete numerical eigenvalue group in both MP and Wald-basis selection ([eigenvalue_groups.py:20](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/eigenvalue_groups.py:20), [projected_wald_projection_basis.py:30](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_projection_basis.py:30)). This can make actual df exceed the nominal target, but the expanded basis is also what forms the statistic, so the reported df still matches it. This addresses rotation invariance within tied eigenspaces; it does not validate the adaptive-rank reference law.
- **Adaptive energy option:** the setting defaults to `None` in gate metadata and orchestration ([annotation_bundle.py:36](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/decomposition/gates/annotation_bundle.py:36), [orchestrator.py:138](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/decomposition/gates/orchestrator.py:138)). When enabled, the kernel selects the prefix reaching the requested fraction of the tested contrast's projected energy, then uses that selected rank for both the statistic and df ([projected_wald_kernel.py:102](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_kernel.py:102)). The option is explicitly set to 0.90 in the pancreas benchmark configuration ([pancreas_benchmark.py:239](/Users/berksakalli/Projects/kl-te-cluster/applications/scrna/pancreas_benchmark.py:239)).

The reference-distribution module states the governing limitation: chi-square
rank `k` follows for an isotropic standardized null conditional on a fixed
projection not selected from the tested contrast; it does not account for
selection of PCA rows and projection dimension from the same data
([projected_wald_reference_distribution.py:10](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_reference_distribution.py:10), [projected_wald_reference_distribution.py:15](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_reference_distribution.py:15)). The unit test for adaptive energy checks selected rank, df, and statistic, while the reference test checks the fixed-vector formula; those are bookkeeping/contract checks rather than calibration evidence ([39_test_projected_wald_statistics.py:63](/Users/berksakalli/Projects/kl-te-cluster/tests/statistics/39_test_projected_wald_statistics.py:63), [49_test_projected_wald_reference_distribution.py:29](/Users/berksakalli/Projects/kl-te-cluster/tests/statistics/49_test_projected_wald_reference_distribution.py:29)).

## Verification and scope

Inspected `pyproject.toml`, `Makefile`, `.venv/pyvenv.cfg`, and the pytest
configuration before running tests. The repository uses Python 3.11.12 in
`.venv`, `uv` 0.9.5, pytest 8.4.2, and the configured `tests/` suite. The
read-only baseline command

```text
.venv/bin/python -m pytest -q tests/statistics/39_test_projected_wald_statistics.py tests/statistics/49_test_projected_wald_reference_distribution.py tests/statistics/42_test_contrast_covariance.py tests/validation/calibration/194_test_reference_law_reporting.py
```

passed all 41 tests. This bounded check verifies implementation contracts; it
does not test the adaptive selected-rank null law. The worktree already
contained extensive unrelated and in-progress changes, which were left intact.
