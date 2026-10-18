# Passthrough p-value normalization audit

## Scope and result

This audit traces the current endotype adaptive-diffusion subspace pipeline
from input matrix through tree geometry and TBS statistics. The tree geometry
and statistical feature coordinates are separate inputs: diffusion and cosine
normalizations determine the tree, while TBS receives the original `data`
matrix for node distributions and tests. Therefore the geometry can change
which node contrasts are tested, but it does not directly rescale a fixed
node's raw distributions when `passthrough` changes.

The inspected working tree contains pre-existing edits in eigenvalue grouping,
projection and statistics files. Findings below describe those live files; this
audit changed only this report. Environment inspection confirmed macOS arm64,
the project `.venv` Python 3.11.12, NumPy 2.3.4, SciPy 1.16.3, pandas 2.3.3,
pydiffmap 0.2.0.1 and pytest 8.4.2. Baseline `make wiki-lint` passed for 203
pages; focused pytest collection discovered 23 tests and did not run them.

## Data to tree geometry

In the current GO workflow, `weight_feature_matrix` returns raw values for
`binary` weighting, or applies TF-IDF with `norm=None`; the cosine operator
then divides each row by its Euclidean norm and forms the row cosine Gram
matrix. It sorts its symmetric eigensystem in descending order and retains
positive modes, extending the target rank through a numerical tie
([adaptive_cosine.py](../../tree_break_selection/space_separation/adaptive_cosine.py#L49),
[adaptive_cosine.py](../../tree_break_selection/space_separation/adaptive_cosine.py#L74)).
For each spectral block, the workflow extracts those sample coordinates,
constructs adaptive diffusion distances, and clusters those distances by
average linkage ([adaptive_diffusion_go_ic_workflow.py](../../applications/endotypes/pipelines/adaptive_diffusion_go_ic_workflow.py#L137),
[adaptive_diffusion_go_ic_workflow.py](../../applications/endotypes/pipelines/adaptive_diffusion_go_ic_workflow.py#L234)).

The diffusion path uses `pydiffmap.kernel.Kernel` for variable bandwidth and
Gaussian affinity construction, but it does not call
`pydiffmap.DiffusionMap`. The returned neighbor affinity is converted to a
dense matrix; `compute_diffusion_coordinates` zeros its diagonal, symmetrizes
with the elementwise maximum of `W` and `W.T`, forms
`S = D^(-1/2) W D^(-1/2)`, and computes coordinates `U lambda^t` from the
symmetric eigensystem ([diffusion.py](../../tree_break_selection/space_separation/diffusion.py#L162),
[diffusion.py](../../tree_break_selection/space_separation/diffusion.py#L278)).
For positive degrees, if `P = D^(-1) W`, an eigenpair `S u = lambda u`
corresponds to right eigenvector `D^(-1/2)u` and left eigenvector
`u.T D^(1/2)` of `P`; these are biorthogonally paired. The current coordinates
use `u`, not the right eigenvector of `P`, so their unweighted Euclidean
distances define the symmetric-coordinate geometry.

The installed pydiffmap 0.2.0.1 source shows that its full `DiffusionMap`
constructs a density estimate `q`, applies its `q^(-alpha)` normalization,
row-normalizes to a Markov matrix, and builds a generator. The production
helper instead uses the kernel affinity followed by symmetric degree
normalization, with no `q^(-alpha)` step ([installed
diffusion_map.py](../../.venv/lib/python3.12/site-packages/pydiffmap/diffusion_map.py#L117),
[installed diffusion_map.py](../../.venv/lib/python3.12/site-packages/pydiffmap/diffusion_map.py#L135),
[installed diffusion_map.py](../../.venv/lib/python3.12/site-packages/pydiffmap/diffusion_map.py#L142)).
This is a demonstrable implementation distinction and a valid alternative
geometry; the inspected code does not establish that it is a defect. Its
effect is on tree topology and distances, not a direct change to a fixed
node's feature normalization.

The production handoff confirms the distinction: `run_current_kl(data,
distances, ...)` passes original `data` and the diffusion-derived distances to
`run_tbs_on_distance`; the runner builds the tree from distances, then calls
`tree.populate_node_divergences(data_df, ...)` and supplies `leaf_data=data_df`
to gate annotation ([adaptive_diffusion_go_ic_workflow.py](../../applications/endotypes/pipelines/adaptive_diffusion_go_ic_workflow.py#L73),
[adaptive_diffusion_go_ic_workflow.py](../../applications/endotypes/pipelines/adaptive_diffusion_go_ic_workflow.py#L263),
[tbs_runner.py](../../benchmarks/shared/runners/tbs_runner.py#L125),
[tbs_runner.py](../../benchmarks/shared/runners/tbs_runner.py#L207)).

## Raw feature distributions and covariance whitening

`FeatureSpace` makes the statistical coordinate contract explicit:
Bernoulli blocks use their single raw coordinate, categorical blocks use a
drop-last simplex chart with `K-1` contrast coordinates, and continuous blocks
use empirical Gaussian covariance ([feature_space.py](../../tree_break_selection/tree/feature_space.py#L13),
[feature_space.py](../../tree_break_selection/tree/feature_space.py#L41)).
Tree leaves store raw rows and internal nodes store leaf-count-weighted
subtree means ([distributions.py](../../tree_break_selection/tree/distributions.py#L119),
[distributions.py](../../tree_break_selection/tree/distributions.py#L160)).
Categorical one-hot column names require an explicit feature-space declaration;
the inference helper rejects them without one ([feature_space.py](../../tree_break_selection/tree/feature_space.py#L283)).

Sibling Wald contrasts are left-minus-right. The pooled null distribution is
the leaf-count-weighted sibling mean, with variance multiplier
`1/n_left + 1/n_right`. Bernoulli variance is `p(1-p)`; the categorical
drop-last covariance is `diag(p) - p p.T`; both are multiplied by that
variance factor and receive the configured small ridge. Continuous features
use the parent's stored empirical covariance multiplied by the same factor
([contrast_covariance.py](../../tree_break_selection/hierarchy_analysis/statistics/contrast_covariance.py#L167),
[contrast_covariance.py](../../tree_break_selection/hierarchy_analysis/statistics/contrast_covariance.py#L757),
[contrast_covariance.py](../../tree_break_selection/hierarchy_analysis/statistics/contrast_covariance.py#L855)).
Whitening solves with the Cholesky factor of these explicit covariance blocks;
it is not columnwise z-scoring of all input columns ([contrast_covariance.py](../../tree_break_selection/hierarchy_analysis/statistics/contrast_covariance.py#L49)).

For continuous data, node scatter is divided by `max(N-1, 1)` when stored.
The default `parent_total` policy uses that parent's empirical covariance,
which includes between-child mean variation. The opt-in
`guarded_within_child` policy uses the degrees-of-freedom-weighted immediate
child covariance only when both child leaf counts pass the guard, otherwise
falls back to parent-total covariance ([distributions.py](../../tree_break_selection/tree/distributions.py#L265),
[distributions.py](../../tree_break_selection/tree/distributions.py#L371)).
These are declared covariance policies. The prior known-branch-skew analysis
also records that parent-total covariance includes between-child mean
variation ([known-branch-skew-diagnostic.md](../../wiki/sources/known-branch-skew-diagnostic.md#L49)).

The node-local spectral basis has a separate normalization path. Raw rows are
first centered at the null distribution and whitened by the feature-family
null covariance. Constant features are removed; the active whitened data are
centered, and the primal covariance or dual Gram matrix divides by `N`, not
`N-1`. The dual branch recovers feature-space eigenvectors from sample-space
eigenvectors. Both are symmetric eigensolutions of the same feature covariance
operator ([contrast_covariance.py](../../tree_break_selection/hierarchy_analysis/statistics/contrast_covariance.py#L236),
[preparation.py](../../tree_break_selection/hierarchy_analysis/decomposition/backends/eigen/preparation.py#L21),
[operators.py](../../tree_break_selection/hierarchy_analysis/decomposition/backends/eigen/operators.py#L8),
[decomposition.py](../../tree_break_selection/hierarchy_analysis/decomposition/backends/eigen/decomposition.py#L18),
[projection.py](../../tree_break_selection/hierarchy_analysis/decomposition/backends/eigen/projection.py#L91)).
The `N` denominator is the covariance spectrum used by the MP threshold; the
`N-1` denominator is the stored continuous sample covariance used in Wald
contrasts. These different estimands are implemented intentionally and are
not, by themselves, a normalization defect.

The test statistic projects the covariance-whitened contrast onto the
orthonormal feature-space PCA rows and uses its squared norm with a
chi-square reference whose degrees of freedom equal the actual projected rank
([projected_wald_kernel.py](../../tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_kernel.py#L135),
[projected_wald_reference_distribution.py](../../tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_reference_distribution.py#L6)).
The p-value documentation limits that reference to a fixed projection
independent of the tested contrast. The current node PCA and, in some modes,
rank selection use the same data as the test. That is a separate calibration
limitation; changing traversal passthrough does not repair it.

## Calibration-sensitive details

The statistical contract standardizes Bernoulli features by their marginal
`p(1-p)` variances, one coordinate per Bernoulli block. This diagonal
whitening does not remove cross-feature covariance when the raw Bernoulli
columns are dependent. The chi-square reference assumes the whitened null is
isotropic on the tested subspace and the projection is fixed independently of
the contrast; empirical-null inflation does not by itself demonstrate that
`R Cov(z) R.T = I` ([feature_space.py](../../tree_break_selection/tree/feature_space.py#L166),
[contrast_covariance.py](../../tree_break_selection/hierarchy_analysis/statistics/contrast_covariance.py#L420),
[projected_wald_reference_distribution.py](../../tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_reference_distribution.py#L6)).
For correlated binary annotations, this is an applicability and calibration
question, not a proven cause of passthrough-run outcomes.

There is also a ridge placement difference between per-observation tangent
whitening and Wald contrast covariance. Tangent whitening uses
`Sigma + ridge*I`; Wald covariance uses
`variance_scale*Sigma + ridge*I` ([contrast_covariance.py](../../tree_break_selection/hierarchy_analysis/statistics/contrast_covariance.py#L914),
[contrast_covariance.py](../../tree_break_selection/hierarchy_analysis/statistics/contrast_covariance.py#L867),
[contrast_covariance.py](../../tree_break_selection/hierarchy_analysis/statistics/contrast_covariance.py#L878)).
For positive ridge these do not generally differ by only the expected scalar
variance factor; with ridge zero they do. The practical difference is confined
to small-variance or nearly singular coordinates where the ridge matters. The
known-branch oracle diagnostic used zero ridge, so it does not test this case.
This is an implementation distinction requiring calibration for near-singular
regimes, not evidence that ridge placement caused observed passthrough p-value
changes.

## What this establishes

- **No reproducible normalization defect found in this bounded source audit.**
  The observed transformations are explicit and distinct: geometry
  normalization for tree construction, feature-family covariance whitening
  for TBS statistics, and `N`-scaled null-whitened covariance for PCA/MP.
- **Diffusion normalization differs from full pydiffmap.** The custom symmetric
  kernel path omits pydiffmap's density `q` correction and Markov generator.
  This is a geometry choice; whether it is intended for every application
  remains an applicability question, not a demonstrated failure.
- **Continuous covariance uses default `parent_total`.** Its between-child
  contribution can change Wald standard errors and p-values relative to the
  opt-in within-child covariance policy. That change comes from covariance
  policy/data, not passthrough itself.
- **Same-data PCA/rank selection remains a calibration limitation.** It is
  separate from representation normalization and was not addressed here.
- **Marginal Bernoulli whitening leaves dependence unmodeled.** Whether that
  matters for a specific correlated binary feature set requires null
  calibration under its joint feature distribution.
- **Ridge placement differs across per-observation and contrast covariance.**
  This can matter when ridge is material relative to covariance eigenvalues;
  this audit does not show a resulting p-value error.
- No arbitrary feature standardization is indicated by this audit: it would
  change the covariance model and could change the null reference law.
