# statistics/

Statistical kernels and multiple-testing logic used by the decomposition gate
pipeline.

## Public Entrypoints

Import each entrypoint from the path shown. `child_parent_divergence`, `multiple_testing`,
`multiple_testing.tree_bh` and `projection.projected_wald` re-export their public names; the
`statistics`, `projection` and `sibling_divergence` package `__init__.py` files do not.

| Entrypoint | Purpose |
| ---------- | ------- |
| `child_parent_divergence.annotate_child_parent_divergence` | Annotate child-parent edge evidence for tree traversal. |
| `sibling_divergence.inflated_projected_wald_annotation.pipeline.annotate_sibling_divergence` | Annotate projected-Wald sibling evidence using empirical inflation and traversal-aligned BH; selected-tail validity remains open. |
| `sibling_divergence.fixed_subspace_annotation.annotate_fixed_subspace_sibling_divergence` | Opt-in sibling gate using fixed covariance-whitened coordinate or feature-block BH aggregation without parent PCA adaptation. |
| `multiple_testing.benjamini_hochberg_correction` | Flat BH correction helper. |
| `multiple_testing.apply_tree_bh_correction` | Tree-aware BH correction for hierarchical edge testing. |
| `independent_gaussian.fit_independent_gaussian_tree` | Fit topology and PCA on discovery measurements with known Gaussian covariance. |
| `independent_gaussian.test_independent_gaussian_tree` | Test independent measurements of the same leaves using the conditional energy reference, then existing corrections and traversal. |

## Directory Map

| Path | Purpose |
| ---- | ------- |
| `contrast_covariance.py` | Canonical contrast-covariance object for projected-Wald statistics. |
| `independent_gaussian.py` | Explicit discovery/test API for independent paired Gaussian measurements and externally known covariance; no same-data empirical inflation. |
| `projection/projected_wald/` | Shared projected-Wald kernel, projection basis, and reference distribution. |
| `projection/selected_gaussian_radial.py` | Internal continuous-Gaussian selected radial prototype with a nuisance-derived projection and exact declared selection intervals; not connected to production gates. |
| `projection/selected_gaussian_hierarchy.py` | Analytic observed merge-sequence interval construction for the internal squared-Euclidean average-linkage Gaussian prototype; not a production tree model. |
| `projection/selected_gaussian_diffusion.py` | Pointwise adaptive-pydiffmap and Hamming-neighbor average-linkage replay along the Gaussian radial path; diagnostic only, with no exact selected interval claim. |
| `projection/spectral/` | Node spectral tasks, Marchenko-Pastur thresholding, and tree-level spectral context. |
| `projection/projection_dimension_estimation/` | Projection dimension estimators. |
| `child_parent_divergence/` | Child-parent projected-Wald tests and tree-level edge annotation. |
| `sibling_divergence/` | Sibling pair collection, projected-Wald testing, fixed-subspace gates, empirical-null inflation, and sibling FDR. |
| `multiple_testing/` | Flat BH and tree-BH helpers. |
| `branch_length_utils.py` | Branch-length validation and normalized tree-time variance utilities. |
| `data_utils.py` | Node-attribute extraction and edge/sibling gate annotation-column initialization shared by the gate statistics. |

## Active Method Flow

1. Build node distributions on the `PosetTree`.
2. Build spectral/projection context for each node.
3. Run branch-time-inflated child-parent projected-Wald tests.
4. Run the configured sibling gate. The default path collects sibling pair
   records, estimates empirical-null inflation, and applies the historical
   plug-in chi-square tail before sibling FDR. The same-selected-hierarchy
   reference law remains unvalidated. The opt-in fixed-subspace path skips parent PCA and
   edge-derived sibling dimensions, then applies either the full fixed
   chi-square statistic or coordinate-wise/feature-block BH p-value
   aggregation before traversal-aligned sibling FDR. A positive
   `sibling_gate_alpha_penalty` divides the sibling FDR alpha before
   annotation, matching the selected-topology penalty used in diagnostics. An
   optional root-stability guard can then fail closed on an unstable open root
   split without changing non-root sibling decisions. Validation-only profiles
   can attach selected-permutation diagnostics with explicit scope: `root`,
   broad `open_internal`, narrower `passthrough_descendant` contexts, or the
   conservative `global_sibling_min_passthrough_descendant` selected-family
   correction. These diagnostics are for understanding selected-tree null
   behavior, not for defining the TBS runtime method. The refined global scope
   reruns selected-family Monte Carlo floor pass-through cases at higher
   resolution before accepting the diagnostic split.
   Pure Bernoulli and pure categorical fixed-coordinate BH p-values use the
   exact vectorized `compute_whitened_wald_contrast` path before the existing
   coordinate-wise BH aggregation.
5. Write sibling gate columns.
6. `TreeDecomposition` consumes the annotation bundle during top-down traversal.

## Contracts

- Feature-space covariance and contrast dimensions come from the active
  `FeatureSpace`; tests should not infer alternate schema names.
- Missing support closes sibling gates without substituting a neutral inflation
  value. With support, `internal_admissible` means the empirical rule is
  available; it does not certify the selected-tail reference law.
- Projection quantities and raw feature coordinates are separate objects.
- Spectral rank targets extend through complete numerical eigenvalue groups,
  using the shared `tree_break_selection/eigenvalue_groups.py` rule. The raw MP
  signal count is distinct from the extended test dimension, and Wald degrees
  of freedom use the returned basis dimension. A supplied PCA basis cannot
  silently omit known tied directions. This geometric requirement does not
  establish a null reference law for contrast-selected rank or selected trees.
- The selected Gaussian radial prototype supports a fixed scalar branch-time
  covariance multiplier through radius rescaling. Recomputed linkage/NNLS
  branch lengths and `branch_length_state` projections are not fixed by the
  current merge-signature event and remain outside its chi-radial law.
- Adaptive pydiffmap replay refits its data-dependent diffusion geometry at
  every radius. Hamming diffusion replay is defined only at binary or one-hot
  reconstructed matrices; a continuous Gaussian radial interval is therefore
  unsupported for that geometry.
- Fixed-subspace sibling gates are opt-in method choices, not automatic
  production promotion.
- Selected-permutation guard scopes are validation diagnostics, not TBS method
  choices. The
  `passthrough_descendant` scope targets descendant splits reached through an
  ordinary closed sibling ancestor and avoids compounding below explicitly
  guard-blocked roots. The
  `global_sibling_min_passthrough_descendant` scope evaluates those
  pass-through descendants against the minimum fixed-subspace sibling p-value
  over every binary parent in a fully reselected permutation null tree.
- DataFrame annotations are a transport format for gate columns; statistical
  kernels should use typed records where available.

## Independent Gaussian measurement contract

The independent API requires two measurement matrices with identical leaf and
feature labels, genuinely independent noise, independent rows, and a common
externally known positive-definite Gaussian covariance. Distinct observation
IDs represent distinct measurements of the same leaves. Renaming IDs or copying
one matrix cannot create independence; neither can splitting leaves into two
unmatched trees. Covariance estimated from these matrices is not a known input.

Discovery fits tree topology and positive PCA candidate directions. Inference
uses only the frozen discovery context and the separate measurement matrix.
`energy_selected_reference.py` conditions on contrast-energy selection through
complete numerical eigenvalue groups. Fixed full candidate rank uses chi-square;
adaptive rank uses a conditional mixture, not chi-square at a fitted df.
Retained rank is recorded separately from the reference-law name. The legacy
`Sibling_Degrees_of_Freedom` transport column carries retained rank and has an
explicit semantics column in this path.

Known covariance removes the empirical-inflation stage for this API. Both edge
and sibling p-values are replaced before existing Tree-BH, sibling BH and
`TreeDecomposition`. For binary groups with common known covariance, normalized
child-parent and sibling contrasts coincide up to sign. The API does not treat
fitted branch lengths as stochastic covariance times.

Ordinary single-matrix callers remain on their existing path because they lack
these inputs. Bernoulli/categorical data, covariance estimation, dependent paired
measurements and diffusion on binary data are outside this Gaussian contract.
Benchmark evidence is in `reports/production_pvalue_replacement_20260928/`.

For discovery-only adaptive Euclidean diffusion with fixed-topology NNLS, pass
`diffusion=IndependentGaussianDiffusion()` to `fit_independent_gaussian_tree`.
This reuses the existing pydiffmap geometry and NNLS implementation, requires a
linkage tree, and fits both only on original discovery measurements. Set
`fit_nnls=False` on the configuration to compare diffusion without NNLS.
Topology, branch lengths, diffusion coordinates and PCA are frozen before
testing. Solver nonconvergence is rejected; solver diagnostics are retained.

Under the independent known-covariance model, NNLS lengths have no covariance
role, so diffusion and diffusion plus NNLS give identical tests on a fixed
topology. This differs from the current single-matrix branch-time model.
The matched 900-pair benchmark, root tails, solver diagnostics and limitations
are in `reports/independent_diffusion_nnls_20260928/`.
