# Compatibility with the independent equations

## Internal distributions

`tree_estimator.py::_precompute_branch_length_internal_states` recursively uses
weights `descendant_leaf_count / max(edge_length, floor)`. The floor is
`max(median_positive_length * 1e-6, machine_epsilon)` when positive lengths
exist. These are diagnostic internal states, not an established ancestral
posterior distribution. The stored node `distribution` stays the empirical
leaf-count barycenter. `spectral_include_internal_barycenters=True` is needed
to add these states to the spectral matrix; the mode string alone does nothing.

In the current method the added states are dependent combinations of the same
observations used in testing. Keeping the MP row count equal to the number of
leaves avoids counting invented independent observations, but does not prove
an MP reference for this augmented covariance or calibrate the selected test.

In the independent diagnostic, let D be discovery measurements, Y independent
inference measurements, Sigma externally known feature covariance and Q(D) an
orthonormal candidate projection. The topology, lengths and internal states
are all functions of D. Under a null sibling contrast,

    z = Sigma^(-1/2) (mean(Y_A) - mean(Y_B)) / sqrt(1/n_A + 1/n_B)

has standard Gaussian feature coordinates conditional on D. Therefore Q(D)z
has independent standard Gaussian candidate coordinates. Adding discovery
internal rows changes Q and its ordered spectrum without adding test degrees
of freedom. The report adapter retains the leaf-only candidate-rank cap,
checks orthonormality, and applies the existing complete-group energy selection
and conditional tail to independent inference measurements.

This establishes compatibility under the stated Gaussian model. It does not
prove whole-tree FDR, validate estimated covariance, or make the current
same-data pipeline independent. The adapter uses known-covariance whitening,
not the current method's estimated local covariance and MP rank selection.

## Branch variance

The edge flag is specifically `edge_branch_length_variance_policy`. When it is
`none`, sibling branch variance is still present. Current sibling covariance
uses `m = 1 + (b_left + b_right)/(2 * mean_positive_branch_length)`; the enabled
edge policy uses `m_edge = 1 + b_edge/mean_positive_branch_length`.

For the independent-row generator the true contrast covariance is fixed by
Sigma and leaf counts. Fitted geometric branch lengths supply no additional
noise in this generator. Conditionally on discovery, let the energy-selected
statistic be T with the existing conditional survival function S. Adding a
scalar multiplier m to its variance gives T_new = T/m, while uniform scaling
preserves the energy-selected rank and complete-group selection event. Hence

    P(T_new >= t | discovery, selected rank) = S(m * t).

At the observed statistic, the correct tail is S(m * T_new) = S(T). Using
S(T_new) is a scale mismatch, not a degrees-of-freedom correction. Since m >= 1,
that uncorrected tail is at least S(T). The benchmark explicitly records this
negative control for sibling multipliers and verifies restoration to numerical
precision. The same algebra applies separately to either edge multiplier.
No full clustering run treats these incompatible tails as valid new equations.

To justify a different branch-dependent reference, a stochastic covariance
model for the measurements and branch lengths is still needed. Independence
of fitting data alone does not turn a geometric length into that covariance.
The older `53_test_selected_gaussian_branch_lengths.py` counterexample also
remains relevant: holding a merge signature fixed does not hold recomputed
branch-length-dependent projections fixed in a same-data radial construction.

## Source locations

- `tree_break_selection/hierarchy_analysis/statistics/projection/spectral/tree_estimator.py`
- `tree_break_selection/hierarchy_analysis/statistics/contrast_covariance.py`
- `tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/pair_testing/collection/record_collection.py`
- `tree_break_selection/hierarchy_analysis/statistics/child_parent_divergence/child_parent_divergence_annotation/tree_testing.py`
- `tree_break_selection/hierarchy_analysis/statistics/independent_gaussian.py`
- `tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/energy_selected_reference.py`
- `tests/statistics/53_test_selected_gaussian_branch_lengths.py`
