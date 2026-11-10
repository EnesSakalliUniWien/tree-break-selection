"""Matched controls for the fitted raw native continuous sibling statistic."""

from numbers import Integral

import numpy as np

from tree_break_selection.eigenvalue_groups import tied_eigenvalue_neighbors
from tree_break_selection.hierarchy_analysis.statistics.branch_length_utils import (
    compute_mean_branch_length,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.pair_testing.collection.pair_observations import (
    extract_sibling_pair_observations,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.pair_testing.wald_statistic.sibling_divergence_test import (
    sibling_divergence_test,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.projection.pair_testing.parent_projection_resolution import (
    resolve_parent_projection_inputs_for_sibling_test,
)
from tree_break_selection.tree.distributions import resolve_node_continuous_covariance_by_block

from .subtree_topology_tests import _native_fit_record


def native_sibling_attribution(values, subtree_indices, *, oracle_covariance=None, gram_rank=None):
    """Refit native components, then vary one stated part of the raw score.

    Unit branch holds the fitted tree, covariance and projection context fixed.
    Oracle covariance uses the supplied externally fixed positive-definite
    matrix in the same native context. Fixed Gram rank replaces the native
    projection by the leading raw subtree Gram directions at a prespecified
    rank while retaining the guarded covariance and branch factor. These are
    score controls; neither projected variant implies an analytic null law.
    """
    data, space, fitted, dimensions, spectral, record = _native_fit_record(values, subtree_indices)
    tree = fitted.extra["tree"]
    projection = resolve_parent_projection_inputs_for_sibling_test(
        record.parent,
        sibling_projection_dimensions_from_edge_comparisons=dimensions,
        parent_principal_component_projections=spectral.principal_component_projections_by_node,
        parent_principal_component_eigenvalues=spectral.principal_component_eigenvalues_by_node,
    )
    contrast = extract_sibling_pair_observations(tree, record.parent, record.left, record.right)[:4]
    guarded = resolve_node_continuous_covariance_by_block(tree, record.parent, space)
    mean = compute_mean_branch_length(tree)
    g = 1.0 if mean is None else 1.0 + record.branch_length_sum / (2 * mean)
    annotations = fitted.extra["annotations"]
    annotated = annotations.loc[record.parent]
    inflation = float(annotated.get("Sibling_Inflation_Factor", np.nan))
    has_display = np.isfinite(inflation) and inflation > 0

    def raw(*, branch=record.branch_length_sum, covariance=guarded, context=projection, adapt=0.90):
        k, basis, eig = context
        return sibling_divergence_test(
            *contrast,
            branch_length_sum=branch,
            mean_branch_length=mean,
            feature_space=space,
            projection_dimension_from_edge_comparisons=k,
            parent_principal_component_projection=basis,
            parent_principal_component_eigenvalues=eig,
            continuous_covariance_by_block=covariance,
            adaptive_projection_dimension_energy_fraction=adapt,
        )

    unit, _, unit_df, _ = raw(branch=0.0)
    result = dict(
        native_raw=float(record.stat),
        native_df=float(record.degrees_of_freedom),
        native_display=float(record.stat / (record.reference_scale * inflation))
        if has_display
        else None,
        annotated_display=float(annotated["Sibling_Test_Statistic"]) if has_display else None,
        scale=float(record.reference_scale),
        inflation=inflation if has_display else None,
        unit_branch_raw=float(unit),
        unit_branch_df=float(unit_df),
        branch_factor=float(g),
        sibling_branch_sum=float(record.branch_length_sum),
        positive_branch_mean=mean,
        child_sizes=[int(contrast[2]), int(contrast[3])],
    )
    if oracle_covariance is not None:
        oracle = np.asarray(oracle_covariance, dtype=float)
        p = data.shape[1]
        if (
            oracle.shape != (p, p)
            or not np.isfinite(oracle).all()
            or not np.allclose(oracle, oracle.T)
        ):
            raise ValueError("oracle_covariance must be a symmetric feature covariance matrix")
        np.linalg.cholesky(oracle)
        name = space.continuous_blocks[0].name
        oracle_raw, _, oracle_df, _ = raw(covariance={name: oracle})
        result.update(
            oracle_covariance_raw=float(oracle_raw), oracle_covariance_df=float(oracle_df)
        )
    if gram_rank is not None:
        if isinstance(gram_rank, bool) or not isinstance(gram_rank, Integral) or gram_rank < 1:
            raise ValueError("gram_rank must be a positive integer")
        clades = tree.compute_descendant_sets()
        x = data.loc[sorted(clades[record.parent])].to_numpy()
        _, singular, basis = np.linalg.svd(x - x.mean(axis=0), full_matrices=False)
        eigenvalues = singular**2
        available = int(
            np.count_nonzero(singular > singular[0] * max(x.shape) * np.finfo(float).eps)
        )
        if gram_rank > available or (
            gram_rank < available and tied_eigenvalue_neighbors(eigenvalues)[gram_rank - 1]
        ):
            raise ValueError("gram_rank exceeds the positive rank or cuts a tied eigenspace")
        gram_raw, _, gram_df, _ = raw(
            context=(gram_rank, basis[:gram_rank], eigenvalues[:gram_rank]), adapt=None
        )
        result.update(gram_rank_raw=float(gram_raw), gram_rank_df=float(gram_df))
    return result
