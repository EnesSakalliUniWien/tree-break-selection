"""Spectral context preparation for edge-divergence annotation."""

from __future__ import annotations

from dataclasses import replace
from time import perf_counter

import networkx as nx
import numpy as np
import pandas as pd

from tree_break_selection.tree.feature_space import FeatureSpace

from ...projection.spectral.spectral_decomposition_result import SpectralDecompositionResult
from ...projection.spectral.tree_estimator import (
    INTERNAL_DISTRIBUTION_EMPIRICAL_BARYCENTER,
    compute_spectral_decomposition,
)
from ...projection.spectral.tree_helpers import precompute_descendants

EDGE_GATE_SPECTRAL_MINIMUM_PROJECTION_DIMENSION = 2


def _validate_spectral_context_outputs(
    test_projection_dimensions_by_node: dict[str, int],
    principal_component_projections_by_node: dict[str, np.ndarray],
    principal_component_eigenvalues_by_node: dict[str, np.ndarray],
) -> None:
    """Validate paired edge-gate PCA outputs before sibling testing."""
    projection_nodes = set(principal_component_projections_by_node)
    eigenvalue_nodes = set(principal_component_eigenvalues_by_node)
    if projection_nodes != eigenvalue_nodes:
        missing_eigenvalues = sorted(projection_nodes - eigenvalue_nodes)
        missing_projections = sorted(eigenvalue_nodes - projection_nodes)
        raise ValueError(
            "Edge-gate spectral context requires matching PCA projection/eigenvalue node keys. "
            f"Missing eigenvalues for {missing_eigenvalues}; "
            f"missing projections for {missing_projections}."
        )

    dimension_nodes = set(test_projection_dimensions_by_node)
    missing_dimensions = sorted(projection_nodes - dimension_nodes)
    if missing_dimensions:
        raise ValueError(
            "Edge-gate spectral context has PCA outputs without test projection dimensions for "
            f"node(s): {missing_dimensions}."
        )
    missing_projection_outputs = sorted(
        node_id
        for node_id, projection_dimension in test_projection_dimensions_by_node.items()
        if int(projection_dimension) > 0 and node_id not in projection_nodes
    )
    if missing_projection_outputs:
        raise ValueError(
            "Edge-gate spectral context has positive test projection dimensions without PCA outputs "
            f"for node(s): {missing_projection_outputs}."
        )

    for node_id in sorted(projection_nodes):
        projection = np.asarray(principal_component_projections_by_node[node_id])
        eigenvalues = np.asarray(principal_component_eigenvalues_by_node[node_id])
        if projection.ndim != 2:
            raise ValueError(
                f"Edge-gate PCA projection for node {node_id!r} must be a 2-D matrix; "
                f"got shape {projection.shape}."
            )
        if eigenvalues.ndim != 1:
            raise ValueError(
                f"Edge-gate PCA eigenvalues for node {node_id!r} must be a 1-D vector; "
                f"got shape {eigenvalues.shape}."
            )
        if projection.shape[0] != eigenvalues.shape[0]:
            raise ValueError(
                f"Edge-gate PCA projection/eigenvalue row count mismatch for node {node_id!r}: "
                f"{projection.shape[0]} projection row(s), {eigenvalues.shape[0]} eigenvalue(s)."
            )
        if projection.shape[0] < int(test_projection_dimensions_by_node[node_id]):
            raise ValueError(
                f"Edge-gate PCA projection row count for node {node_id!r} must be at "
                "least its test projection dimension."
            )


def compute_child_parent_spectral_context(
    tree: nx.DiGraph,
    leaf_data: pd.DataFrame,
    *,
    minimum_projection_dimension: int = EDGE_GATE_SPECTRAL_MINIMUM_PROJECTION_DIMENSION,
    projection_basis_dimension: int | None = None,
    feature_space: FeatureSpace | None = None,
    include_internal_barycenters: bool = False,
    internal_distribution_mode: str = INTERNAL_DISTRIBUTION_EMPIRICAL_BARYCENTER,
    basis_row_indices_by_node: dict[str, list[int]] | None = None,
) -> SpectralDecompositionResult:
    """Prepare Marchenko-Pastur spectral context for the edge gate."""
    start_sec = perf_counter()
    spectral_decomposition = compute_spectral_decomposition(
        tree,
        leaf_data,
        feature_space=feature_space,
        minimum_projection_dimension=int(minimum_projection_dimension),
        projection_basis_dimension=projection_basis_dimension,
        include_internal_barycenters=bool(include_internal_barycenters),
        internal_distribution_mode=str(internal_distribution_mode),
        basis_row_indices_by_node=basis_row_indices_by_node,
    )

    _validate_spectral_context_outputs(
        spectral_decomposition.test_projection_dimensions_by_node,
        spectral_decomposition.principal_component_projections_by_node,
        spectral_decomposition.principal_component_eigenvalues_by_node,
    )
    stage_timings = dict(spectral_decomposition.stage_timings)
    stage_timings["spectral_context_sec"] = float(perf_counter() - start_sec)

    return replace(spectral_decomposition, stage_timings=stage_timings)


def split_edge_test_rows(
    tree: nx.DiGraph,
    leaf_data: pd.DataFrame,
) -> tuple[dict[str, list[int]], dict[str, list[int]]]:
    """Split every child's leaf rows into a test half and a basis half.

    A PCA basis fitted to the leaves whose contrast it projects is not fixed:
    it turns toward that contrast and breaks the chi-square law of the edge
    test. Each child's rows, in leaf-label order, alternate between halves,
    starting with the test half, so a single-leaf child is always tested. A
    parent's basis rows are its children's basis halves. The split ignores
    the data values and the row order of leaf_data, so under the null the
    basis is independent of the test-half contrast.

    Returns
    -------
    tuple
        Test rows by child node, and basis rows by parent node.
    """
    leaf_rows_by_node = precompute_descendants(
        tree, {label: row for row, label in enumerate(leaf_data.index)}
    )
    test_rows_by_child: dict[str, list[int]] = {}
    basis_rows_by_parent: dict[str, list[int]] = {}
    for parent_id, child_id in tree.edges():
        child_rows = sorted(leaf_rows_by_node[child_id], key=lambda row: str(leaf_data.index[row]))
        test_rows_by_child[child_id] = child_rows[0::2]
        basis_rows_by_parent.setdefault(parent_id, []).extend(child_rows[1::2])
    return test_rows_by_child, basis_rows_by_parent


__all__ = [
    "EDGE_GATE_SPECTRAL_MINIMUM_PROJECTION_DIMENSION",
    "compute_child_parent_spectral_context",
    "split_edge_test_rows",
]
