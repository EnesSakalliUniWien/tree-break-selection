"""Typed tree-level output for local spectral decomposition."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class SpectralDecompositionResult:
    """Per-node spectral contracts consumed by projected-Wald gates."""

    test_projection_dimensions_by_node: dict[str, int]
    raw_mp_signal_counts_by_node: dict[str, int]
    effective_independent_rows_by_node: dict[str, int]
    mp_threshold_rows_by_node: dict[str, int]
    principal_component_projections_by_node: dict[str, np.ndarray]
    principal_component_eigenvalues_by_node: dict[str, np.ndarray]
    descendant_leaf_row_counts_by_node: dict[str, int]
    internal_distribution_row_counts_by_node: dict[str, int]
    spectral_matrix_row_counts_by_node: dict[str, int]
    full_component_eigenvalues_by_node: dict[str, np.ndarray]
    active_feature_counts_by_node: dict[str, int]
    stage_timings: dict[str, float]


__all__ = ["SpectralDecompositionResult"]
