"""Top-level sibling Wald test orchestration."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
from numpy.typing import NDArray

from tree_break_selection.tree.feature_space import FeatureSpace

from ....contrast_covariance import compute_whitened_wald_contrast
from ....projection.projected_wald.projected_wald_kernel import run_projected_wald_kernel
from ...projection.pair_testing.projection_dimension import (
    resolve_sibling_projection_dimension,
)


def sibling_divergence_test(
    left_distribution: np.ndarray,
    right_distribution: np.ndarray,
    left_sample_size: float,
    right_sample_size: float,
    *,
    branch_length_sum: float,
    mean_branch_length: float | None,
    projection_dimension_from_edge_comparisons: int | None = None,
    parent_principal_component_projection: np.ndarray | None = None,
    parent_principal_component_eigenvalues: np.ndarray | None = None,
    feature_space: FeatureSpace | None = None,
    continuous_covariance_by_block: Mapping[str, NDArray[np.floating]] | None = None,
    adaptive_projection_dimension_energy_fraction: float | None = None,
) -> tuple[float, float, float, float]:
    """Two-sample Wald test for sibling divergence.

    The contrast variance is always scaled by the sibling branch time,
    ``1 + branch_length_sum / (2 * mean_branch_length)``, whatever
    ``edge_branch_length_variance_policy`` says; that policy governs only the
    child-parent edge test. This is the only sibling test with a branch-time
    factor; the fixed-subspace gates use none. Pass ``branch_length_sum=0.0``
    for a tree without branch lengths (``mean_branch_length`` is then ``None``).
    """
    z_scores = compute_whitened_wald_contrast(
        left_distribution,
        right_distribution,
        float(left_sample_size),
        float(right_sample_size),
        comparison="sibling",
        feature_space=feature_space,
        continuous_covariance_by_block=continuous_covariance_by_block,
        tree_time=branch_length_sum,
        tree_time_normalizer=None
        if mean_branch_length is None
        else 2.0 * float(mean_branch_length),
    )

    if not np.isfinite(z_scores).all():
        raise ValueError(
            "Sibling Wald z-scores must be finite; "
            f"found {int(np.sum(~np.isfinite(z_scores)))} non-finite component(s)."
        )

    sibling_projection_dimension, _ = resolve_sibling_projection_dimension(
        projection_dimension_from_edge_comparisons=projection_dimension_from_edge_comparisons,
    )

    z_scores = z_scores.astype(np.float64, copy=False)

    result = run_projected_wald_kernel(
        z_scores,
        spectral_k=sibling_projection_dimension,
        pca_projection=parent_principal_component_projection,
        pca_eigenvalues=parent_principal_component_eigenvalues,
        adaptive_dimension_energy_fraction=(adaptive_projection_dimension_energy_fraction),
    )

    return (
        result.statistic,
        result.reference_scale,
        result.degrees_of_freedom,
        result.p_value,
    )


__all__ = ["sibling_divergence_test"]
