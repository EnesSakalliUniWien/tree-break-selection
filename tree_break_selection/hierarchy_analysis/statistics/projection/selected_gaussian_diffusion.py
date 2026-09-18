"""Pointwise diffusion-hierarchy replay along a Gaussian radial path.

Adaptive diffusion refits its neighbor graph, kernel bandwidths, eigenspace,
and diffusion distances at every radius.  Hamming diffusion is defined only
where the reconstructed matrix remains binary or one-hot.  These replays
therefore diagnose the selected-hierarchy boundary but do not construct an
exact radial interval union or justify a selected chi-radial tail.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.cluster.hierarchy import linkage

from ....space_separation import (
    DiffusionGeometry,
    adaptive_diffusion_geometry,
    hamming_knn_diffusion_geometry,
)
from .selected_gaussian_hierarchy import HierarchyMergeSignature
from .selected_gaussian_radial import GaussianSelectedProjectionPath


@dataclass(frozen=True)
class DiffusionAverageLinkageReplay:
    """One numerical diffusion fit and average-linkage replay at a fixed radius."""

    radius: float
    distance_condensed: NDArray[np.float64]
    merge_signature: HierarchyMergeSignature
    geometry_metadata: dict[str, object]


def _reconstructed_matrix(
    path: GaussianSelectedProjectionPath,
    *,
    radius: float,
) -> tuple[float, NDArray[np.float64]]:
    if not isinstance(path, GaussianSelectedProjectionPath):
        raise TypeError("Diffusion replay requires a GaussianSelectedProjectionPath.")
    if path.nuisance_matrix.shape[0] < 2:
        raise ValueError("Diffusion hierarchy replay requires at least two leaves.")
    radius_value = float(radius)
    return radius_value, path.reconstruct(radius_value)


# Reconstructing the observed matrix through the whitened radial factors leaves
# roundoff near machine epsilon, while genuinely continuous radial candidates
# leave binary support by macroscopic margins.
_BINARY_SUPPORT_ABSOLUTE_TOLERANCE = 1e-9


def _require_binary_radial_candidate(
    feature_matrix: NDArray[np.float64],
) -> NDArray[np.float64]:
    rounded = np.rint(feature_matrix)
    if (
        not np.isin(rounded, (0.0, 1.0)).all()
        or float(np.abs(feature_matrix - rounded).max())
        > _BINARY_SUPPORT_ABSOLUTE_TOLERANCE
    ):
        raise ValueError(
            "Hamming diffusion radial replay requires binary or one-hot reconstructed "
            "values; continuous Gaussian radial candidates are unsupported."
        )
    return rounded


def _merge_signature_from_linkage(
    linkage_matrix: NDArray[np.float64],
    *,
    n_leaves: int,
) -> HierarchyMergeSignature:
    members_by_cluster: dict[int, tuple[int, ...]] = {
        leaf_index: (leaf_index,) for leaf_index in range(n_leaves)
    }
    signature = []
    for step_index, row in enumerate(linkage_matrix):
        left_id = int(row[0])
        right_id = int(row[1])
        left_members, right_members = sorted(
            (members_by_cluster[left_id], members_by_cluster[right_id])
        )
        signature.append((left_members, right_members))
        members_by_cluster[n_leaves + step_index] = tuple(
            sorted((*left_members, *right_members))
        )
    return tuple(signature)


def _replay_geometry(
    geometry: DiffusionGeometry,
    *,
    radius: float,
    n_leaves: int,
) -> DiffusionAverageLinkageReplay:
    linkage_matrix = linkage(geometry.distance_condensed, method="average")
    return DiffusionAverageLinkageReplay(
        radius=radius,
        distance_condensed=np.asarray(geometry.distance_condensed, dtype=np.float64),
        merge_signature=_merge_signature_from_linkage(
            linkage_matrix,
            n_leaves=n_leaves,
        ),
        geometry_metadata=dict(geometry.metadata),
    )


def replay_adaptive_pydiffmap_average_linkage(
    path: GaussianSelectedProjectionPath,
    *,
    radius: float,
    k_neighbors: int,
    diffusion_time: int,
    n_components: int,
    metric: str,
    bandwidth_type: str | float | None,
    epsilon: str | float,
) -> DiffusionAverageLinkageReplay:
    """Refit adaptive pydiffmap geometry and replay average linkage at one radius."""
    radius_value, feature_matrix = _reconstructed_matrix(path, radius=radius)
    if str(metric).lower() == "hamming":
        feature_matrix = _require_binary_radial_candidate(feature_matrix)
    geometry = adaptive_diffusion_geometry(
        feature_matrix,
        k_neighbors=k_neighbors,
        diffusion_time=diffusion_time,
        n_components=n_components,
        metric=metric,
        bandwidth_type=bandwidth_type,
        epsilon=epsilon,
    )
    return _replay_geometry(
        geometry,
        radius=radius_value,
        n_leaves=feature_matrix.shape[0],
    )


def replay_hamming_knn_average_linkage(
    path: GaussianSelectedProjectionPath,
    *,
    radius: float,
    k_neighbors: int,
    diffusion_time: int,
    n_components: int,
) -> DiffusionAverageLinkageReplay:
    """Refit Hamming-neighbor diffusion and replay average linkage at one radius."""
    radius_value, feature_matrix = _reconstructed_matrix(path, radius=radius)
    feature_matrix = _require_binary_radial_candidate(feature_matrix)
    geometry = hamming_knn_diffusion_geometry(
        feature_matrix,
        k_neighbors=k_neighbors,
        diffusion_time=diffusion_time,
        n_components=n_components,
    )
    return _replay_geometry(
        geometry,
        radius=radius_value,
        n_leaves=feature_matrix.shape[0],
    )


__all__ = [
    "DiffusionAverageLinkageReplay",
    "replay_adaptive_pydiffmap_average_linkage",
    "replay_hamming_knn_average_linkage",
]
