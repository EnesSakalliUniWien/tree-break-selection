"""Reusable space-separation methods used across applications and benchmarks."""

from .adaptive_cosine import (
    SpectralBlock,
    adaptive_spectral_blocks,
    coordinates_for_block,
    cosine_eigendecomposition,
    separate_adaptive_cosine_space,
    weight_feature_matrix,
)
from .diffusion import (
    FEATURE_FAMILY_DIFFUSION_SETTINGS,
    FEATURE_FAMILY_DIFFUSION_TREE_DISTANCE_METRIC,
    DiffusionGeometry,
    FeatureFamilyDiffusionSettings,
    adaptive_diffusion_geometry,
    block_adaptive_diffusion_geometry,
    block_diffusion_geometry,
    compute_diffusion_coordinates,
    feature_family_diffusion_geometry,
    feature_family_diffusion_settings,
    feature_family_diffusion_unsupported_reason,
    hamming_knn_diffusion_geometry,
    resolve_adaptive_epsilon,
    resolve_neighbor_search_k,
)
from .leading_axes import decompose_leading_axes_space

__all__ = [
    "decompose_leading_axes_space",
    "DiffusionGeometry",
    "FEATURE_FAMILY_DIFFUSION_SETTINGS",
    "FEATURE_FAMILY_DIFFUSION_TREE_DISTANCE_METRIC",
    "FeatureFamilyDiffusionSettings",
    "SpectralBlock",
    "adaptive_spectral_blocks",
    "adaptive_diffusion_geometry",
    "block_adaptive_diffusion_geometry",
    "block_diffusion_geometry",
    "compute_diffusion_coordinates",
    "coordinates_for_block",
    "feature_family_diffusion_geometry",
    "feature_family_diffusion_settings",
    "feature_family_diffusion_unsupported_reason",
    "cosine_eigendecomposition",
    "hamming_knn_diffusion_geometry",
    "separate_adaptive_cosine_space",
    "resolve_adaptive_epsilon",
    "resolve_neighbor_search_k",
    "weight_feature_matrix",
]
