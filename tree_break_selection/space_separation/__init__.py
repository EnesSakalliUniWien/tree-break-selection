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
    DiffusionGeometry,
    adaptive_diffusion_geometry,
    block_adaptive_diffusion_geometry,
    block_diffusion_geometry,
    compute_diffusion_coordinates,
    hamming_knn_diffusion_geometry,
    resolve_adaptive_epsilon,
    resolve_neighbor_search_k,
)
from .leading_axes import decompose_leading_axes_space

__all__ = [
    "decompose_leading_axes_space",
    "DiffusionGeometry",
    "SpectralBlock",
    "adaptive_spectral_blocks",
    "adaptive_diffusion_geometry",
    "block_adaptive_diffusion_geometry",
    "block_diffusion_geometry",
    "compute_diffusion_coordinates",
    "coordinates_for_block",
    "cosine_eigendecomposition",
    "hamming_knn_diffusion_geometry",
    "separate_adaptive_cosine_space",
    "resolve_adaptive_epsilon",
    "resolve_neighbor_search_k",
    "weight_feature_matrix",
]
