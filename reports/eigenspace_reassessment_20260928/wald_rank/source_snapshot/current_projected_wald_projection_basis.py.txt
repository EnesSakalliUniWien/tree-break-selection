"""Projection basis construction from PCA components."""

from __future__ import annotations

import numpy as np

from tree_break_selection.eigenvalue_groups import extend_rank_through_ties


def _resolve_pca_component(
    pca_projection: np.ndarray,
    pca_eigenvalues: np.ndarray | None,
    target_dim: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Select at least target_dim rows without cutting a known numerical tie."""
    pca_basis = np.asarray(pca_projection, dtype=np.float64)
    if pca_basis.shape[0] < target_dim:
        raise ValueError(
            f"PCA projection has {pca_basis.shape[0]} row(s), but {target_dim} "
            "projection dimension(s) were requested."
        )
    if pca_eigenvalues is None:
        raise ValueError("PCA projection requires matching eigenvalues.")
    eigenvalues = np.asarray(pca_eigenvalues, dtype=np.float64)
    if eigenvalues.shape[0] < target_dim:
        raise ValueError(
            f"PCA eigenvalues have {eigenvalues.shape[0]} value(s), but {target_dim} "
            "projection dimension(s) were requested."
        )
    target_dim = extend_rank_through_ties(eigenvalues, target_dim)
    if pca_basis.shape[0] < target_dim:
        raise ValueError(
            f"PCA projection has {pca_basis.shape[0]} row(s), but {target_dim} "
            "projection dimension(s) are needed to retain the complete numerical tie group."
        )
    return pca_basis[:target_dim], eigenvalues[:target_dim]


def build_pca_projection_basis(
    k: int,
    *,
    pca_projection: np.ndarray | None = None,
    pca_eigenvalues: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Construct parent PCA rows, extending target k through numerical ties."""

    target_projection_dim = int(k)

    if pca_projection is None:
        raise ValueError("Projected Wald tests require a PCA projection basis.")

    return _resolve_pca_component(pca_projection, pca_eigenvalues, target_projection_dim)


__all__ = [
    "build_pca_projection_basis",
]
