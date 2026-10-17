"""Synthetic continuous covariance profiles for calibration programs."""

from __future__ import annotations

from typing import Literal

import numpy as np

ContinuousCovarianceProfile = Literal["identity", "ar1_0.6", "ill_conditioned"]


def continuous_covariance(
    *,
    dimension: int,
    profile: ContinuousCovarianceProfile,
) -> np.ndarray:
    """Construct one supported synthetic continuous covariance profile."""

    if dimension <= 0:
        raise ValueError("dimension must be positive.")
    if profile == "identity":
        return np.eye(dimension, dtype=np.float64)
    if profile == "ar1_0.6":
        indices = np.arange(dimension)
        return 0.6 ** np.abs(indices[:, None] - indices[None, :])
    if profile == "ill_conditioned":
        eigenvalues = np.geomspace(1.0, 1e-3, num=dimension)
        return np.diag(eigenvalues)
    raise ValueError(f"Unsupported covariance profile: {profile!r}.")


__all__ = [
    "ContinuousCovarianceProfile",
    "continuous_covariance",
]
