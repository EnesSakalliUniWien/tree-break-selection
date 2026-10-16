"""Mass-weighted distributional movement diagnostics for tree edges."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike


@dataclass(frozen=True)
class EdgeDistributionalActionSummary:
    parent_mass: float
    child_mass: float
    child_parent_mass_fraction: float
    squared_displacement: float
    action: float


def _as_1d_float_vector(value: ArrayLike, *, name: str) -> np.ndarray:
    vector = np.asarray(value, dtype=float)
    if vector.ndim != 1:
        raise ValueError(f"{name} must be a one-dimensional vector.")
    if not np.all(np.isfinite(vector)):
        raise ValueError(f"{name} must contain only finite values.")
    return vector


def _as_positive_count(value: int | float, *, name: str) -> float:
    count = float(value)
    if not np.isfinite(count) or count <= 0.0:
        raise ValueError(f"{name} must be positive and finite.")
    return count


def edge_distributional_action_summary(
    parent_mean: ArrayLike,
    child_mean: ArrayLike,
    parent_leaf_count: int | float,
    child_leaf_count: int | float,
) -> EdgeDistributionalActionSummary:
    """Return edge action together with parent and child distribution masses."""
    parent = _as_1d_float_vector(parent_mean, name="parent_mean")
    child = _as_1d_float_vector(child_mean, name="child_mean")
    if parent.shape != child.shape:
        raise ValueError("parent_mean and child_mean must have the same shape.")
    parent_count = _as_positive_count(parent_leaf_count, name="parent_leaf_count")
    child_count = _as_positive_count(child_leaf_count, name="child_leaf_count")
    delta = child - parent
    squared_displacement = float(np.dot(delta, delta))
    return EdgeDistributionalActionSummary(
        parent_mass=parent_count,
        child_mass=child_count,
        child_parent_mass_fraction=float(child_count / parent_count),
        squared_displacement=squared_displacement,
        action=float(child_count * squared_displacement),
    )
