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


@dataclass(frozen=True)
class SplitDistributionalActionSummary:
    child_masses: tuple[float, ...]
    parent_mass: float
    parent_mean: np.ndarray
    child_parent_mass_fractions: tuple[float, ...]
    squared_child_parent_displacements: tuple[float, ...]
    child_edge_actions: tuple[float, ...]
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


def split_distributional_action_summary(
    child_means: list[ArrayLike] | tuple[ArrayLike, ...],
    child_masses: list[int | float] | tuple[int | float, ...],
) -> SplitDistributionalActionSummary:
    """Return the barycentric split action for one internal node.

    This is the between-child term in the parallel-axis decomposition:
    ``sum_i m_i ||mu_i - mu_parent||^2``.
    """
    if len(child_means) != len(child_masses):
        raise ValueError("child_means and child_masses must have the same length.")
    if len(child_means) == 0:
        raise ValueError("At least one child distribution is required.")
    means = tuple(
        _as_1d_float_vector(mean, name=f"child_means[{index}]")
        for index, mean in enumerate(child_means)
    )
    first_shape = means[0].shape
    if any(mean.shape != first_shape for mean in means):
        raise ValueError("All child means must have the same shape.")
    masses = tuple(
        _as_positive_count(mass, name=f"child_masses[{index}]")
        for index, mass in enumerate(child_masses)
    )
    parent_mass = float(sum(masses))
    parent_mean = sum(mass * mean for mass, mean in zip(masses, means, strict=True))
    parent_mean = np.asarray(parent_mean, dtype=float) / parent_mass
    squared_displacements = tuple(
        float(np.dot(mean - parent_mean, mean - parent_mean)) for mean in means
    )
    edge_actions = tuple(
        float(mass * squared_displacement)
        for mass, squared_displacement in zip(masses, squared_displacements, strict=True)
    )
    return SplitDistributionalActionSummary(
        child_masses=masses,
        parent_mass=parent_mass,
        parent_mean=parent_mean,
        child_parent_mass_fractions=tuple(mass / parent_mass for mass in masses),
        squared_child_parent_displacements=squared_displacements,
        child_edge_actions=edge_actions,
        action=float(sum(edge_actions)),
    )
