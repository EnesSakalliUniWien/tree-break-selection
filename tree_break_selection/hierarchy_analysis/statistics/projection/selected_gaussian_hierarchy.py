"""Exact squared-Euclidean average-linkage cells along a Gaussian radial path.

For this explicitly limited hierarchy model, every pairwise squared distance is
quadratic in the radial coordinate. Average-linkage updates preserve that
quadratic form, so the event selecting one complete observed merge sequence is
an intersection of finitely many quadratic inequalities. This module solves
those inequalities analytically and certifies the resulting interval interiors
by replay. It does not model SciPy average linkage over ordinary Euclidean
distances or any production tree-distance transformation.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext
from fractions import Fraction

import numpy as np
from numpy.typing import NDArray

from .selected_gaussian_radial import (
    GaussianSelectedProjectionPath,
    RadialSelectionInterval,
)

LeafCluster = tuple[int, ...]
HierarchyMerge = tuple[LeafCluster, LeafCluster]
HierarchyMergeSignature = tuple[HierarchyMerge, ...]
DistancePolynomial = NDArray[np.object_]


@dataclass(frozen=True)
class SquaredEuclideanAverageLinkageRegion:
    """Observed merge-sequence cell on the nonnegative Gaussian radial path."""

    observed_merge_signature: HierarchyMergeSignature
    selection_intervals: tuple[RadialSelectionInterval, ...]
    constraint_count: int
    nonnegative_boundary_radii: tuple[float, ...]


@dataclass(frozen=True)
class _HierarchyRun:
    merge_signature: HierarchyMergeSignature
    selection_constraints: tuple[DistancePolynomial, ...]


def _validate_path(path: GaussianSelectedProjectionPath) -> int:
    if not isinstance(path, GaussianSelectedProjectionPath):
        raise TypeError(
            "Squared-Euclidean average-linkage replay requires a "
            "GaussianSelectedProjectionPath."
        )
    n_leaves = int(path.nuisance_matrix.shape[0])
    if n_leaves < 2:
        raise ValueError("Hierarchy replay requires at least two leaves.")
    return n_leaves


def _pair_key(left: int, right: int) -> tuple[int, int]:
    return (left, right) if left < right else (right, left)


def _initial_squared_distance_polynomials(
    path: GaussianSelectedProjectionPath,
) -> dict[tuple[int, int], DistancePolynomial]:
    n_leaves = _validate_path(path)
    # Preserve the supplied binary floating-point factors exactly while forming
    # coefficients. Extra precision only at the root solver is too late: tiny
    # squared distances can already have vanished in coefficient subtraction.
    rational = np.vectorize(lambda value: Fraction(float(value)), otypes=[object])
    unit_contrast = rational(path.contrast) / Fraction(path.contrast_norm)
    covariance_factor = rational(path.covariance_factor)
    raw_fixed_contrast = covariance_factor @ rational(path.orthogonal_whitened_contrast)
    raw_radial_direction = (
        covariance_factor @ rational(path.projection.T) @ rational(path.projected_direction)
    )
    path_at_zero = rational(path.nuisance_matrix) + np.outer(
        unit_contrast,
        raw_fixed_contrast,
    )
    # Every row shares this direction, scaled only by its unit-contrast value.
    # Constructing it directly preserves exact slope-energy cancellations that
    # reconstruct(1) - reconstruct(0) can perturb at machine precision.
    radial_direction_squared_norm = raw_radial_direction @ raw_radial_direction
    polynomials: dict[tuple[int, int], DistancePolynomial] = {}
    for left in range(n_leaves):
        for right in range(left + 1, n_leaves):
            intercept_difference = path_at_zero[left] - path_at_zero[right]
            contrast_difference = unit_contrast[left] - unit_contrast[right]
            slope_difference = contrast_difference * raw_radial_direction
            polynomials[(left, right)] = np.array(
                [
                    intercept_difference @ intercept_difference,
                    2 * (intercept_difference @ slope_difference),
                    contrast_difference
                    * contrast_difference
                    * radial_direction_squared_norm,
                ],
                dtype=object,
            )
    return polynomials


def _evaluate_polynomial(polynomial: DistancePolynomial, radius: float) -> Fraction:
    exact_radius = Fraction(float(radius))
    return (polynomial[2] * exact_radius + polynomial[1]) * exact_radius + polynomial[0]


def _ordered_member_pair(
    left_members: LeafCluster,
    right_members: LeafCluster,
) -> HierarchyMerge:
    if left_members <= right_members:
        return left_members, right_members
    return right_members, left_members


def _run_squared_euclidean_average_linkage(
    path: GaussianSelectedProjectionPath,
    *,
    radius: float,
    collect_constraints: bool,
) -> _HierarchyRun:
    n_leaves = _validate_path(path)
    radius_value = float(radius)
    if not np.isfinite(radius_value) or radius_value < 0.0:
        raise ValueError(
            "Hierarchy replay radius must be finite and non-negative; "
            f"got {radius!r}."
        )

    active_clusters: dict[int, LeafCluster] = {
        leaf_id: (leaf_id,) for leaf_id in range(n_leaves)
    }
    distances = _initial_squared_distance_polynomials(path)
    merge_signature: list[HierarchyMerge] = []
    constraints: list[DistancePolynomial] = []

    for step_index in range(n_leaves - 1):
        candidate_pairs = sorted(
            pair
            for pair in distances
            if pair[0] in active_clusters and pair[1] in active_clusters
        )
        if not candidate_pairs:
            raise ValueError(f"Hierarchy replay has no active pair at step {step_index}.")
        selected_pair = min(
            candidate_pairs,
            key=lambda pair: (_evaluate_polynomial(distances[pair], radius_value), pair),
        )
        selected_polynomial = distances[selected_pair]
        if collect_constraints:
            constraints.extend(
                selected_polynomial - distances[competitor]
                for competitor in candidate_pairs
                if competitor != selected_pair
            )

        left, right = selected_pair
        left_members = active_clusters[left]
        right_members = active_clusters[right]
        merge_signature.append(_ordered_member_pair(left_members, right_members))
        new_members = tuple(sorted((*left_members, *right_members)))
        new_cluster_id = n_leaves + step_index

        other_cluster_ids = [
            cluster_id for cluster_id in active_clusters if cluster_id not in selected_pair
        ]
        for other in other_cluster_ids:
            left_distance = distances[_pair_key(left, other)]
            right_distance = distances[_pair_key(right, other)]
            distances[_pair_key(new_cluster_id, other)] = (
                len(left_members) * left_distance + len(right_members) * right_distance
            ) / len(new_members)

        for pair in list(distances):
            if left in pair or right in pair:
                del distances[pair]
        del active_clusters[left]
        del active_clusters[right]
        active_clusters[new_cluster_id] = new_members

    return _HierarchyRun(
        merge_signature=tuple(merge_signature),
        selection_constraints=tuple(constraints),
    )


def replay_squared_euclidean_average_linkage(
    path: GaussianSelectedProjectionPath,
    *,
    radius: float,
) -> HierarchyMergeSignature:
    """Replay deterministic squared-Euclidean average linkage at one radius."""
    return _run_squared_euclidean_average_linkage(
        path,
        radius=radius,
        collect_constraints=False,
    ).merge_signature


def _nonnegative_real_roots(
    polynomial: DistancePolynomial,
) -> tuple[float, ...]:
    constant, linear, quadratic = polynomial
    if quadratic == 0:
        if linear == 0:
            return ()
        root = -constant / linear
        return (float(root),) if root >= 0 else ()

    discriminant = linear * linear - 4 * quadratic * constant
    if discriminant < 0:
        return ()
    if discriminant == 0:
        root = -linear / (2 * quadratic)
        return (float(root),) if root >= 0 else ()

    def decimal(value: Fraction) -> Decimal:
        return Decimal(value.numerator) / Decimal(value.denominator)

    with localcontext() as context:
        context.prec = 80
        square_root = decimal(discriminant).sqrt()
        signed_root = square_root if linear >= 0 else -square_root
        stable_numerator = -(decimal(linear) + signed_root) / 2
        candidates = (
            stable_numerator / decimal(quadratic),
            decimal(constant) / stable_numerator,
        )
        return tuple(sorted(float(root) for root in candidates if root >= 0))


def _unique_sorted_radii(radii: list[float]) -> tuple[float, ...]:
    unique: list[float] = []
    for radius in sorted(radii):
        if unique and radius == unique[-1]:
            continue
        unique.append(float(radius))
    return tuple(unique)


def _constraint_holds(
    polynomial: DistancePolynomial,
    radius: float,
) -> bool:
    return bool(_evaluate_polynomial(polynomial, radius) <= 0.0)


def _all_constraints_hold(
    constraints: tuple[DistancePolynomial, ...],
    radius: float,
) -> bool:
    return all(_constraint_holds(polynomial, radius) for polynomial in constraints)


def _selected_intervals_from_constraints(
    constraints: tuple[DistancePolynomial, ...],
    boundary_radii: tuple[float, ...],
) -> tuple[RadialSelectionInterval, ...]:
    positive_boundaries = tuple(radius for radius in boundary_radii if radius > 0.0)
    finite_edges = (0.0, *positive_boundaries)
    selected: list[RadialSelectionInterval] = []
    for lower, upper in zip(finite_edges, finite_edges[1:]):
        if upper <= lower:
            continue
        midpoint = lower + 0.5 * (upper - lower)
        if _all_constraints_hold(constraints, midpoint):
            selected.append(RadialSelectionInterval(lower, upper))

    tail_lower = finite_edges[-1]
    tail_probe = tail_lower + max(1.0, abs(tail_lower))
    if _all_constraints_hold(constraints, tail_probe):
        selected.append(RadialSelectionInterval(tail_lower, np.inf))

    merged: list[RadialSelectionInterval] = []
    for interval in selected:
        if merged and merged[-1].upper == interval.lower:
            merged[-1] = RadialSelectionInterval(merged[-1].lower, interval.upper)
        else:
            merged.append(interval)
    return tuple(merged)


def _interval_probe(interval: RadialSelectionInterval) -> float:
    if interval.upper == np.inf:
        return float(interval.lower + max(1.0, abs(interval.lower)))
    return float(interval.lower + 0.5 * (interval.upper - interval.lower))


def construct_squared_euclidean_average_linkage_region(
    path: GaussianSelectedProjectionPath,
) -> SquaredEuclideanAverageLinkageRegion:
    """Construct and replay-certify the exact observed hierarchy interval union."""
    observed = _run_squared_euclidean_average_linkage(
        path,
        radius=path.observed_radius,
        collect_constraints=True,
    )
    boundary_radii = _unique_sorted_radii(
        [
            root
            for constraint in observed.selection_constraints
            for root in _nonnegative_real_roots(constraint)
        ]
    )
    intervals = _selected_intervals_from_constraints(
        observed.selection_constraints,
        boundary_radii,
    )
    if not intervals or not any(
        interval.lower <= path.observed_radius <= interval.upper
        for interval in intervals
    ):
        raise ValueError(
            "Observed hierarchy selection event has no positive-probability radial "
            "interval containing the observed radius."
        )
    for interval in intervals:
        replayed_signature = replay_squared_euclidean_average_linkage(
            path,
            radius=_interval_probe(interval),
        )
        if replayed_signature != observed.merge_signature:
            raise ValueError(
                "Analytic hierarchy interval failed merge-sequence replay certification."
            )

    return SquaredEuclideanAverageLinkageRegion(
        observed_merge_signature=observed.merge_signature,
        selection_intervals=intervals,
        constraint_count=len(observed.selection_constraints),
        nonnegative_boundary_radii=boundary_radii,
    )


__all__ = [
    "HierarchyMergeSignature",
    "SquaredEuclideanAverageLinkageRegion",
    "construct_squared_euclidean_average_linkage_region",
    "replay_squared_euclidean_average_linkage",
]
