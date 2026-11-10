"""Descriptive geometry of truth groups and explicitly supplied subtree splits."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from itertools import combinations

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class StatisticalProfiles:
    """Group-level, pairwise group, and explicit subtree profiles."""

    groups: pd.DataFrame
    pairs: pd.DataFrame
    subtrees: pd.DataFrame


@dataclass(frozen=True)
class _Sample:
    points: np.ndarray
    centroid: np.ndarray
    covariance: np.ndarray


def _sample(points: np.ndarray) -> _Sample:
    centroid = points.mean(axis=0)
    residual = points - centroid
    covariance = residual.T @ residual
    if not np.isfinite(covariance).all() or not np.isfinite(centroid).all():
        raise ValueError("Centroids and covariance must be finite.")
    return _Sample(points, centroid, covariance)


def _spectrum(eigenvalues: np.ndarray) -> tuple[tuple[float, ...], int, float]:
    values = np.maximum(eigenvalues[::-1], 0)
    if values[0] == 0:
        return tuple(values), 0, np.nan
    scaled = values / values[0]
    rank = int((scaled > 1e-10).sum())
    weights = scaled[scaled > 0] / scaled.sum()
    return tuple(values), rank, float(np.exp(-np.sum(weights * np.log(weights))))


def _group_profile(sample: _Sample, total_n: int) -> dict[str, object]:
    n = len(sample.points)
    covariance = sample.covariance / (n - 1) if n > 1 else None
    spectrum, rank, effective = (
        _spectrum(np.linalg.eigvalsh(covariance))
        if covariance is not None
        else ((np.nan,) * len(sample.centroid), np.nan, np.nan)
    )
    return {
        "size": n,
        "share": n / total_n,
        "effective_n": n,
        "covariance_eigenvalues": spectrum,
        "covariance_rank": rank,
        "effective_rank": effective,
    }


def _contrast_profile(left: _Sample, right: _Sample) -> dict[str, object]:
    a, b = len(left.points), len(right.points)
    delta = right.centroid - left.centroid
    distance = float(np.linalg.norm(delta))
    if not np.isfinite(distance):
        raise ValueError("Centroid distance must be finite.")
    mahalanobis, overlap, covariance_difference = np.nan, np.nan, np.nan
    if distance > 0:
        direction = delta / distance
        scores_left = np.einsum("ij,j->i", left.points - left.centroid, direction) - distance / 2
        scores_right = np.einsum("ij,j->i", right.points - right.centroid, direction) + distance / 2
        tolerance = 1e-12 * max(distance, np.abs(scores_left).max(), np.abs(scores_right).max())
        overlap = float(((scores_left > tolerance).mean() + (scores_right < -tolerance).mean()) / 2)
        overlap += float(
            ((np.abs(scores_left) <= tolerance).mean() + (np.abs(scores_right) <= tolerance).mean())
            / 4
        )
    if min(a, b) > 1:
        covariance_a = left.covariance / (a - 1)
        covariance_b = right.covariance / (b - 1)
        denominator = np.linalg.norm(covariance_a) + np.linalg.norm(covariance_b)
        if denominator > 0:
            covariance_difference = float(np.linalg.norm(covariance_a - covariance_b) / denominator)
    if a + b > 2:
        pooled = left.covariance / (a + b - 2) + right.covariance / (a + b - 2)
        eigenvalues, vectors = np.linalg.eigh(pooled)
        spectrum, rank, effective = _spectrum(eigenvalues)
        status = "rank_deficient"
        if rank == len(delta):
            mahalanobis = float(np.linalg.norm((delta @ vectors) / np.sqrt(eigenvalues)))
            status = "full_rank"
    else:
        spectrum, rank, effective = (np.nan,) * len(delta), np.nan, np.nan
        status = "no_covariance_degrees_of_freedom"
    return {
        "size": a + b,
        "left_size": a,
        "right_size": b,
        "smaller_share": min(a, b) / (a + b),
        "effective_n": a * b / (a + b),
        "centroid_distance": distance,
        "mahalanobis_distance": mahalanobis,
        "mahalanobis_status": status,
        "centroid_overlap_share": overlap,
        "pooled_covariance_eigenvalues": spectrum,
        "pooled_covariance_rank": rank,
        "pooled_effective_rank": effective,
        "relative_covariance_difference": covariance_difference,
    }


def _positions(values: Sequence[int], n: int) -> np.ndarray:
    rows = np.asarray(values)
    if (
        rows.ndim != 1
        or not rows.size
        or rows.dtype.kind not in "iu"
        or (rows < 0).any()
        or (rows >= n).any()
        or len(np.unique(rows)) != len(rows)
    ):
        raise ValueError("Subtree positions must be nonempty, unique, in-range integers.")
    return rows


def build_statistical_profiles(
    features: np.ndarray | pd.DataFrame,
    *,
    true_labels: np.ndarray | pd.Series | None = None,
    subtrees: Mapping[str, tuple[Sequence[int], Sequence[int]]] | None = None,
    geometry_disagreement: Mapping[str, float] | None = None,
) -> StatisticalProfiles:
    """Profile truth groups and explicit row-position subtree splits separately.

    Effective rank is the exponential entropy of the covariance spectrum.
    Contrast effective n is 1/(1/n_left+1/n_right). Overlap is the average
    in-sample opposite-centroid share across sides, with half weight for
    equidistant points (relative projection tolerance 1e-12). It is
    undefined at coincident centroids. Mahalanobis
    distance is unavailable if pooled covariance is numerically rank-deficient
    (eigenvalues at most 1e-10 times the largest are excluded from rank).
    """
    x = np.asarray(features, dtype=float)
    if x.ndim != 2 or min(x.shape) < 1 or not np.isfinite(x).all():
        raise ValueError("Expected a nonempty finite rows-by-features matrix.")
    if isinstance(features, pd.DataFrame) and not features.index.is_unique:
        raise ValueError("Feature index must be unique.")
    subtrees = {} if subtrees is None else subtrees
    geometry = {} if geometry_disagreement is None else geometry_disagreement
    if not geometry.keys() <= subtrees.keys() or not all(np.isfinite(v) for v in geometry.values()):
        raise ValueError("Supply finite geometry disagreement only for explicit subtrees.")

    group_rows, pair_rows, subtree_rows = [], [], []
    if true_labels is not None:
        if isinstance(true_labels, pd.Series) and isinstance(features, pd.DataFrame):
            if (
                not true_labels.index.is_unique
                or len(true_labels) != len(x)
                or not features.index.isin(true_labels.index).all()
            ):
                raise ValueError("Feature and truth indexes must match uniquely.")
            true_labels = true_labels.reindex(features.index)
        labels = np.asarray(true_labels)
        if labels.shape != (len(x),) or pd.isna(labels).any():
            raise ValueError("Truth labels must align to all feature rows without missing values.")
        codes, unique = pd.factorize(labels, sort=True)
        samples = {group: _sample(x[codes == index]) for index, group in enumerate(unique)}
        nearest = {}
        for left, right in combinations(samples, 2):
            row = _contrast_profile(samples[left], samples[right])
            pair_rows.append({"left_group": left, "right_group": right, **row})
            for own, other in ((left, right), (right, left)):
                if (
                    own not in nearest
                    or row["centroid_distance"] < nearest[own][1]["centroid_distance"]
                ):
                    nearest[own] = (other, row)
        for group, sample in samples.items():
            other, neighbor = nearest.get(group, (None, {}))
            group_rows.append(
                {
                    "group": group,
                    **_group_profile(sample, len(x)),
                    "nearest_group": other,
                    **{
                        f"nearest_{field}": neighbor.get(field, np.nan)
                        for field in (
                            "centroid_distance",
                            "mahalanobis_distance",
                            "centroid_overlap_share",
                            "relative_covariance_difference",
                        )
                    },
                }
            )

    for name, (left, right) in subtrees.items():
        a, b = _positions(left, len(x)), _positions(right, len(x))
        if np.intersect1d(a, b).size:
            raise ValueError("Explicit subtree sides must be disjoint.")
        row = _contrast_profile(_sample(x[a]), _sample(x[b]))
        subtree_rows.append(
            {"subtree": name, **row, "geometry_disagreement": geometry.get(name, np.nan)}
        )

    return StatisticalProfiles(
        pd.DataFrame(group_rows) if group_rows else pd.DataFrame(columns=["group"]),
        pd.DataFrame(pair_rows)
        if pair_rows
        else pd.DataFrame(columns=["left_group", "right_group"]),
        pd.DataFrame(subtree_rows) if subtree_rows else pd.DataFrame(columns=["subtree"]),
    )
