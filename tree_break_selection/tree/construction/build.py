"""Directory-independent dispatch for supported tree-construction methods."""

from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

from tree_break_selection.tree.construction.defaults import SUPPORTED_LINKAGE_METHODS
from tree_break_selection.tree.construction.hierarchical import (
    tree_from_linkage,
)
from tree_break_selection.tree.poset_tree import PosetTree


@dataclass(frozen=True)
class TreeBuildResult:
    """Tree plus the linkage matrix it was built from."""

    tree: PosetTree
    linkage_matrix: np.ndarray


def _validate_condensed_distance(
    distance_condensed: np.ndarray,
    *,
    sample_count: int,
) -> np.ndarray:
    """Return a validated condensed distance vector."""
    distances = np.asarray(distance_condensed, dtype=float)
    expected_pairs = sample_count * (sample_count - 1) // 2
    if distances.ndim != 1 or len(distances) != expected_pairs:
        raise ValueError(
            "distance_condensed must contain exactly one value per unordered sample pair; "
            f"expected {expected_pairs}, got shape={distances.shape!r}."
        )
    if not np.isfinite(distances).all():
        raise ValueError("distance_condensed contains non-finite values.")
    if np.any(distances < 0.0):
        raise ValueError("distance_condensed contains negative values.")
    return distances


def _canonical_sample_order(labels: list[object]) -> list[int]:
    """Return a row-order-independent sample order for deterministic tie breaks."""
    if all(isinstance(label, str) for label in labels):

        def natural_string_key(position: int) -> tuple[tuple[tuple[int, object], ...], str]:
            label = labels[position]
            assert isinstance(label, str)
            chunks = tuple(
                (0, int(chunk)) if chunk.isdigit() else (1, chunk)
                for chunk in re.split(r"(\d+)", label)
            )
            return chunks, label

        return sorted(range(len(labels)), key=natural_string_key)

    try:
        return sorted(range(len(labels)), key=labels.__getitem__)
    except TypeError as exc:
        raise ValueError(
            "Sample labels must be mutually orderable, or all strings, for deterministic "
            "tree construction. Normalize mixed-type labels before building the tree."
        ) from exc


def _reorder_condensed_distance(
    distances: np.ndarray,
    order: list[int],
) -> np.ndarray:
    """Reorder a condensed distance vector without changing its values."""
    if order == list(range(len(order))):
        return distances
    square = squareform(distances)
    reordered = square[np.ix_(order, order)]
    return squareform(reordered, checks=False)


def build_tree(
    data: pd.DataFrame,
    distance_condensed: np.ndarray,
    *,
    linkage_method: str,
) -> TreeBuildResult:
    """Build a rooted linkage tree, rooted at the final merge, from condensed distances."""
    if len(data) < 2:
        raise ValueError("Tree construction requires at least two samples.")
    if not data.index.is_unique:
        raise ValueError("Tree construction requires unique sample labels.")
    if linkage_method not in SUPPORTED_LINKAGE_METHODS:
        raise ValueError(
            f"Unsupported linkage_method {linkage_method!r}; allowed={SUPPORTED_LINKAGE_METHODS!r}."
        )

    canonical_order = _canonical_sample_order(data.index.tolist())
    canonical_data = data.iloc[canonical_order]
    distances = _validate_condensed_distance(
        distance_condensed,
        sample_count=len(data),
    )
    distances = _reorder_condensed_distance(distances, canonical_order)
    linkage_matrix = linkage(distances, method=linkage_method)
    tree = tree_from_linkage(
        linkage_matrix,
        leaf_names=canonical_data.index.tolist(),
    )
    return TreeBuildResult(tree=tree, linkage_matrix=linkage_matrix)


__all__ = [
    "TreeBuildResult",
    "build_tree",
]
