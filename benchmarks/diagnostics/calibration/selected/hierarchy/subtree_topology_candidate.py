"""Bounded, all-pair additive distance ranking of alternatives within one subtree."""

from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import is_valid_linkage, linkage
from scipy.linalg import cholesky, solve_triangular
from scipy.optimize import nnls
from scipy.spatial.distance import pdist, squareform

from tree_break_selection.tree.optimized_branch_lengths import (
    BRANCH_LENGTH_TARGET_MEAN_SQUARED_STANDARDIZED_EUCLIDEAN,
    _mean_squared_standardized_targets,
)

from .subtree_topology_geometry import (
    MAX_TOPOLOGY_LEAVES,
    AdditiveDistanceFit,
    Cut,
    Tree,
    _spectral_coordinate_spans,
    rooted_binary_topologies,
    topology_linkage,
    tree_edge_splits,
)

MAX_CANDIDATE_LEAVES = 256
FIT_TIE_RTOL = 1e-10


@dataclass(frozen=True, eq=False)
class CandidateFit:
    topology: Tree
    sources: tuple[str, ...]
    fit: AdditiveDistanceFit
    zero_internal_edge_count: int
    equivalent_to_selected_metric: bool


@dataclass(frozen=True, eq=False)
class SubtreeTopologySelection:
    selected_topology: Tree
    selected_linkage: np.ndarray
    current_topology: Tree
    candidates: tuple[CandidateFit, ...]
    candidate_count: int
    spectral_span_count: int
    coverage: str
    pair_count: int
    target_metric: str
    elapsed_sec: float


def topology_from_linkage(matrix: np.ndarray) -> Tree:
    """Read the unranked topology from a SciPy linkage over local row positions."""
    matrix = np.asarray(matrix, dtype=float)
    if matrix.ndim != 2 or matrix.shape[1] != 4 or not np.isfinite(matrix).all():
        raise ValueError("Linkage must be a finite matrix with four columns")
    integral = matrix[:, [0, 1, 3]]
    if np.any(integral != np.floor(integral)):
        raise ValueError("Linkage child IDs and leaf counts must be integers")
    is_valid_linkage(matrix, throw=True)
    nodes: list[Tree] = list(range(len(matrix) + 1))
    smallest = list(range(len(nodes)))
    counts = [1] * len(nodes)
    heights = [0.0] * len(nodes)
    for (a, b), (_, _, height, count) in zip(matrix[:, :2].astype(int), matrix, strict=True):
        if count != counts[a] + counts[b]:
            raise ValueError("Linkage leaf count must equal the sum of its child counts")
        if height < max(heights[a], heights[b]):
            raise ValueError("Linkage merge height cannot be below a child height")
        if smallest[b] < smallest[a]:
            a, b = b, a
        nodes.append((nodes[a], nodes[b]))
        smallest.append(smallest[a])
        counts.append(int(count))
        heights.append(height)
    return nodes[-1]


def _relabel(tree: Tree, labels: np.ndarray) -> Tree:
    if isinstance(tree, int):
        return int(labels[tree])
    return (_relabel(tree[0], labels), _relabel(tree[1], labels))


def _local_ward(coordinates: np.ndarray) -> np.ndarray:
    scale = float(np.max(np.abs(coordinates)))
    distances = pdist(coordinates / scale)
    maximum = float(distances.max())
    z = linkage(np.round(distances / maximum, decimals=12), method="ward")
    z[:, 2] *= maximum * scale
    if not np.isfinite(z).all():
        raise ValueError("Ward heights must fit in finite float64")
    return z


def _relabel_linkage(matrix: np.ndarray, order: np.ndarray) -> np.ndarray:
    z = matrix.copy()
    children = z[:, :2].astype(int)
    is_leaf = children < len(order)
    z[:, :2][is_leaf] = order[children[is_leaf]]
    return z


def _fit_from_split_gram(
    splits: frozenset[Cut], normalized_target: np.ndarray, scale: float
) -> AdditiveDistanceFit:
    """Fit the exact all-pair NNLS objective without a pair-by-edge matrix."""
    ordered = tuple(sorted(splits, key=lambda split: split.left))
    n = len(normalized_target)
    member = np.column_stack([np.isin(np.arange(n), split.left) for split in ordered]).astype(float)
    inside = member.sum(axis=0)
    both = member.T @ member
    neither = n - inside[:, None] - inside[None, :] + both
    only_left = inside[:, None] - both
    gram = both * neither + only_left * only_left.T
    rhs = np.sum(member * (normalized_target @ (1 - member)), axis=0)
    # Unique splits eliminate the duplicate root columns. Distinct splits of a
    # binary tree have linearly independent pair-incidence columns.
    factor = cholesky(gram, lower=False)
    transformed = solve_triangular(factor, rhs, trans="T")
    lengths, _ = nnls(factor, transformed)
    weighted_member = member * lengths
    totals = weighted_member.sum(axis=1)
    fitted = totals[:, None] + totals[None, :] - 2 * (weighted_member @ member.T)
    left, right = np.triu_indices(n, k=1)
    fitted = np.maximum(0.0, fitted[left, right])
    residual = fitted - normalized_target[left, right]
    with np.errstate(over="raise", invalid="raise"):
        return AdditiveDistanceFit(
            ordered,
            lengths * scale,
            fitted * scale,
            float(np.sqrt(np.mean(residual**2)) * scale),
            float(np.max(np.abs(residual)) * scale),
        )


def select_subtree_topology(
    points: np.ndarray,
    *,
    current_topology: Tree,
    target_metric: str = "euclidean",
    leaf_labels=None,
    exhaustive: bool = False,
    current_linkage: np.ndarray | None = None,
) -> SubtreeTopologySelection:
    """Rank current and spectral Ward topologies on one shared full-subtree target.

    Search is bounded to 256 leaves and uses every pair. `exhaustive=True` adds
    all rooted binary trees only through six leaves. Supply unique sortable
    leaf labels to keep Ward ties and residual tie-breaking stable on row
    reorder. Returned topologies and linkages use the caller's local row order.
    Original linkage heights are retained for the current tree if supplied.
    Exhaustive-only candidates have topology heights, not fitted branch times.
    """
    start = perf_counter()
    x = np.asarray(points, dtype=float)
    if (
        x.ndim != 2
        or not 2 <= len(x) <= MAX_CANDIDATE_LEAVES
        or x.shape[1] < 1
        or not np.isfinite(x).all()
    ):
        raise ValueError(f"Expected finite coordinates for 2 through {MAX_CANDIDATE_LEAVES} leaves")
    n = len(x)
    if exhaustive and n > MAX_TOPOLOGY_LEAVES:
        raise ValueError(
            f"Exhaustive topology enumeration requires at most {MAX_TOPOLOGY_LEAVES} leaves"
        )
    pending = [current_topology]
    seen = set()
    while pending:
        node = pending.pop()
        if isinstance(node, tuple) and len(node) == 2:
            pending.extend(node)
        elif type(node) is int and 0 <= node < n and node not in seen:
            seen.add(node)
        else:
            raise ValueError("Current topology must contain each local leaf exactly once")
    if len(seen) != n:
        raise ValueError("Current topology must contain each local leaf exactly once")
    topology_matrix = topology_linkage(current_topology)
    current = topology_from_linkage(topology_matrix)
    z_current = (
        topology_matrix if current_linkage is None else np.asarray(current_linkage, dtype=float)
    )
    if topology_from_linkage(z_current) != current:
        raise ValueError("Current linkage must encode current_topology")

    labels = tuple(range(n)) if leaf_labels is None else tuple(leaf_labels)
    if len(labels) != n or len(set(labels)) != n:
        raise ValueError("leaf_labels must contain one unique sortable label per row")
    order = np.asarray(sorted(range(n), key=labels.__getitem__), dtype=np.intp)
    inverse = np.argsort(order)
    canonical_current = topology_from_linkage(topology_linkage(_relabel(current, inverse)))

    if target_metric == "euclidean":
        target = pdist(x)
    elif target_metric == BRANCH_LENGTH_TARGET_MEAN_SQUARED_STANDARDIZED_EUCLIDEAN:
        left, right = np.triu_indices(n, k=1)
        target = _mean_squared_standardized_targets(pd.DataFrame(x), left, right)
    else:
        raise ValueError(f"Unsupported subtree target metric: {target_metric!r}")
    if not np.isfinite(target).all():
        raise ValueError("Full-subtree distance target must be finite")
    scale = float(target.max()) or 1.0
    normalized_target = squareform(target / scale)

    sources = {canonical_current: ["current"]}
    linkages = {canonical_current: z_current}
    span_count = 0
    for group, _, coordinates in _spectral_coordinate_spans(x[order]):
        z = _local_ward(coordinates)
        tree = topology_from_linkage(z)
        sources.setdefault(tree, []).append(f"spectral:{group}")
        if tree not in linkages:
            linkages[tree] = _relabel_linkage(z, order)
        span_count += 1
    if exhaustive:
        for tree in rooted_binary_topologies(n):
            sources.setdefault(tree, []).append("exhaustive")

    trees = {key: topology_from_linkage(topology_linkage(_relabel(key, order))) for key in sources}
    by_splits = {}
    fits = {}
    for key, tree in trees.items():
        splits = tree_edge_splits(tree)
        if splits not in by_splits:
            by_splits[splits] = _fit_from_split_gram(splits, normalized_target, scale)
        fits[key] = by_splits[splits]

    tolerance = FIT_TIE_RTOL * float(target.max())
    best = min(fit.rmse for fit in fits.values())
    tied = [key for key, fit in fits.items() if fit.rmse <= best + tolerance]
    selected = canonical_current if canonical_current in tied else min(tied, key=repr)
    selected_fit = fits[selected]
    candidates = []
    for key in sorted(fits, key=repr):
        fit = fits[key]
        zero_internal = sum(
            len(split.left) > 1 and len(split.right) > 1 and length <= tolerance
            for split, length in zip(fit.edge_splits, fit.branch_lengths, strict=True)
        )
        same_metric = bool(
            np.max(np.abs(fit.fitted_distances - selected_fit.fitted_distances)) <= tolerance
        )
        candidates.append(
            CandidateFit(trees[key], tuple(sources[key]), fit, zero_internal, same_metric)
        )

    selected_linkage = linkages.get(selected)
    if selected_linkage is None:
        selected_linkage = topology_linkage(trees[selected])
    return SubtreeTopologySelection(
        selected_topology=trees[selected],
        selected_linkage=selected_linkage,
        current_topology=current,
        candidates=tuple(candidates),
        candidate_count=len(candidates),
        spectral_span_count=span_count,
        coverage="exhaustive_rooted" if exhaustive else "bounded_spectral_ward",
        pair_count=len(target),
        target_metric=target_metric,
        elapsed_sec=perf_counter() - start,
    )
