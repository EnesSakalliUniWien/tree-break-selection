"""Bounded geometry diagnostics on a fixed labeled subtree.

Leaf labels are row positions. Cut enumeration, Ward trees in spectral spans,
and the combinatorial topology family are separate objects.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from functools import cache
from itertools import combinations, product

import numpy as np
from scipy.cluster.hierarchy import cophenet, linkage
from scipy.optimize import nnls
from scipy.spatial.distance import pdist, squareform

from tree_break_selection.eigenvalue_groups import tied_eigenvalue_neighbors
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.split_count_rule_out import (
    is_separable,
)

MAX_CUT_LEAVES = 10
MAX_TOPOLOGY_LEAVES = 6

type Tree = int | tuple[Tree, Tree]


@dataclass(frozen=True)
class Cut:
    """Unordered split with the smallest leaf on the left."""

    left: tuple[int, ...]
    right: tuple[int, ...]


@dataclass(frozen=True, eq=False)
class SpectralSubspace:
    """Descending Gram eigenblock ranges and the resulting sample geometry."""

    eigenblocks: tuple[tuple[int, int], ...]
    rank: int
    projector: np.ndarray
    coordinates: np.ndarray


@dataclass(frozen=True, eq=False)
class SpectralGeometry:
    subspace: SpectralSubspace
    cuts: tuple[Cut, ...]
    distances: np.ndarray
    ward_linkage: np.ndarray
    ward_tree: Tree


@dataclass(frozen=True, eq=False)
class AdditiveDistanceFit:
    """Nonnegative lengths by unique edge split and fitted condensed distances."""

    edge_splits: tuple[Cut, ...]
    branch_lengths: np.ndarray
    fitted_distances: np.ndarray
    rmse: float
    max_abs_residual: float


def unordered_cuts(n: int) -> tuple[Cut, ...]:
    """All nonempty complementary splits, with row zero on the left."""
    if not 2 <= n <= MAX_CUT_LEAVES:
        raise ValueError(f"Cut enumeration requires 2 through {MAX_CUT_LEAVES} leaves.")
    return tuple(
        Cut(
            (0, *(j for j in range(1, n) if mask & (1 << (j - 1)))),
            tuple(j for j in range(1, n) if not mask & (1 << (j - 1))),
        )
        for mask in range((1 << (n - 1)) - 1)
    )


def spectral_subspaces(points: np.ndarray) -> tuple[SpectralSubspace, ...]:
    """Return each positive whole Gram block and each distinct leading block union.

    Each projector acts on row space. The coordinates retain Euclidean geometry
    within the projected feature space. SVD avoids explicitly squaring the
    condition number or allocating the feature Gram matrix.
    """
    x = np.asarray(points, dtype=np.float64)
    if x.ndim == 2 and not 2 <= len(x) <= MAX_CUT_LEAVES:
        raise ValueError(f"Subtree geometry requires 2 through {MAX_CUT_LEAVES} leaves.")
    return tuple(
        SpectralSubspace(group, basis.shape[1], basis @ basis.T, coordinates)
        for group, basis, coordinates in _spectral_coordinate_spans(x)
    )


def _spectral_coordinate_spans(
    points: np.ndarray,
) -> Iterator[tuple[tuple[tuple[int, int], ...], np.ndarray, np.ndarray]]:
    """Generate whole spectral spans without retaining dense row projectors."""
    x = np.asarray(points, dtype=np.float64)
    if x.ndim != 2 or min(x.shape) < 1 or not np.isfinite(x).all():
        raise ValueError("Expected a finite, nonempty rows-by-features matrix.")
    if len(x) < 2:
        raise ValueError("Subtree geometry requires at least two leaves.")

    # Scale before centering and SVD to keep the Gram spectrum finite.
    scale = float(np.max(np.abs(x)))
    if scale == 0:
        return
    normalized = x / scale
    normalized -= normalized[0].copy()
    u, singular, _ = np.linalg.svd(normalized - normalized.mean(axis=0), full_matrices=False)
    threshold = np.finfo(float).eps * max(x.shape) * singular[0]
    rank = min(len(x) - 1, int(np.count_nonzero(singular > threshold)))
    if rank == 0:
        return

    values = (singular[:rank] / singular[0]) ** 2
    ends = [*(np.flatnonzero(~tied_eigenvalue_neighbors(values)) + 1).tolist(), rank]
    blocks = tuple(zip([0, *ends[:-1]], ends, strict=True))
    groups = [(block,) for block in blocks] + [blocks[:k] for k in range(2, len(blocks) + 1)]
    for group in groups:
        start, stop = group[0][0], group[-1][1]
        row_basis = u[:, start:stop]
        coordinates = (row_basis * singular[start:stop]) * scale
        if not np.isfinite(coordinates).all():
            raise ValueError("Projected coordinates must fit in finite float64.")
        yield group, row_basis, coordinates


def spectral_geometry(points: np.ndarray) -> tuple[SpectralGeometry, ...]:
    """Enumerate strict affine cuts and compute a fixed-row-order Ward tree per span."""
    subspaces = spectral_subspaces(points)
    records = []
    for subspace in subspaces:
        x = subspace.coordinates
        # Column rescaling preserves affine feasibility and avoids scaling the LP
        # by the Gram eigenvalues or measurement units.
        lp_points = x / np.max(np.abs(x), axis=0)
        cuts = tuple(
            cut
            for cut in unordered_cuts(len(x))
            if is_separable(lp_points, np.isin(np.arange(len(x)), cut.left))
        )
        scale = np.max(np.abs(x))
        scaled_distances = pdist(x / scale)
        distances = scaled_distances * scale
        if not np.isfinite(distances).all():
            raise ValueError("Projected distances must fit in finite float64.")
        # shortcut: round relative input distances at 1e-12 for deterministic numerical
        # ties in fixed row order, replace if exact near-tie Ward decisions are required.
        max_distance = scaled_distances.max()
        ward = linkage(np.round(scaled_distances / max_distance, 12), method="ward")
        ward[:, 2] *= max_distance * scale
        if not np.isfinite(ward).all():
            raise ValueError("Ward linkage heights must fit in finite float64.")
        nodes: list[Tree] = list(range(len(x)))
        first_leaf = list(range(len(x)))
        for a, b in ward[:, :2].astype(int):
            if first_leaf[b] < first_leaf[a]:
                a, b = b, a
            nodes.append((nodes[a], nodes[b]))
            first_leaf.append(first_leaf[a])
        records.append(SpectralGeometry(subspace, cuts, distances, ward, nodes[-1]))
    return tuple(records)


def rooted_binary_topologies(n: int) -> tuple[Tree, ...]:
    """Enumerate all unordered rooted full binary trees of at most six labeled leaves."""
    if not 2 <= n <= MAX_TOPOLOGY_LEAVES:
        raise ValueError(f"Topology enumeration requires 2 through {MAX_TOPOLOGY_LEAVES} leaves.")

    def recurse(labels: tuple[int, ...]) -> tuple[Tree, ...]:
        if len(labels) == 1:
            return (labels[0],)
        trees = []
        for cut in unordered_cuts(len(labels)):
            left = recurse(tuple(labels[i] for i in cut.left))
            right = recurse(tuple(labels[i] for i in cut.right))
            trees.extend(product(left, right))
        return tuple(trees)

    return recurse(tuple(range(n)))


def ward_topology_family(
    points_or_condensed: np.ndarray, *, rtol: float = 1e-10
) -> frozenset[Tree]:
    """All Ward trees reachable by tied merges on two through six leaves.

    Accepts Euclidean coordinates or condensed Euclidean distances. At each
    step include every increase in within-cluster SSE within `rtol` times the
    largest original squared pair distance of the minimum. Zero uses exact
    floating-point comparisons.
    """
    x = np.asarray(points_or_condensed, dtype=float)
    if not np.isfinite(rtol) or rtol < 0:
        raise ValueError("rtol must be finite and nonnegative.")
    if not x.size or not np.isfinite(x).all() or x.ndim not in (1, 2):
        raise ValueError("Expected finite coordinates or condensed Euclidean distances.")
    if x.ndim == 1 and (x < 0).any():
        raise ValueError("Euclidean distances must be nonnegative.")
    scale = np.max(np.abs(x))
    if scale > 0:
        x = x / scale
    if x.ndim == 2 and not 2 <= len(x) <= MAX_TOPOLOGY_LEAVES:
        raise ValueError(
            f"Ward topology enumeration requires 2 through {MAX_TOPOLOGY_LEAVES} leaves."
        )
    distances = pdist(x) if x.ndim == 2 else x
    if not 1 <= len(distances) <= MAX_TOPOLOGY_LEAVES * (MAX_TOPOLOGY_LEAVES - 1) // 2:
        raise ValueError(
            f"Ward topology enumeration requires 2 through {MAX_TOPOLOGY_LEAVES} leaves."
        )
    max_distance = distances.max()
    normalized = distances / max_distance if max_distance > 0 else distances
    squared = squareform(normalized**2)
    n = len(squared)
    if x.ndim == 1:
        centered = np.eye(n) - 1 / n
        if np.linalg.eigvalsh(-0.5 * centered @ squared @ centered)[0] < -1e-12:
            raise ValueError("Condensed distances must be Euclidean.")

    @cache
    def merge_cost(left: tuple[int, ...], right: tuple[int, ...]) -> float:
        a, b = len(left), len(right)
        cross = squared[np.ix_(left, right)].sum()
        within_a = sum(squared[i, j] for i, j in combinations(left, 2))
        within_b = sum(squared[i, j] for i, j in combinations(right, 2))
        # (ab/(a+b)) times the squared centroid distance, via pair distances.
        return max(0.0, (a * b * cross - b**2 * within_a - a**2 * within_b) / (a * b * (a + b)))

    @cache
    def expand(forest: tuple[tuple[tuple[int, ...], Tree], ...]) -> frozenset[Tree]:
        if len(forest) == 1:
            return frozenset((forest[0][1],))
        costs = [
            (merge_cost(forest[i][0], forest[j][0]), i, j)
            for i, j in combinations(range(len(forest)), 2)
        ]
        threshold = min(cost for cost, _, _ in costs) + rtol
        trees = set()
        for cost, i, j in costs:
            if cost <= threshold:
                left, right = forest[i], forest[j]
                merged = (tuple(sorted(left[0] + right[0])), (left[1], right[1]))
                trees.update(expand(forest[:i] + (merged,) + forest[i + 1 : j] + forest[j + 1 :]))
        return frozenset(trees)

    return expand(tuple(((i,), i) for i in range(n)))


def _clades(tree: Tree) -> list[frozenset[int]]:
    if isinstance(tree, int):
        return [frozenset((tree,))]
    left, right = map(_clades, tree)
    return [*left, *right, left[-1] | right[-1]]


def tree_clades(tree: Tree) -> frozenset[frozenset[int]]:
    """Nontrivial clades for rooted topology comparison."""
    return frozenset(clade for clade in _clades(tree)[:-1] if len(clade) > 1)


def tree_edge_splits(tree: Tree) -> frozenset[Cut]:
    """Unordered edge splits including pendant edges, deduplicating the root split."""
    clades = _clades(tree)
    leaves = clades.pop()
    first = min(leaves)
    splits = set()
    for clade in clades:
        left, right = (clade, leaves - clade) if first in clade else (leaves - clade, clade)
        splits.add(Cut(tuple(sorted(left)), tuple(sorted(right))))
    return frozenset(splits)


def fit_additive_distances(tree: Tree, condensed_distances: np.ndarray) -> AdditiveDistanceFit:
    """Fit all leaf-pair distances to unique edge splits on a small fixed tree.

    Targets use the supplied distance units. Complementary root edges have
    the same incidence column, so their lengths are represented by one sum.
    Zero fitted internal lengths give unresolved splits in the fitted metric.
    """
    clades = _clades(tree)
    n = len(clades[-1])
    if not 2 <= n <= MAX_TOPOLOGY_LEAVES:
        raise ValueError(
            f"Additive distance fitting requires 2 through {MAX_TOPOLOGY_LEAVES} leaves."
        )
    if clades[-1] != frozenset(range(n)) or sum(len(clade) == 1 for clade in clades) != n:
        raise ValueError("Expected unique leaf labels 0 through n-1.")
    observed = np.asarray(condensed_distances, dtype=float)
    if (
        observed.shape != (n * (n - 1) // 2,)
        or not np.isfinite(observed).all()
        or (observed < 0).any()
    ):
        raise ValueError("Expected finite nonnegative condensed distances for the supplied leaves.")

    splits = tuple(sorted(tree_edge_splits(tree), key=lambda split: split.left))
    membership = np.column_stack([np.isin(np.arange(n), split.left) for split in splits])
    left, right = np.triu_indices(n, k=1)
    design = (membership[left] != membership[right]).astype(float)
    scale = float(observed.max()) or 1.0
    normalized = observed / scale
    lengths, _ = nnls(design, normalized)
    fitted = design @ lengths
    residuals = fitted - normalized
    lengths, fitted = lengths * scale, fitted * scale
    if not np.isfinite(fitted).all():
        raise ValueError("Fitted distances must fit in finite float64.")
    return AdditiveDistanceFit(
        splits,
        lengths,
        fitted,
        float(np.sqrt(np.mean(residuals**2))) * scale,
        float(np.max(np.abs(residuals))) * scale,
    )


def rooted_robinson_foulds(left: Tree, right: Tree) -> int:
    """Number of nontrivial rooted clades present in only one of the trees."""
    if _clades(left)[-1] != _clades(right)[-1]:
        raise ValueError("Rooted RF requires the same leaf labels in both trees.")
    return len(tree_clades(left) ^ tree_clades(right))


def cophenetic_distances(tree_or_linkage: Tree | np.ndarray) -> np.ndarray:
    """Condensed cophenetic distances in row order.

    An abstract topology assigns a node the height of its longest path to a
    leaf in edge units. A linkage matrix preserves its supplied heights.
    """
    if isinstance(tree_or_linkage, np.ndarray):
        return cophenet(tree_or_linkage)

    return cophenet(topology_linkage(tree_or_linkage))


def topology_linkage(tree: Tree) -> np.ndarray:
    """Encode a topology using monotone heights in edge units."""
    leaves = _clades(tree)[-1]
    n = len(leaves)
    if n < 2 or leaves != frozenset(range(n)):
        raise ValueError("Cophenetic topology leaves must be numbered 0 through n-1.")
    merges = []

    def build(node: Tree) -> tuple[int, int, int]:
        if isinstance(node, int):
            return node, 0, 1
        a, ha, na = build(node[0])
        b, hb, nb = build(node[1])
        parent = n + len(merges)
        height, count = 1 + max(ha, hb), na + nb
        merges.append((a, b, height, count))
        return parent, height, count

    _, _, count = build(tree)
    if count != n:
        raise ValueError("Cophenetic topology leaves must be unique.")
    return np.asarray(merges, dtype=float)
