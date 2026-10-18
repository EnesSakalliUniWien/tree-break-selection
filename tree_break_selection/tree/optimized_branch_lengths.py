"""Geometry-based branch-length optimization for rooted ``PosetTree`` objects.

The default linkage constructor stores ultrametric merge-height differences as
``branch_length``.  Those heights are useful topology diagnostics, but they are
not calibrated stochastic time.  This module provides a native fixed-topology
non-negative least-squares fit for continuous data:

    D_ij ~= sum_{e in path(i, j)} tau_e,  tau_e >= 0

where ``D_ij`` is either the mean squared standardized Euclidean distance
across all input features or the unnormalized squared Euclidean distance
from an already prepared geometry. The target identifier
``mean_squared_standardized_euclidean`` explicitly includes that normalization.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from time import perf_counter

import networkx as nx
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import lsq_linear

from tree_break_selection.tree.topology import tree_root

BRANCH_LENGTH_OPTIMIZATION_LINKAGE_ULTRAMETRIC = "linkage_ultrametric"
BRANCH_LENGTH_OPTIMIZATION_FIXED_TOPOLOGY_NNLS = "fixed_topology_nnls"
BRANCH_LENGTH_OPTIMIZATION_METHODS = (
    BRANCH_LENGTH_OPTIMIZATION_LINKAGE_ULTRAMETRIC,
    BRANCH_LENGTH_OPTIMIZATION_FIXED_TOPOLOGY_NNLS,
)

BRANCH_LENGTH_TARGET_MEAN_SQUARED_STANDARDIZED_EUCLIDEAN = "mean_squared_standardized_euclidean"
BRANCH_LENGTH_TARGET_SQUARED_EUCLIDEAN = "squared_euclidean"
BRANCH_LENGTH_TARGET_METRICS = (
    BRANCH_LENGTH_TARGET_MEAN_SQUARED_STANDARDIZED_EUCLIDEAN,
    BRANCH_LENGTH_TARGET_SQUARED_EUCLIDEAN,
)


@dataclass(frozen=True)
class BranchLengthOptimizationResult:
    """Summary of one branch-length optimization pass."""

    method: str
    target_metric: str
    status: str
    n_leaves: int
    n_edges: int
    n_pairs_total: int
    n_pairs_used: int
    design_nnz: int
    design_density: float
    n_zero_design_columns: int
    zero_design_column_fraction: float
    pair_sample_size: int | None
    random_state: int
    applied_to_tree: bool
    apply_nonconverged: bool
    elapsed_sec: float
    cost: float
    optimality: float
    iterations: int
    target_mean: float
    target_median: float
    fitted_mean: float
    fitted_median: float
    residual_rmse: float
    residual_mae: float
    residual_rmse_to_target_mean: float
    residual_mae_to_target_mean: float
    branch_length_mean: float
    branch_length_median: float
    branch_length_min: float
    branch_length_max: float
    solver_message: str
    target_definition: str
    target_normalization: str
    target_normalization_divisor: int

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def validate_branch_length_optimization_method(value: str) -> str:
    """Return a valid branch-length optimization method id."""
    method = str(value)
    if method not in BRANCH_LENGTH_OPTIMIZATION_METHODS:
        raise ValueError(
            "branch_length_optimization_method must be one of "
            f"{BRANCH_LENGTH_OPTIMIZATION_METHODS!r}; got {value!r}."
        )
    return method


def _leaf_nodes_for_data(tree: nx.DiGraph, data_index: pd.Index) -> list[object]:
    label_to_node: dict[str, object] = {}
    for node_id, attrs in tree.nodes(data=True):
        if attrs["is_leaf"]:
            label = str(attrs["label"])
            if label in label_to_node:
                raise ValueError(f"Duplicate tree leaf label {label!r}.")
            label_to_node[label] = node_id

    leaf_nodes: list[object] = []
    missing_labels: list[str] = []
    for label in data_index.astype(str):
        node_id = label_to_node.get(str(label))
        if node_id is None:
            missing_labels.append(str(label))
        else:
            leaf_nodes.append(node_id)
    if missing_labels:
        preview = ", ".join(missing_labels[:5])
        raise ValueError(f"Tree is missing {len(missing_labels)} data leaf label(s): {preview}.")
    return leaf_nodes


_PAIR_CHUNK_SIZE = 4096


def _edge_indices(tree: nx.DiGraph) -> dict[tuple[object, object], int]:
    return {(parent, child): index for index, (parent, child) in enumerate(tree.edges())}


def _leaf_edge_incidence(
    tree: nx.DiGraph,
    leaf_nodes: list[object],
    edge_index_by_pair: dict[tuple[object, object], int],
) -> sparse.csr_matrix:
    """Return a leaf-by-edge 0/1 matrix marking the edges on each leaf's root path."""
    root = tree_root(tree)
    indptr = [0]
    edge_ids: list[int] = []
    for leaf in leaf_nodes:
        current = leaf
        while current != root:
            parents = list(tree.predecessors(current))
            if len(parents) != 1:
                raise ValueError(
                    f"Node {current!r} must have exactly one parent; found {len(parents)}."
                )
            parent = parents[0]
            edge_ids.append(edge_index_by_pair[(parent, current)])
            current = parent
        indptr.append(len(edge_ids))
    return sparse.csr_matrix(
        (np.ones(len(edge_ids), dtype=np.float64), edge_ids, indptr),
        shape=(len(leaf_nodes), len(edge_index_by_pair)),
    )


def _sample_leaf_pairs(
    n_leaves: int,
    *,
    pair_sample_size: int | None,
    random_state: int,
) -> tuple[np.ndarray, np.ndarray, int]:
    total_pairs = n_leaves * (n_leaves - 1) // 2
    if pair_sample_size is None or int(pair_sample_size) >= total_pairs:
        left, right = np.triu_indices(n_leaves, k=1)
        return left.astype(np.int32), right.astype(np.int32), total_pairs

    target_size = int(pair_sample_size)
    if target_size <= 0:
        raise ValueError(f"pair_sample_size must be positive; got {pair_sample_size!r}.")
    rng = np.random.default_rng(int(random_state))
    pairs: set[tuple[int, int]] = set()
    batch_size = max(4096, target_size)
    while len(pairs) < target_size:
        left = rng.integers(0, n_leaves, size=batch_size)
        right = rng.integers(0, n_leaves, size=batch_size)
        different = left != right
        left = left[different]
        right = right[different]
        low = np.minimum(left, right)
        high = np.maximum(left, right)
        for pair in zip(low.tolist(), high.tolist(), strict=True):
            pairs.add((int(pair[0]), int(pair[1])))
            if len(pairs) >= target_size:
                break

    ordered_pairs = np.array(sorted(pairs), dtype=np.int32)
    return ordered_pairs[:, 0], ordered_pairs[:, 1], total_pairs


def _edge_crossing_pairs(
    leaf_edge_incidence: sparse.csr_matrix,
    edges: np.ndarray,
    edge_by_index: dict[int, tuple[object, object]],
) -> tuple[np.ndarray, np.ndarray]:
    """Return one leaf pair across each edge: its first leaf below and first leaf outside."""
    below = leaf_edge_incidence.tocsc()
    all_leaves = np.arange(below.shape[0])
    pairs: set[tuple[int, int]] = set()
    for edge in edges.tolist():
        inside = below.indices[below.indptr[edge] : below.indptr[edge + 1]]
        outside = np.setdiff1d(all_leaves, inside)
        if inside.size == 0 or outside.size == 0:
            raise ValueError(
                f"Edge {edge_by_index[edge]!r} has data leaves on one side only; "
                "no leaf pair crosses it, so its length cannot be fit."
            )
        pairs.add((int(inside.min()), int(outside[0])))
    ordered = np.array(sorted(pairs), dtype=np.int32)
    return ordered[:, 0], ordered[:, 1]


def _split_root_pair_by_linkage(
    tree: nx.DiGraph,
    lengths: np.ndarray,
    edge_index_by_pair: dict[tuple[object, object], int],
) -> np.ndarray:
    """Split the fitted root-pair sum in the ratio of the root edges' linkage lengths.

    Every leaf pair crosses both root edges or neither, so NNLS identifies only
    their sum. The pre-fit linkage lengths carry the tree builder's root
    placement; missing, non-finite or all-zero linkage lengths split evenly.
    """
    root = tree_root(tree)
    children = list(tree.successors(root))
    if len(children) != 2:
        return lengths
    indices = [edge_index_by_pair[(root, child)] for child in children]
    weights = np.array(
        [
            float(
                tree.edges[root, child].get(
                    "linkage_branch_length", tree.edges[root, child].get("branch_length", 0.0)
                )
            )
            for child in children
        ]
    )
    if not np.all(np.isfinite(weights)) or np.any(weights < 0.0) or weights.sum() <= 0.0:
        weights = np.ones(2)
    split = lengths.copy()
    split[indices] = lengths[indices].sum() * weights / weights.sum()
    return split


def _mean_squared_standardized_targets(
    data: pd.DataFrame,
    left: np.ndarray,
    right: np.ndarray,
) -> np.ndarray:
    """Return mean squared standardized Euclidean distances for sample pairs.

    Each feature is centered and scaled by its sample standard deviation
    (``ddof=1``); non-finite or non-positive scales are replaced by one. Sum
    the squared standardized differences and divide by the total number of
    input features, including constant features.
    """
    matrix = data.to_numpy(dtype=np.float64, copy=True)
    if matrix.ndim != 2 or matrix.shape[0] < 2 or matrix.shape[1] < 1:
        raise ValueError("Branch-length optimization requires a 2-D continuous data matrix.")
    center = np.mean(matrix, axis=0, keepdims=True)
    scale = np.std(matrix, axis=0, ddof=1, keepdims=True)
    scale[~np.isfinite(scale) | (scale <= 0.0)] = 1.0
    standardized = (matrix - center) / scale
    diffs = standardized[left] - standardized[right]
    return np.sum(diffs * diffs, axis=1) / float(standardized.shape[1])


def _squared_euclidean_targets(
    data: pd.DataFrame,
    left: np.ndarray,
    right: np.ndarray,
) -> np.ndarray:
    matrix = data.to_numpy(dtype=np.float64, copy=False)
    if matrix.ndim != 2 or matrix.shape[0] < 2 or matrix.shape[1] < 1:
        raise ValueError("Branch-length optimization requires a 2-D geometry matrix.")
    diffs = matrix[left] - matrix[right]
    return np.sum(diffs * diffs, axis=1)


def _path_incidence_matrix(
    leaf_edge_incidence: sparse.csr_matrix,
    left: np.ndarray,
    right: np.ndarray,
) -> sparse.csr_matrix:
    """Mark the edges on each leaf pair's path: the symmetric difference of root paths."""
    # Chunked so the shared root-path prefixes of a chunk are the only transient copies.
    blocks = []
    for start in range(0, len(left), _PAIR_CHUNK_SIZE):
        stop = start + _PAIR_CHUNK_SIZE
        block = abs(leaf_edge_incidence[left[start:stop]] - leaf_edge_incidence[right[start:stop]])
        block.eliminate_zeros()
        blocks.append(block)
    return sparse.vstack(blocks, format="csr")


def fit_fixed_topology_nnls_branch_lengths(
    tree: nx.DiGraph,
    continuous_data: pd.DataFrame,
    *,
    target_metric: str = BRANCH_LENGTH_TARGET_MEAN_SQUARED_STANDARDIZED_EUCLIDEAN,
    pair_sample_size: int | None = 100_000,
    random_state: int = 0,
    solver_tolerance: float = 1e-6,
    max_iterations: int | None = None,
    apply_nonconverged: bool = False,
) -> BranchLengthOptimizationResult:
    """Fit non-negative edge lengths on a fixed tree topology.

    The function mutates ``tree`` by replacing each edge's ``branch_length``
    with the fitted value.  The previous value is preserved as
    ``linkage_branch_length`` when present. Only the sum of the root's two edge
    lengths is identified; it is split in the ratio of their linkage lengths.
    An edge no sampled pair crosses gets one crossing pair, so every edge
    length is fit; ``n_zero_design_columns`` counts those edges. ``continuous_data`` is interpreted
    as raw continuous observations for the standardized target and as a
    prepared Euclidean embedding for the unstandardized target. The
    ``mean_squared_standardized_euclidean`` target divides the sum of squared
    sample-standardized differences by the total input feature count;
    ``squared_euclidean`` uses the unnormalized sum in the supplied geometry.
    """
    start_sec = perf_counter()
    if target_metric not in BRANCH_LENGTH_TARGET_METRICS:
        raise ValueError(
            "Unsupported branch-length target metric "
            f"{target_metric!r}; expected one of {BRANCH_LENGTH_TARGET_METRICS!r}."
        )
    if not isinstance(continuous_data, pd.DataFrame):
        raise TypeError("continuous_data must be a pandas DataFrame.")
    if continuous_data.isna().to_numpy().any():
        raise ValueError("continuous_data must not contain missing values.")

    leaf_nodes = _leaf_nodes_for_data(tree, continuous_data.index)
    edge_index_by_pair = _edge_indices(tree)
    n_edges = len(edge_index_by_pair)
    if n_edges == 0:
        raise ValueError("Branch-length optimization requires a non-empty tree.")

    leaf_edge_incidence = _leaf_edge_incidence(tree, leaf_nodes, edge_index_by_pair)
    edge_by_index = {index: edge for edge, index in edge_index_by_pair.items()}
    left, right, total_pairs = _sample_leaf_pairs(
        len(leaf_nodes),
        pair_sample_size=pair_sample_size,
        random_state=random_state,
    )
    design = _path_incidence_matrix(leaf_edge_incidence, left, right)
    # An edge no sampled pair crosses has a zero design column; the solver would
    # leave it at its starting value. Add one crossing pair for each such edge.
    uncovered_edges = np.flatnonzero(design.getnnz(axis=0) == 0)
    n_zero_design_columns = int(uncovered_edges.size)
    if n_zero_design_columns:
        extra_left, extra_right = _edge_crossing_pairs(
            leaf_edge_incidence, uncovered_edges, edge_by_index
        )
        left = np.concatenate([left, extra_left])
        right = np.concatenate([right, extra_right])
        design = sparse.vstack(
            [design, _path_incidence_matrix(leaf_edge_incidence, extra_left, extra_right)],
            format="csr",
        )
    if target_metric == BRANCH_LENGTH_TARGET_MEAN_SQUARED_STANDARDIZED_EUCLIDEAN:
        targets = _mean_squared_standardized_targets(continuous_data, left, right)
        target_definition = (
            "mean squared standardized Euclidean distance "
            "(sample standard deviation with ddof=1; "
            "non-finite or non-positive scales replaced by 1)"
        )
        target_normalization = "divide_by_total_feature_count"
        target_normalization_divisor = int(continuous_data.shape[1])
    else:
        targets = _squared_euclidean_targets(continuous_data, left, right)
        target_definition = "squared Euclidean distance in the supplied geometry"
        target_normalization = "none"
        target_normalization_divisor = 1
    design_nnz = int(design.nnz)
    design_density = float(design_nnz / (design.shape[0] * design.shape[1]))

    solution = lsq_linear(
        design,
        targets,
        bounds=(0.0, np.inf),
        method="trf",
        tol=float(solver_tolerance),
        max_iter=max_iterations,
        lsmr_tol="auto",
        verbose=0,
    )
    lengths = _split_root_pair_by_linkage(
        tree, np.asarray(solution.x, dtype=float), edge_index_by_pair
    )
    fitted = np.asarray(design @ lengths, dtype=float)
    residuals = fitted - targets
    solver_converged = bool(solution.success)
    applied_to_tree = solver_converged or bool(apply_nonconverged)

    if applied_to_tree:
        for index, length in enumerate(lengths):
            parent, child = edge_by_index[index]
            attrs = tree.edges[parent, child]
            if "branch_length" in attrs and "linkage_branch_length" not in attrs:
                attrs["linkage_branch_length"] = attrs["branch_length"]
            attrs["branch_length"] = float(length)
            attrs["branch_length_source"] = BRANCH_LENGTH_OPTIMIZATION_FIXED_TOPOLOGY_NNLS

    tree.graph["branch_length_optimization"] = {
        "method": BRANCH_LENGTH_OPTIMIZATION_FIXED_TOPOLOGY_NNLS,
        "target_metric": target_metric,
        "target_definition": target_definition,
        "target_normalization": target_normalization,
        "target_normalization_divisor": target_normalization_divisor,
        "status": "ok" if solver_converged else "solver_not_converged",
        "pair_sample_size": None if pair_sample_size is None else int(pair_sample_size),
        "random_state": int(random_state),
        "applied_to_tree": bool(applied_to_tree),
        "apply_nonconverged": bool(apply_nonconverged),
    }
    target_mean = float(np.mean(targets))
    residual_rmse = float(np.sqrt(np.mean(residuals * residuals)))
    residual_mae = float(np.mean(np.abs(residuals)))
    residual_rmse_to_target_mean = (
        float(residual_rmse / target_mean)
        if np.isfinite(target_mean) and target_mean > 0.0
        else float("nan")
    )
    residual_mae_to_target_mean = (
        float(residual_mae / target_mean)
        if np.isfinite(target_mean) and target_mean > 0.0
        else float("nan")
    )

    return BranchLengthOptimizationResult(
        method=BRANCH_LENGTH_OPTIMIZATION_FIXED_TOPOLOGY_NNLS,
        target_metric=target_metric,
        status="ok" if solver_converged else "solver_not_converged",
        n_leaves=len(leaf_nodes),
        n_edges=n_edges,
        n_pairs_total=total_pairs,
        n_pairs_used=len(left),
        design_nnz=design_nnz,
        design_density=design_density,
        n_zero_design_columns=n_zero_design_columns,
        zero_design_column_fraction=float(n_zero_design_columns / n_edges),
        pair_sample_size=None if pair_sample_size is None else int(pair_sample_size),
        random_state=int(random_state),
        applied_to_tree=bool(applied_to_tree),
        apply_nonconverged=bool(apply_nonconverged),
        elapsed_sec=perf_counter() - start_sec,
        cost=float(solution.cost),
        optimality=float(solution.optimality),
        iterations=int(solution.nit),
        target_mean=target_mean,
        target_median=float(np.median(targets)),
        fitted_mean=float(np.mean(fitted)),
        fitted_median=float(np.median(fitted)),
        residual_rmse=residual_rmse,
        residual_mae=residual_mae,
        residual_rmse_to_target_mean=residual_rmse_to_target_mean,
        residual_mae_to_target_mean=residual_mae_to_target_mean,
        branch_length_mean=float(np.mean(lengths)),
        branch_length_median=float(np.median(lengths)),
        branch_length_min=float(np.min(lengths)),
        branch_length_max=float(np.max(lengths)),
        solver_message=str(solution.message),
        target_definition=target_definition,
        target_normalization=target_normalization,
        target_normalization_divisor=target_normalization_divisor,
    )


__all__ = [
    "BRANCH_LENGTH_OPTIMIZATION_FIXED_TOPOLOGY_NNLS",
    "BRANCH_LENGTH_OPTIMIZATION_LINKAGE_ULTRAMETRIC",
    "BRANCH_LENGTH_OPTIMIZATION_METHODS",
    "BRANCH_LENGTH_TARGET_METRICS",
    "BRANCH_LENGTH_TARGET_SQUARED_EUCLIDEAN",
    "BRANCH_LENGTH_TARGET_MEAN_SQUARED_STANDARDIZED_EUCLIDEAN",
    "BranchLengthOptimizationResult",
    "fit_fixed_topology_nnls_branch_lengths",
    "validate_branch_length_optimization_method",
]
