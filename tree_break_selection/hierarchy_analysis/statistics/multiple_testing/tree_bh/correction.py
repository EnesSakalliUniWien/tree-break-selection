"""Tree-BH hierarchical FDR correction implementation."""

from __future__ import annotations

from collections import deque

import networkx as nx
import numpy as np

from ..base import benjamini_hochberg_correction
from .helpers import find_root_nodes, group_child_indices_by_parent
from .models import ChildParentEdgeTreeBHResult


def _validate_tree_bh_inputs(
    p_values: np.ndarray,
    child_ids: list[str],
    alpha: float,
) -> tuple[np.ndarray, float]:
    """Validate Tree-BH boundary inputs before hierarchical traversal."""
    p_values_array = np.asarray(p_values, dtype=float)
    alpha_value = float(alpha)

    if p_values_array.ndim != 1:
        raise ValueError(f"Tree-BH p-values must be a 1-D array; got shape {p_values_array.shape}.")
    if p_values_array.shape[0] != len(child_ids):
        raise ValueError(
            "Tree-BH p-values must align one-to-one with child_ids: "
            f"{p_values_array.shape[0]} p-value(s), {len(child_ids)} child id(s)."
        )
    if not np.isfinite(alpha_value) or not (0.0 < alpha_value <= 1.0):
        raise ValueError(f"Tree-BH alpha must be finite and in (0, 1]; got {alpha!r}.")
    if np.any(~np.isfinite(p_values_array)):
        raise ValueError("Tree-BH p-values must be finite before correction.")
    if np.any((p_values_array < 0.0) | (p_values_array > 1.0)):
        raise ValueError("Tree-BH p-values must lie in [0, 1].")

    return p_values_array, alpha_value


def apply_tree_bh_correction(
    tree: nx.DiGraph,
    p_values: np.ndarray,
    child_ids: list[str],
    alpha: float,
) -> ChildParentEdgeTreeBHResult:
    """Apply Tree-BH hierarchical FDR correction.

    The procedure works top-down through the tree. Within each reached sibling
    group, it applies BH at a sibling-group-specific level obtained by scaling
    the base alpha by ancestor rejection proportions. Stored corrected p-values
    are the BH adjustments within each reached family; rejection flags compare
    them with that family level, rather than the base alpha.
    """
    p_values, alpha = _validate_tree_bh_inputs(p_values, child_ids, alpha)
    hypothesis_count = len(p_values)
    rejected = np.zeros(hypothesis_count, dtype=bool)
    corrected_p_values = np.ones(hypothesis_count, dtype=float)
    tested = np.zeros(hypothesis_count, dtype=bool)

    if hypothesis_count:
        indices_by_parent = group_child_indices_by_parent(tree, child_ids)
        # A family's level is its parent family's level times that family's rejection share.
        family_alpha = {root: alpha for root in find_root_nodes(tree)}
        frontier = deque(family_alpha)
        while frontier:
            parent_id = frontier.popleft()
            indices = indices_by_parent.get(parent_id)
            if not indices:
                continue
            family_rejected, family_corrected = benjamini_hochberg_correction(
                p_values[indices], alpha=family_alpha[parent_id]
            )
            rejected[indices] = family_rejected
            corrected_p_values[indices] = family_corrected
            tested[indices] = True
            child_family_alpha = family_alpha[parent_id] * family_rejected.sum() / len(indices)
            for index in np.asarray(indices)[family_rejected]:
                child_id = str(child_ids[index])
                family_alpha[child_id] = child_family_alpha
                frontier.append(child_id)

    return ChildParentEdgeTreeBHResult(
        child_parent_edge_null_rejected_by_tree_bh=rejected,
        child_parent_edge_corrected_p_values_by_tree_bh=corrected_p_values,
        child_parent_edge_tested_by_tree_bh=tested,
    )


__all__ = ["apply_tree_bh_correction"]
