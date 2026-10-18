"""Data models for the Tree-BH hierarchical multiple-testing procedure."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ChildParentEdgeTreeBHResult:
    """Tree-BH correction output for child-parent edge hypotheses.

    Stored corrected p-values are BH adjustments within each reached family.
    Rejection flags use the ancestor-scaled family levels, so comparing these
    values with the base alpha does not reproduce hierarchical rejections.
    Untested edges retain the placeholder value 1.
    """

    child_parent_edge_null_rejected_by_tree_bh: np.ndarray
    child_parent_edge_corrected_p_values_by_tree_bh: np.ndarray
    child_parent_edge_tested_by_tree_bh: np.ndarray


__all__ = ["ChildParentEdgeTreeBHResult"]
