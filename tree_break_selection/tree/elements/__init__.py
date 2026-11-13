"""Set based identity for tree nodes: ``Partition``, ``PartitionSet`` and split comparison.

Copied from BranchArchitect (``brancharchitect/elements`` at commit a9c8829, MIT, same
author); only the import paths and the ``Self`` import differ. ``FrozenPartitionSet`` was
not copied.
"""

from tree_break_selection.tree.elements.partition import Partition
from tree_break_selection.tree.elements.partition_set import PartitionSet
from tree_break_selection.tree.elements.split_semantics import (
    same_rooted_split,
    same_unrooted_split_within_pivot,
    unrooted_difference_within_pivot,
)

__all__ = [
    "Partition",
    "PartitionSet",
    "same_rooted_split",
    "same_unrooted_split_within_pivot",
    "unrooted_difference_within_pivot",
]
