"""Tree-construction interface grouped by methodological responsibility."""

from tree_break_selection.tree.construction.build import (
    TreeBuildDiagnostics,
    TreeBuildResult,
    build_tree,
)
from tree_break_selection.tree.construction.defaults import (
    DEFAULT_BINARY_TREE_DISTANCE_METRIC,
    DEFAULT_TREE_LINKAGE_METHOD,
    SUPPORTED_LINKAGE_METHODS,
    TREE_CONSENSUS_STRATEGY_PRIORITY,
    TreeLinkageMethod,
)
from tree_break_selection.tree.construction.hierarchical import (
    tree_from_linkage,
)

__all__ = [
    "DEFAULT_BINARY_TREE_DISTANCE_METRIC",
    "DEFAULT_TREE_LINKAGE_METHOD",
    "SUPPORTED_LINKAGE_METHODS",
    "TREE_CONSENSUS_STRATEGY_PRIORITY",
    "TreeBuildDiagnostics",
    "TreeBuildResult",
    "TreeLinkageMethod",
    "build_tree",
    "tree_from_linkage",
]
