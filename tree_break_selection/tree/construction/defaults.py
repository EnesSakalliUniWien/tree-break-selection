"""Tree-construction method names owned by the tree module.

These values describe the supported construction seams. They are constants,
not mutable runtime configuration; callers with a different geometry or
topology policy must pass it explicitly.
"""

from __future__ import annotations

from typing import Literal

TreeLinkageMethod = Literal[
    "average",
    "complete",
    "weighted",
    "single",
    "centroid",
    "median",
    "ward",
]

DEFAULT_BINARY_TREE_DISTANCE_METRIC = "hamming"
DEFAULT_TREE_LINKAGE_METHOD = "average"

SUPPORTED_LINKAGE_METHODS: tuple[TreeLinkageMethod, ...] = (
    "average",
    "complete",
    "weighted",
    "single",
    "centroid",
    "median",
    "ward",
)
TREE_CONSENSUS_PRIORITY_ORDER: tuple[str, ...] = (
    "weighted",
    "ward",
    "complete",
    "average",
    "centroid",
    "median",
    "single",
)
TREE_CONSENSUS_STRATEGY_PRIORITY: dict[str, int] = {
    strategy: index for index, strategy in enumerate(TREE_CONSENSUS_PRIORITY_ORDER)
}

__all__ = [
    "DEFAULT_BINARY_TREE_DISTANCE_METRIC",
    "DEFAULT_TREE_LINKAGE_METHOD",
    "SUPPORTED_LINKAGE_METHODS",
    "TREE_CONSENSUS_PRIORITY_ORDER",
    "TREE_CONSENSUS_STRATEGY_PRIORITY",
    "TreeLinkageMethod",
]
