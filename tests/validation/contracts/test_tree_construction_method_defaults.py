from __future__ import annotations

from benchmarks.shared import tree_consensus
from benchmarks.validation.sweeps import family_metric_nnls_grid
from tree_break_selection.tree.construction import (
    DEFAULT_TREE_LINKAGE_METHOD,
    SUPPORTED_LINKAGE_METHODS,
)


def test_tree_strategies_use_the_supported_linkage_methods() -> None:
    assert set(tree_consensus.TREE_CONSENSUS_STRATEGY_PRIORITY) == set(SUPPORTED_LINKAGE_METHODS)
    assert family_metric_nnls_grid.parse_args([]).topologies == SUPPORTED_LINKAGE_METHODS


def test_family_metric_topology_parameters_use_the_linkage_method() -> None:
    assert family_metric_nnls_grid._topology_parameters(DEFAULT_TREE_LINKAGE_METHOD) == {
        "tree_linkage_method": DEFAULT_TREE_LINKAGE_METHOD,
    }
