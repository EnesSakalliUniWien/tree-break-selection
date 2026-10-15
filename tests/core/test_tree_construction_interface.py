from __future__ import annotations

import networkx as nx
import numpy as np
import pandas as pd
import pytest
from scipy.spatial.distance import pdist

from tree_break_selection.core_utils.tree_utils import bottom_up_nodes
from tree_break_selection.tree.construction import build_tree
from tree_break_selection.tree.poset_tree import PosetTree


def _data() -> pd.DataFrame:
    return pd.DataFrame(
        [[0.0], [0.1], [2.0], [2.1]],
        index=[10, 11, 20, 21],
        columns=["x"],
    )


def test_linkage_builder_preserves_input_index_labels() -> None:
    data = _data()

    result = build_tree(
        data,
        pdist(data.to_numpy()),
        linkage_method="average",
    )

    assert result.linkage_matrix is not None
    assert result.tree.get_leaves() == [10, 11, 20, 21]
    assert result.diagnostics.leaf_count == 4
    assert result.diagnostics.internal_node_count == 3
    assert result.diagnostics.edge_count == 6
    assert result.diagnostics.root_child_count == 2


def test_linkage_builder_records_exact_distance_tie_burden() -> None:
    data = pd.DataFrame(
        [[0.0], [0.0], [1.0], [1.0]],
        index=["a", "b", "c", "d"],
        columns=["x"],
    )

    result = build_tree(
        data,
        pdist(data.to_numpy()),
        linkage_method="average",
    )

    assert result.diagnostics.distance_pair_count == 6
    assert result.diagnostics.zero_distance_pair_count == 2
    assert result.diagnostics.repeated_distance_pair_count == 6
    assert result.diagnostics.repeated_distance_value_count == 2
    assert result.diagnostics.repeated_linkage_height_count == 2


def test_linkage_canonicalization_uses_natural_numeric_label_order() -> None:
    data = pd.DataFrame(
        [[3.0], [1.0], [0.0], [2.0]],
        index=[10, 2, 1, 3],
        columns=["x"],
    )

    result = build_tree(
        data,
        pdist(data.to_numpy()),
        linkage_method="average",
    )

    assert [
        result.tree.nodes[node]["label"]
        for node in result.tree
        if result.tree.out_degree(node) == 0
    ] == [1, 2, 3, 10]


def test_linkage_canonicalization_uses_natural_numbered_string_order() -> None:
    data = pd.DataFrame(
        [[3.0], [1.0], [0.0], [2.0]],
        index=["S10", "S2", "S1", "S3"],
        columns=["x"],
    )

    result = build_tree(
        data,
        pdist(data.to_numpy()),
        linkage_method="average",
    )

    assert [
        result.tree.nodes[node]["label"]
        for node in result.tree
        if result.tree.out_degree(node) == 0
    ] == ["S1", "S2", "S3", "S10"]


def test_linkage_tie_resolution_is_invariant_to_dataframe_row_order() -> None:
    data = pd.DataFrame(
        [[1.0], [0.0], [1.0], [0.0]],
        index=["d", "a", "c", "b"],
        columns=["x"],
    )
    shuffled = data.loc[["b", "c", "a", "d"]]

    first = build_tree(
        data,
        pdist(data.to_numpy()),
        linkage_method="average",
    )
    second = build_tree(
        shuffled,
        pdist(shuffled.to_numpy()),
        linkage_method="average",
    )

    assert first.tree.get_leaves() == ["a", "b", "c", "d"]
    assert first.tree.compute_descendant_sets(use_labels=True) == (
        second.tree.compute_descendant_sets(use_labels=True)
    )


def test_builder_rejects_non_binary_result_at_construction_seam(monkeypatch) -> None:
    data = pd.DataFrame(
        [[0.0], [1.0], [2.0]],
        index=["a", "b", "c"],
        columns=["x"],
    )
    tree = PosetTree()
    tree.add_node("root", is_leaf=False, label="root")
    for label in data.index:
        tree.add_node(label, is_leaf=True, label=label)
        tree.add_edge("root", label, branch_length=1.0)
    tree.graph["root"] = "root"

    monkeypatch.setattr(
        "tree_break_selection.tree.construction.build.tree_from_linkage",
        lambda *_args, **_kwargs: tree,
    )

    with pytest.raises(ValueError, match="exactly two children"):
        build_tree(
            data,
            np.array([1.0, 2.0, 1.0]),
            linkage_method="average",
        )


def test_bottom_up_order_does_not_depend_on_node_insertion_order() -> None:
    first = nx.DiGraph()
    first.add_edges_from([("root", "right"), ("root", "left")])
    second = nx.DiGraph()
    second.add_edges_from([("root", "left"), ("root", "right")])

    first_order = list(bottom_up_nodes(first))
    second_order = list(bottom_up_nodes(second))

    assert first_order == second_order
    assert first_order[-1] == "root"
