from __future__ import annotations

import networkx as nx
import pandas as pd
from scipy.spatial.distance import pdist

from tree_break_selection.tree.construction import build_tree
from tree_break_selection.tree.topology import bottom_up_nodes


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


def test_bottom_up_order_does_not_depend_on_node_insertion_order() -> None:
    first = nx.DiGraph()
    first.add_edges_from([("root", "right"), ("root", "left")])
    second = nx.DiGraph()
    second.add_edges_from([("root", "left"), ("root", "right")])

    first_order = list(bottom_up_nodes(first))
    second_order = list(bottom_up_nodes(second))

    assert first_order == second_order
    assert first_order[-1] == "root"
