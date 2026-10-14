"""Targeted tests for TreeDecomposition annotation handling."""

from __future__ import annotations

import pandas as pd
from tree_break_selection.hierarchy_analysis.tree_decomposition import TreeDecomposition
from tree_break_selection.tree.poset_tree import PosetTree


def test_tree_decomposition_preserves_non_string_node_ids_in_annotations():
    tree = PosetTree()
    tree.add_node(0, is_leaf=False)
    tree.add_node(1, is_leaf=True, label="left")
    tree.add_node(2, is_leaf=True, label="right")
    tree.add_edge(0, 1)
    tree.add_edge(0, 2)

    annotations = pd.DataFrame(
        {
            "Child_Parent_Divergence_Significant": [False, True, True],
            "Sibling_BH_Different": [True, False, False],
            "Sibling_Divergence_Skipped": [False, False, False],
        },
        index=[0, 1, 2],
    )

    result = TreeDecomposition(tree=tree, annotations_df=annotations).decompose_tree()

    assert result["num_clusters"] == 2
    assert [cluster["root_node"] for cluster in result["cluster_assignments"].values()] == [1, 2]
