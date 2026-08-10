"""Shared tree data for gate evaluation, traversal, and spectral transport tests."""

from __future__ import annotations

import networkx as nx
import numpy as np
import pandas as pd
from tree_break_selection.tree.poset_tree import PosetTree


def _annotate_tree_structure(tree: nx.DiGraph, leaves: set[str]) -> None:
    for node in tree.nodes:
        tree.nodes[node]["is_leaf"] = node in leaves
        tree.nodes[node]["label"] = node
        tree.nodes[node]["distribution"] = np.array([0.5], dtype=float)


def _make_binary_tree() -> PosetTree:
    tree = PosetTree()
    tree.add_edges_from(
        [
            ("root", "L"),
            ("root", "R"),
            ("L", "L1"),
            ("L", "L2"),
            ("R", "R1"),
            ("R", "R2"),
        ]
    )
    _annotate_tree_structure(tree, {"L1", "L2", "R1", "R2"})
    return tree


def _make_deep_tree() -> PosetTree:
    tree = PosetTree()
    tree.add_edges_from(
        [
            ("root", "A"),
            ("root", "B"),
            ("A", "A1"),
            ("A", "A2"),
            ("B", "C"),
            ("B", "D"),
            ("C", "C1"),
            ("C", "C2"),
            ("D", "D1"),
            ("D", "D2"),
        ]
    )
    _annotate_tree_structure(tree, {"A1", "A2", "C1", "C2", "D1", "D2"})
    return tree


def _make_annotations(
    tree: nx.DiGraph,
    *,
    edge_divergent: dict[str, bool],
    sibling_different: dict[str, bool],
    sibling_skipped: dict[str, bool] | None = None,
) -> pd.DataFrame:
    if sibling_skipped is None:
        sibling_skipped = {node: False for node in tree.nodes}

    edge_p_values = {node: 0.01 if edge_divergent[node] else 0.90 for node in tree.nodes}
    sibling_p_values = {
        node: 0.02 if sibling_different[node] and not sibling_skipped[node] else 0.80
        for node in tree.nodes
    }

    return pd.DataFrame(
        {
            "Child_Parent_Divergence_Significant": pd.Series(edge_divergent, dtype=bool),
            "Child_Parent_Divergence_P_Value": pd.Series(edge_p_values, dtype=float),
            "Child_Parent_Divergence_P_Value_BH": pd.Series(edge_p_values, dtype=float),
            "Child_Parent_Divergence_Tested": True,
            "Child_Parent_Divergence_Ancestor_Blocked": False,
            "Sibling_BH_Different": pd.Series(sibling_different, dtype=bool),
            "Sibling_Divergence_Skipped": pd.Series(sibling_skipped, dtype=bool),
            "Sibling_Divergence_P_Value": pd.Series(sibling_p_values, dtype=float),
            "Sibling_Divergence_P_Value_Corrected": pd.Series(sibling_p_values, dtype=float),
            "Sibling_Test_Statistic": pd.Series(
                {
                    node: 8.0 if sibling_different[node] and not sibling_skipped[node] else 0.5
                    for node in tree.nodes
                },
                dtype=float,
            ),
            "Sibling_Degrees_of_Freedom": 1.0,
            "Sibling_Test_Method": "synthetic_projected_wald",
        }
    ).reindex(list(tree.nodes))
