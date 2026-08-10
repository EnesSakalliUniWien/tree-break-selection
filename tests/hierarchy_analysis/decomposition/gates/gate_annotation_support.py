"""Shared test data for gate profiles, guards, and orchestration."""

from __future__ import annotations

import networkx as nx
import numpy as np
import pandas as pd


def _build_small_tree_with_leaf_data() -> tuple[nx.DiGraph, pd.DataFrame, pd.DataFrame]:
    tree = nx.DiGraph()
    tree.graph["root"] = "root"
    tree.add_edge("root", "A", branch_length=0.25)
    tree.add_edge("root", "B", branch_length=0.20)
    tree.add_edge("cal", "C", branch_length=0.10)
    tree.add_edge("cal", "D", branch_length=0.10)

    root_dist = np.array([0.50, 0.50, 0.50, 0.50, 0.50, 0.50], dtype=np.float64)
    a_dist = np.array([0.12, 0.12, 0.12, 0.88, 0.88, 0.88], dtype=np.float64)
    b_dist = np.array([0.88, 0.88, 0.88, 0.12, 0.12, 0.12], dtype=np.float64)
    cal_dist = np.array([0.50, 0.50, 0.50, 0.50, 0.50, 0.50], dtype=np.float64)
    c_dist = np.array([0.49, 0.49, 0.49, 0.51, 0.51, 0.51], dtype=np.float64)
    d_dist = np.array([0.51, 0.51, 0.51, 0.49, 0.49, 0.49], dtype=np.float64)

    for node, dist, leaf_count, is_leaf in (
        ("root", root_dist, 200, False),
        ("A", a_dist, 100, True),
        ("B", b_dist, 100, True),
        ("cal", cal_dist, 200, False),
        ("C", c_dist, 100, True),
        ("D", d_dist, 100, True),
    ):
        tree.nodes[node]["distribution"] = dist
        tree.nodes[node]["leaf_count"] = leaf_count
        tree.nodes[node]["is_leaf"] = is_leaf
        tree.nodes[node]["label"] = node

    annotations_df = pd.DataFrame(
        {
            "leaf_count": {
                "root": 200,
                "A": 100,
                "B": 100,
                "cal": 200,
                "C": 100,
                "D": 100,
            }
        }
    )
    leaf_data = pd.DataFrame(
        [
            [0, 0, 0, 1, 1, 1],
            [1, 1, 1, 0, 0, 0],
            [0.49, 0.49, 0.49, 0.51, 0.51, 0.51],
            [0.51, 0.51, 0.51, 0.49, 0.49, 0.49],
        ],
        index=["A", "B", "C", "D"],
        dtype=np.float64,
    )
    return tree, annotations_df, leaf_data


def _selected_root(selected_p_value: float):
    def runner(*_args, **_kwargs):
        return {
            "root_observed_p_value": 0.001,
            "root_selective_p_value": selected_p_value,
            "root_selective_null_min_p_value": 0.01,
            "root_selective_null_q05_p_value": 0.05,
        }

    return runner
