from __future__ import annotations

import networkx as nx
import numpy as np
import pandas as pd
import pytest

from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.pair_testing.collection import (
    record_collection,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.pair_testing.collection.child_parent_edge_metadata import (
    determine_whether_sibling_pair_is_edge_blocked,
    determine_whether_sibling_pair_is_null_like,
    estimate_sibling_null_weight_from_child_parent_edges,
    extract_child_parent_edge_significance_by_node,
    extract_child_parent_edge_testing_status_by_node,
    resolve_sibling_calibration_dependency_group,
)


def test_child_parent_edge_metadata_preserves_tree_node_id_identity() -> None:
    annotations = pd.DataFrame(
        {
            "Child_Parent_Divergence_Significant": [False, True, False],
            "Child_Parent_Divergence_P_Value_BH": [1.0, 0.03, 0.8],
            "Child_Parent_Divergence_Tested": [True, False, True],
            "Child_Parent_Divergence_Ancestor_Blocked": [False, True, False],
        },
        index=[0, 1, 2],
    )

    significance_by_node = extract_child_parent_edge_significance_by_node(annotations)
    tested_by_node, blocked_by_node = extract_child_parent_edge_testing_status_by_node(annotations)

    assert significance_by_node == {0: False, 1: True, 2: False}
    assert tested_by_node == {0: True, 1: False, 2: True}
    assert blocked_by_node == {0: False, 1: True, 2: False}

    assert (
        determine_whether_sibling_pair_is_edge_blocked(
            1,
            2,
            child_parent_edge_tested_by_node=tested_by_node,
            child_parent_edge_ancestor_blocked_by_node=blocked_by_node,
        )
        is True
    )
    assert (
        determine_whether_sibling_pair_is_null_like(
            1,
            2,
            child_parent_edge_significance_by_node=significance_by_node,
        )
        is False
    )


def test_sibling_null_weight_uses_joint_child_edge_weights() -> None:
    edge_p_values = {"L": 0.25, "R": 0.8}
    edge_significant = {"L": False, "R": False}
    edge_tested = {"L": True, "R": True}
    edge_blocked = {"L": False, "R": False}

    weight = estimate_sibling_null_weight_from_child_parent_edges(
        "L",
        "R",
        child_parent_edge_p_values_by_node=edge_p_values,
        child_parent_edge_significance_by_node=edge_significant,
        child_parent_edge_tested_by_node=edge_tested,
        child_parent_edge_ancestor_blocked_by_node=edge_blocked,
    )

    assert weight == 0.25 * 0.8


def test_sibling_null_weight_keeps_small_significant_child_edge_evidence() -> None:
    weight = estimate_sibling_null_weight_from_child_parent_edges(
        "L",
        "R",
        child_parent_edge_p_values_by_node={"L": 0.01, "R": 0.8},
        child_parent_edge_significance_by_node={"L": True, "R": False},
        child_parent_edge_tested_by_node={"L": True, "R": True},
        child_parent_edge_ancestor_blocked_by_node={"L": False, "R": False},
    )

    assert weight == 0.01 * 0.8


def test_sibling_null_weight_maps_significant_zero_p_value_to_positive_underflow_bound() -> None:
    weight = estimate_sibling_null_weight_from_child_parent_edges(
        "L",
        "R",
        child_parent_edge_p_values_by_node={"L": 0.0, "R": 0.8},
        child_parent_edge_significance_by_node={"L": True, "R": False},
        child_parent_edge_tested_by_node={"L": True, "R": True},
        child_parent_edge_ancestor_blocked_by_node={"L": False, "R": False},
    )

    assert 0.0 < weight < 1e-160


def test_sibling_null_weight_joint_product_is_underflow_stable() -> None:
    weight = estimate_sibling_null_weight_from_child_parent_edges(
        "L",
        "R",
        child_parent_edge_p_values_by_node={"L": 0.0, "R": 0.0},
        child_parent_edge_significance_by_node={"L": True, "R": True},
        child_parent_edge_tested_by_node={"L": True, "R": True},
        child_parent_edge_ancestor_blocked_by_node={"L": False, "R": False},
    )

    assert weight > 0.0


def test_sibling_null_weight_treats_tree_bh_stopped_edges_as_ancestor_null_evidence() -> None:
    weight = estimate_sibling_null_weight_from_child_parent_edges(
        "L",
        "R",
        child_parent_edge_p_values_by_node={"L": float("nan"), "R": 0.8},
        child_parent_edge_significance_by_node={"L": False, "R": False},
        child_parent_edge_tested_by_node={"L": False, "R": True},
        child_parent_edge_ancestor_blocked_by_node={"L": True, "R": False},
    )

    assert weight == 0.8


def test_nested_blocked_parents_resolve_to_the_same_stopping_event() -> None:
    tree = nx.DiGraph(
        [
            ("root", "A"),
            ("root", "B"),
            ("A", "C"),
            ("A", "D"),
            ("C", "E"),
            ("C", "F"),
        ]
    )
    tested = {"A": True, "B": True, "C": False, "D": False, "E": False, "F": False}
    significant = {
        "A": False,
        "B": False,
        "C": False,
        "D": False,
        "E": False,
        "F": False,
    }

    direct_group = resolve_sibling_calibration_dependency_group(
        tree,
        "A",
        is_edge_blocked=True,
        child_parent_edge_tested_by_node=tested,
        child_parent_edge_significance_by_node=significant,
    )
    nested_group = resolve_sibling_calibration_dependency_group(
        tree,
        "C",
        is_edge_blocked=True,
        child_parent_edge_tested_by_node=tested,
        child_parent_edge_significance_by_node=significant,
    )

    assert direct_group == "root"
    assert nested_group == "root"


def test_singleton_peel_chain_records_share_one_stopping_event(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # root and C are singleton peels, so Tree-BH tests {x, C, y, F, f1, f2} as one family.
    tree = nx.DiGraph([("root", "x"), ("root", "C"), ("C", "y"), ("C", "F")])
    tree.add_edges_from([("F", "f1"), ("F", "f2")])
    leaf_counts = {"root": 4, "C": 3, "F": 2, "x": 1, "y": 1, "f1": 1, "f2": 1}
    for node, leaf_count in leaf_counts.items():
        tree.nodes[node]["distribution"] = np.array([0.5])
        tree.nodes[node]["leaf_count"] = leaf_count
    annotations = pd.DataFrame(
        {
            "Child_Parent_Divergence_Test_Statistic": 0.0,
            "Child_Parent_Divergence_P_Value": 0.9,
            "Child_Parent_Divergence_P_Value_BH": 0.9,
            "Child_Parent_Divergence_Significant": False,
            "Child_Parent_Divergence_df": 1.0,
            "Child_Parent_Divergence_Tested": True,
            "Child_Parent_Divergence_Ancestor_Blocked": False,
        },
        index=list(leaf_counts),
    )
    monkeypatch.setattr(
        record_collection, "sibling_divergence_test", lambda *a, **k: (1, 1, 1, 0.5)
    )
    parents = ["root", "C", "F"]

    records, _ = record_collection.collect_sibling_pair_records(
        tree,
        annotations,
        sibling_projection_dimensions_from_edge_comparisons=dict.fromkeys(parents, 1),
        parent_principal_component_projections={parent: np.eye(1) for parent in parents},
        parent_principal_component_eigenvalues={parent: np.ones(1) for parent in parents},
    )

    assert {record.parent: record.calibration_dependency_group for record in records} == {
        "root": "root",
        "C": "root",
        "F": "root",
    }
