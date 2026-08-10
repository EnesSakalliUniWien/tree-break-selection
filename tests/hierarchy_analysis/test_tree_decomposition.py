"""Contracts for live tree-decomposition traversal."""

from __future__ import annotations

import networkx as nx
import pandas as pd
from tree_break_selection.hierarchy_analysis.tree_decomposition import TreeDecomposition

from .gate_support import (
    _make_annotations,
    _make_binary_tree,
    _make_deep_tree,
)


def _decompose_with_annotations(
    tree: nx.DiGraph,
    annotations_df: pd.DataFrame,
    *,
    passthrough: bool,
) -> dict[str, object]:
    decomposer = TreeDecomposition(
        tree=tree,
        annotations_df=annotations_df,
        passthrough=passthrough,
        trace_level="full",
    )
    return decomposer.decompose_tree()


def _cluster_leaf_sets(decomposition_results: dict[str, object]) -> list[set[str]]:
    cluster_assignments = decomposition_results["cluster_assignments"]
    return [set(cluster["leaves"]) for cluster in cluster_assignments.values()]


class TestTreeDecompositionTraversal:
    def test_decompose_tree_splits_to_leaf_clusters(self) -> None:
        tree = _make_binary_tree()
        annotations_df = _make_annotations(
            tree,
            edge_divergent={node: True for node in tree.nodes},
            sibling_different={node: True for node in tree.nodes},
        )

        result = _decompose_with_annotations(
            tree,
            annotations_df,
            passthrough=False,
        )

        cluster_leaf_sets = sorted(_cluster_leaf_sets(result), key=lambda leaves: min(leaves))
        assert cluster_leaf_sets == [{"L1"}, {"L2"}, {"R1"}, {"R2"}]

    def test_decompose_tree_default_trace_level_is_compact(self) -> None:
        tree = _make_binary_tree()
        annotations_df = _make_annotations(
            tree,
            edge_divergent={node: True for node in tree.nodes},
            sibling_different={node: False for node in tree.nodes},
        )

        decomposer = TreeDecomposition(
            tree=tree,
            annotations_df=annotations_df,
            passthrough=False,
        )

        result = decomposer.decompose_tree()

        assert "traversal_trace" not in result
        assert "full_edge_traversal_trace" not in result
        assert result["traversal_counters"] == {
            "live_nodes_visited": 1,
            "live_internal_tuples": 1,
            "live_split_count": 0,
            "live_pass_through_count": 0,
            "live_boundary_count": 1,
            "live_passthrough_candidate_count": 0,
            "live_passthrough_support_blocked_count": 0,
        }

    def test_decompose_tree_passthrough_reaches_descendant_split(self) -> None:
        tree = _make_deep_tree()
        edge_divergent = {node: True for node in tree.nodes}
        sibling_different = {node: False for node in tree.nodes}
        sibling_different["B"] = True

        annotations_df = _make_annotations(
            tree,
            edge_divergent=edge_divergent,
            sibling_different=sibling_different,
        )

        result = _decompose_with_annotations(
            tree,
            annotations_df,
            passthrough=True,
        )

        cluster_leaf_sets = sorted(_cluster_leaf_sets(result), key=lambda leaves: min(leaves))
        assert cluster_leaf_sets == [{"A1", "A2"}, {"C1", "C2"}, {"D1", "D2"}]
        root_trace = result["traversal_trace"][0]
        assert root_trace["passthrough_candidate"] is True
        assert root_trace["passthrough_decision_reason"] == "pass_through"
        assert result["traversal_counters"]["live_passthrough_candidate_count"] == 1

    def test_decompose_tree_without_passthrough_merges_at_root(self) -> None:
        tree = _make_deep_tree()
        edge_divergent = {node: True for node in tree.nodes}
        sibling_different = {node: False for node in tree.nodes}
        sibling_different["B"] = True

        annotations_df = _make_annotations(
            tree,
            edge_divergent=edge_divergent,
            sibling_different=sibling_different,
        )

        result = _decompose_with_annotations(
            tree,
            annotations_df,
            passthrough=False,
        )

        cluster_leaf_sets = _cluster_leaf_sets(result)
        assert cluster_leaf_sets == [{"A1", "A2", "C1", "C2", "D1", "D2"}]

    def test_full_edge_traversal_walks_past_sibling_closed_root(self) -> None:
        tree = _make_deep_tree()
        edge_divergent = {node: True for node in tree.nodes}
        sibling_different = {node: False for node in tree.nodes}
        sibling_different["B"] = True

        annotations_df = _make_annotations(
            tree,
            edge_divergent=edge_divergent,
            sibling_different=sibling_different,
        )

        result = _decompose_with_annotations(
            tree,
            annotations_df,
            passthrough=False,
        )

        live_trace = result["traversal_trace"]
        full_trace = result["full_edge_traversal_trace"]
        assert [row["node_id"] for row in live_trace] == ["root"]
        assert {row["node_id"] for row in full_trace} == set(tree.nodes)
        assert full_trace[0]["node_id"] == "root"
        assert full_trace[0]["actual_decision"] == "boundary"
        assert full_trace[0]["actual_visited"] is True
        assert full_trace[0]["edge_traversal_action"] == "continue"
        assert any(row["node_id"] == "B" and not row["actual_visited"] for row in full_trace)

    def test_full_edge_traversal_stops_when_child_edges_close(self) -> None:
        tree = _make_deep_tree()
        edge_divergent = {node: False for node in tree.nodes}
        edge_divergent["A"] = True
        edge_divergent["B"] = True
        sibling_different = {node: False for node in tree.nodes}

        annotations_df = _make_annotations(
            tree,
            edge_divergent=edge_divergent,
            sibling_different=sibling_different,
        )

        result = _decompose_with_annotations(
            tree,
            annotations_df,
            passthrough=False,
        )

        full_trace = result["full_edge_traversal_trace"]
        assert [row["node_id"] for row in full_trace] == ["root", "A", "B"]
        stop_reasons = {row["node_id"]: row["edge_traversal_stop_reason"] for row in full_trace}
        assert stop_reasons == {
            "root": "edge_open_continue",
            "A": "edge_closed",
            "B": "edge_closed",
        }

    def test_full_edge_traversal_records_branch_lengths_and_counters(self) -> None:
        tree = _make_deep_tree()
        tree.edges["root", "A"]["branch_length"] = 1.25
        tree.edges["root", "B"]["branch_length"] = 2.5
        edge_divergent = {node: False for node in tree.nodes}
        edge_divergent["A"] = True
        edge_divergent["B"] = True
        sibling_different = {node: False for node in tree.nodes}

        annotations_df = _make_annotations(
            tree,
            edge_divergent=edge_divergent,
            sibling_different=sibling_different,
        )

        result = _decompose_with_annotations(
            tree,
            annotations_df,
            passthrough=False,
        )

        root_row = result["full_edge_traversal_trace"][0]
        assert root_row["left_child"] == "A"
        assert root_row["right_child"] == "B"
        assert root_row["left_edge_test_tuple"] == ("root", "A")
        assert root_row["right_edge_test_tuple"] == ("root", "B")
        assert root_row["sibling_test_tuple"] == ("root", "A", "B")
        assert root_row["left_edge_p_value"] == 0.01
        assert root_row["left_edge_p_value_bh"] == 0.01
        assert root_row["left_edge_tested"] is True
        assert root_row["right_edge_p_value"] == 0.01
        assert root_row["right_edge_p_value_bh"] == 0.01
        assert root_row["right_edge_tested"] is True
        assert root_row["sibling_p_value"] == 0.80
        assert root_row["sibling_p_value_corrected"] == 0.80
        assert root_row["sibling_test_statistic"] == 0.5
        assert root_row["sibling_degrees_of_freedom"] == 1.0
        assert root_row["sibling_test_method"] == "synthetic_projected_wald"
        assert root_row["descendant_leaf_signature"] == (
            "A1",
            "A2",
            "C1",
            "C2",
            "D1",
            "D2",
        )
        assert root_row["passthrough_enabled"] is False
        assert root_row["passthrough_candidate"] is False
        assert root_row["passthrough_decision_reason"] == "passthrough_disabled"
        assert root_row["left_branch_length"] == 1.25
        assert root_row["right_branch_length"] == 2.5
        assert root_row["left_branch_length_missing"] is False
        assert root_row["right_branch_length_missing"] is False

        counters = result["traversal_counters"]
        assert counters["live_nodes_visited"] == len(result["traversal_trace"])
        assert counters["live_internal_tuples"] == 1
        assert counters["live_boundary_count"] == 1
        assert counters["full_edge_nodes_visited"] == len(result["full_edge_traversal_trace"])
        assert counters["full_edge_internal_tuples"] == 3
        assert counters["full_edge_continue_count"] == 1
        assert counters["full_edge_stop_count"] == 2
        assert counters["full_edge_closed_stop_count"] == 2
        assert counters["full_edge_passthrough_candidate_count"] == 0

    def test_decompose_tree_passthrough_ignores_blocked_descendant_signal(self) -> None:
        tree = _make_deep_tree()
        edge_divergent = {node: True for node in tree.nodes}
        edge_divergent["C"] = False
        edge_divergent["D"] = False
        sibling_different = {node: False for node in tree.nodes}
        sibling_different["B"] = True

        annotations_df = _make_annotations(
            tree,
            edge_divergent=edge_divergent,
            sibling_different=sibling_different,
        )

        result = _decompose_with_annotations(
            tree,
            annotations_df,
            passthrough=True,
        )

        cluster_leaf_sets = _cluster_leaf_sets(result)
        assert cluster_leaf_sets == [{"A1", "A2", "C1", "C2", "D1", "D2"}]

    def test_decompose_tree_selected_family_guard_blocks_passthrough(self) -> None:
        tree = _make_deep_tree()
        edge_divergent = {node: True for node in tree.nodes}
        sibling_different = {node: False for node in tree.nodes}
        sibling_different["B"] = True
        annotations_df = _make_annotations(
            tree,
            edge_divergent=edge_divergent,
            sibling_different=sibling_different,
        )
        annotations_df["Selective_Permutation_Guard_Would_Block"] = False
        annotations_df["Selective_Permutation_Guard_Blocked"] = False
        annotations_df.loc["root", "Selective_Permutation_Guard_Would_Block"] = True

        decomposer = TreeDecomposition(
            tree=tree,
            annotations_df=annotations_df,
            passthrough=True,
            selected_family_passthrough_guard=True,
            trace_level="full",
        )

        result = decomposer.decompose_tree()

        assert _cluster_leaf_sets(result) == [{"A1", "A2", "C1", "C2", "D1", "D2"}]
        root_trace = result["traversal_trace"][0]
        assert root_trace["decision"] == "boundary"
        assert root_trace["passthrough_candidate"] is True
        assert root_trace["passthrough_supported"] is False
        assert root_trace["passthrough_bottleneck"] == ("selected_family_passthrough_guard_blocked")
        assert root_trace["passthrough_decision_reason"] == "passthrough_support_blocked"
        assert result["traversal_counters"]["live_passthrough_support_blocked_count"] == 1
