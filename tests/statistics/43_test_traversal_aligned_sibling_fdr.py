from __future__ import annotations

import networkx as nx
import numpy as np
import pandas as pd
import pytest

from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.inflated_projected_wald_annotation.fdr_annotation import (
    apply_traversal_aligned_sibling_bh_results,
    init_sibling_annotation_df,
)


def _binary_three_level_tree() -> nx.DiGraph:
    tree = nx.DiGraph()
    tree.add_edge("root", "A")
    tree.add_edge("root", "B")
    tree.add_edge("A", "A1")
    tree.add_edge("A", "A2")
    tree.add_edge("B", "B1")
    tree.add_edge("B", "B2")
    return tree


def _annotations(
    *,
    root_children_significant: bool = True,
    a_children_significant: bool = True,
    b_children_significant: bool = True,
) -> pd.DataFrame:
    edge_flags = {
        "root": False,
        "A": root_children_significant,
        "B": root_children_significant,
        "A1": a_children_significant,
        "A2": a_children_significant,
        "B1": b_children_significant,
        "B2": b_children_significant,
    }
    annotations_df = pd.DataFrame(
        {
            "Child_Parent_Divergence_Significant": pd.Series(
                edge_flags,
                dtype=bool,
            )
        }
    )
    return init_sibling_annotation_df(annotations_df)


def test_root_sibling_test_is_not_penalized_by_unreached_descendants() -> None:
    tree = _binary_three_level_tree()
    annotations_df = _annotations()

    out = apply_traversal_aligned_sibling_bh_results(
        tree,
        annotations_df,
        ["root", "A", "B"],
        [(10.0, 1.0, 0.006), (1.0, 1.0, 0.90), (1.0, 1.0, 0.90)],
        0.01,
    )

    assert bool(out.loc["root", "Sibling_BH_Different"])
    assert np.isclose(out.loc["root", "Sibling_Divergence_P_Value_Corrected"], 0.006)
    assert not bool(out.loc["A", "Sibling_BH_Different"])
    assert not bool(out.loc["B", "Sibling_BH_Different"])


def test_descendants_remain_testable_for_passthrough_when_parent_sibling_is_same() -> None:
    tree = _binary_three_level_tree()
    annotations_df = _annotations()

    out = apply_traversal_aligned_sibling_bh_results(
        tree,
        annotations_df,
        ["root", "A", "B"],
        [(1.0, 1.0, 0.90), (10.0, 1.0, 0.001), (1.0, 1.0, 0.90)],
        0.01,
    )

    assert not bool(out.loc["root", "Sibling_BH_Different"])
    assert bool(out.loc["root", "Sibling_BH_Same"])
    assert bool(out.loc["A", "Sibling_BH_Different"])
    assert np.isclose(out.loc["A", "Sibling_Divergence_P_Value_Corrected"], 0.002)
    assert not bool(out.loc["B", "Sibling_BH_Different"])


def test_nodes_with_closed_edge_gate_keep_raw_sibling_stat_but_no_fdr_decision() -> None:
    tree = _binary_three_level_tree()
    annotations_df = _annotations(a_children_significant=False)

    out = apply_traversal_aligned_sibling_bh_results(
        tree,
        annotations_df,
        ["root", "A", "B"],
        [(10.0, 1.0, 0.001), (20.0, 1.0, 0.0), (10.0, 1.0, 0.001)],
        0.01,
    )

    assert bool(out.loc["root", "Sibling_BH_Different"])
    assert np.isclose(out.loc["A", "Sibling_Divergence_P_Value"], 0.0)
    assert pd.isna(out.loc["A", "Sibling_Divergence_P_Value_Corrected"])
    assert not bool(out.loc["A", "Sibling_BH_Different"])
    assert not bool(out.loc["A", "Sibling_BH_Same"])
    assert bool(out.loc["B", "Sibling_BH_Different"])


def test_frontier_passes_through_a_null_like_singleton_peel() -> None:
    tree = nx.DiGraph([("root", "leaf"), ("root", "C"), ("C", "C1"), ("C", "C2")])
    edge_flags = {"root": False, "leaf": False, "C": False, "C1": True, "C2": True}
    annotations_df = init_sibling_annotation_df(
        pd.DataFrame({"Child_Parent_Divergence_Significant": pd.Series(edge_flags, dtype=bool)})
    )

    out = apply_traversal_aligned_sibling_bh_results(
        tree, annotations_df, ["C"], [(10.0, 1.0, 0.001)], 0.01, skipped_parents=["root"]
    )

    assert bool(out.loc["C", "Sibling_BH_Different"])


def test_frontier_passes_through_a_null_like_small_peel_only_with_a_leaf_limit() -> None:
    tree = nx.DiGraph(
        [("root", "S"), ("root", "C"), ("S", "s1"), ("S", "s2"), ("C", "C1"), ("C", "C2")]
        + [("C1", "a1"), ("C1", "a2"), ("C2", "b1"), ("C2", "b2")]
    )
    edge_flags = {node: node in {"C1", "C2"} for node in tree.nodes}

    def run(**kwargs: int) -> pd.DataFrame:
        annotations_df = init_sibling_annotation_df(
            pd.DataFrame({"Child_Parent_Divergence_Significant": pd.Series(edge_flags, dtype=bool)})
        )
        return apply_traversal_aligned_sibling_bh_results(
            tree,
            annotations_df,
            ["C"],
            [(10.0, 1.0, 0.001)],
            0.01,
            skipped_parents=["root"],
            **kwargs,
        )

    # S has 2 leaves (< 3) and C has 4: with the limit the frontier continues in C.
    assert not bool(run().loc["C", "Sibling_BH_Different"])
    assert bool(run(small_peel_leaf_limit=3).loc["C", "Sibling_BH_Different"])


def test_frontier_passes_through_a_closed_parent_without_testing_it() -> None:
    tree = _binary_three_level_tree()
    annotations_df = _annotations()

    out = apply_traversal_aligned_sibling_bh_results(
        tree,
        annotations_df,
        ["A", "B"],
        [(10.0, 1.0, 0.001), (1.0, 1.0, 0.90)],
        0.01,
        closed_parents=["root"],
    )

    assert not bool(out.loc["root", "Sibling_BH_Different"])
    assert pd.isna(out.loc["root", "Sibling_Divergence_P_Value_Corrected"])
    # Only A and B form the depth-1 family: 0.001 * 2 / 1.
    assert np.isclose(out.loc["A", "Sibling_Divergence_P_Value_Corrected"], 0.002)
    assert bool(out.loc["A", "Sibling_BH_Different"])


@pytest.mark.parametrize(
    ("method_label", "expected_calibration"),
    [
        ("fixed_global_chi_square", "fixed_subspace_chi_square"),
        ("fixed_coordinate_bh", "fixed_subspace_bh"),
        ("fixed_coordinate_by", "fixed_subspace_by"),
        ("fixed_coordinate_bonferroni", "fixed_subspace_bonferroni"),
        ("fixed_block_bh", "fixed_subspace_bh"),
        ("fixed_block_simes_bh", "fixed_subspace_simes_bh"),
        ("global_chi_square", "data_independent_chi_square"),
        ("coordinate_bonferroni", "data_independent_coordinate_bonferroni"),
        ("coordinate_bh", "data_independent_coordinate_bh"),
        ("block_bonferroni", "data_independent_block_bonferroni"),
        ("block_bh", "data_independent_block_bh"),
        ("zero_dimensional_sibling_record", "degenerate_zero_dimensional"),
        (
            "context_weighted_supported_empirical_null_inflation:local_kernel",
            "empirical_null_inflation",
        ),
    ],
)
def test_sibling_gate_calibration_is_explicit_for_supported_methods(
    method_label: str,
    expected_calibration: str,
) -> None:
    tree = _binary_three_level_tree()
    annotations_df = _annotations()

    out = apply_traversal_aligned_sibling_bh_results(
        tree,
        annotations_df,
        ["root"],
        [(1.0, 1.0, 0.5)],
        0.01,
        method_labels=[method_label],
    )

    assert out.loc["root", "Sibling_Gate_P_Value_Calibration"] == expected_calibration


def test_unknown_sibling_gate_method_is_rejected() -> None:
    tree = _binary_three_level_tree()
    annotations_df = _annotations()

    with pytest.raises(ValueError, match="Unsupported sibling test method label"):
        apply_traversal_aligned_sibling_bh_results(
            tree,
            annotations_df,
            ["root"],
            [(1.0, 1.0, 0.5)],
            0.01,
            method_labels=["unknown_method"],
        )
