from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from tree_break_selection.hierarchy_analysis.decomposition.gates.orchestrator import (
    run_gate_annotation_pipeline,
)
from tree_break_selection.hierarchy_analysis.statistics.child_parent_divergence import (
    annotate_child_parent_divergence,
)
from tree_break_selection.hierarchy_analysis.statistics.projection.spectral.spectral_decomposition_result import (
    SpectralDecompositionResult,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.projection.gate_inputs.projection_dimensions import (
    derive_sibling_projection_dimensions_from_child_edge_comparisons,
)
from tree_break_selection.tree.poset_tree import PosetTree


def _build_cherry_tree() -> tuple[PosetTree, pd.DataFrame, pd.DataFrame]:
    tree = PosetTree()
    tree.add_node(
        "top",
        is_leaf=False,
        distribution=np.array([0.20, 0.20, 0.20, 0.80, 0.80, 0.80], dtype=float),
        label="top",
        leaf_count=400,
    )
    tree.add_node(
        "root",
        is_leaf=False,
        distribution=np.array([0.50, 0.50, 0.50, 0.50, 0.50, 0.50], dtype=float),
        label="root",
        leaf_count=200,
    )
    tree.add_node(
        "A",
        is_leaf=True,
        distribution=np.array([0.12, 0.12, 0.12, 0.88, 0.88, 0.88], dtype=float),
        label="A",
        leaf_count=100,
    )
    tree.add_node(
        "B",
        is_leaf=True,
        distribution=np.array([0.88, 0.88, 0.88, 0.12, 0.12, 0.12], dtype=float),
        label="B",
        leaf_count=100,
    )
    tree.add_node(
        "cal",
        is_leaf=False,
        distribution=np.array([0.50, 0.50, 0.50, 0.50, 0.50, 0.50], dtype=float),
        label="cal",
        leaf_count=200,
    )
    tree.add_node(
        "C",
        is_leaf=True,
        distribution=np.array([0.495, 0.495, 0.495, 0.505, 0.505, 0.505], dtype=float),
        label="C",
        leaf_count=100,
    )
    tree.add_node(
        "D",
        is_leaf=True,
        distribution=np.array([0.505, 0.505, 0.505, 0.495, 0.495, 0.495], dtype=float),
        label="D",
        leaf_count=100,
    )
    tree.add_edge("top", "root", branch_length=0.10)
    tree.add_edge("top", "cal", branch_length=0.10)
    tree.add_edge("root", "A", branch_length=0.25)
    tree.add_edge("root", "B", branch_length=0.20)
    tree.add_edge("cal", "C", branch_length=0.10)
    tree.add_edge("cal", "D", branch_length=0.10)
    tree.graph["root"] = "top"

    annotations_df = pd.DataFrame(
        {
            "leaf_count": {
                "top": 400,
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
            [0.495, 0.495, 0.495, 0.505, 0.505, 0.505],
            [0.505, 0.505, 0.505, 0.495, 0.495, 0.495],
        ],
        index=["A", "B", "C", "D"],
        dtype=float,
    )
    return tree, annotations_df, leaf_data


def _build_mixed_tree() -> tuple[PosetTree, pd.DataFrame, pd.DataFrame]:
    tree = PosetTree()
    tree.add_node(
        "top",
        is_leaf=False,
        distribution=np.array([0.20, 0.20, 0.20, 0.80, 0.80, 0.80], dtype=float),
        label="top",
        leaf_count=500,
    )
    tree.add_node(
        "root",
        is_leaf=False,
        distribution=np.array([0.50, 0.50, 0.50, 0.50, 0.50, 0.50], dtype=float),
        label="root",
        leaf_count=300,
    )
    tree.add_node(
        "I",
        is_leaf=False,
        distribution=np.array([0.50, 0.50, 0.50, 0.50, 0.50, 0.50], dtype=float),
        label="I",
        leaf_count=200,
    )
    tree.add_node(
        "L1",
        is_leaf=True,
        distribution=np.array([0.10, 0.20, 0.20, 0.80, 0.90, 0.80], dtype=float),
        label="L1",
        leaf_count=100,
    )
    tree.add_node(
        "L2",
        is_leaf=True,
        distribution=np.array([0.90, 0.80, 0.80, 0.20, 0.10, 0.20], dtype=float),
        label="L2",
        leaf_count=100,
    )
    tree.add_node(
        "L3",
        is_leaf=True,
        distribution=np.array([0.80, 0.80, 0.80, 0.20, 0.20, 0.20], dtype=float),
        label="L3",
        leaf_count=100,
    )
    tree.add_node(
        "cal",
        is_leaf=False,
        distribution=np.array([0.50, 0.50, 0.50, 0.50, 0.50, 0.50], dtype=float),
        label="cal",
        leaf_count=200,
    )
    tree.add_node(
        "C",
        is_leaf=True,
        distribution=np.array([0.45, 0.45, 0.45, 0.55, 0.55, 0.55], dtype=float),
        label="C",
        leaf_count=100,
    )
    tree.add_node(
        "D",
        is_leaf=True,
        distribution=np.array([0.55, 0.55, 0.55, 0.45, 0.45, 0.45], dtype=float),
        label="D",
        leaf_count=100,
    )
    tree.add_edge("top", "root")
    tree.add_edge("top", "cal")
    tree.add_edge("root", "I")
    tree.add_edge("root", "L3")
    tree.add_edge("I", "L1")
    tree.add_edge("I", "L2")
    tree.add_edge("cal", "C")
    tree.add_edge("cal", "D")
    tree.graph["root"] = "top"

    annotations_df = pd.DataFrame(
        {
            "leaf_count": {
                "top": 500,
                "root": 300,
                "I": 200,
                "L1": 100,
                "L2": 100,
                "L3": 100,
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
            [1, 1, 1, 1, 0, 0],
            [0.45, 0.45, 0.45, 0.55, 0.55, 0.55],
            [0.55, 0.55, 0.55, 0.45, 0.45, 0.45],
        ],
        index=["L1", "L2", "L3", "C", "D"],
        dtype=float,
    )
    return tree, annotations_df, leaf_data


def test_cherry_with_leaf_data_uses_parent_dimension_for_leaf_pair_parent() -> None:
    tree, annotations_df, leaf_data = _build_cherry_tree()

    _edge_df, spectral_context = annotate_child_parent_divergence(
        tree,
        annotations_df.copy(),
        leaf_data=leaf_data,
    )
    test_projection_dimensions_by_node = spectral_context.test_projection_dimensions_by_node
    assert test_projection_dimensions_by_node["A"] == 0
    assert test_projection_dimensions_by_node["B"] == 0
    assert test_projection_dimensions_by_node["C"] == 0
    assert test_projection_dimensions_by_node["D"] == 0
    assert test_projection_dimensions_by_node["root"] > 0
    assert test_projection_dimensions_by_node["cal"] > 0
    assert test_projection_dimensions_by_node["top"] > 0

    sibling_projection_dimensions_from_edge_comparisons = (
        derive_sibling_projection_dimensions_from_child_edge_comparisons(
            tree,
            spectral_context=spectral_context,
        )
    )
    assert sibling_projection_dimensions_from_edge_comparisons == {
        "root": test_projection_dimensions_by_node["root"],
        "cal": test_projection_dimensions_by_node["cal"],
        "top": min(
            test_projection_dimensions_by_node["top"],
            round(
                np.sqrt(
                    test_projection_dimensions_by_node["root"]
                    * test_projection_dimensions_by_node["cal"]
                )
            ),
        ),
    }

    bundle = run_gate_annotation_pipeline(tree, annotations_df.copy(), leaf_data=leaf_data)
    assert (
        bundle.edge_gate_result.spectral_context.test_projection_dimensions_by_node
        == test_projection_dimensions_by_node
    )
    # One data row per child leaves no rows for the edge basis, so the leaf-pair parent stays closed.
    assert bundle.annotated_df.loc[["A", "B"], "Child_Parent_Divergence_df"].eq(0).all()
    assert bundle.annotated_df.loc["root", "Sibling_Divergence_Skipped"]


def test_mixed_parent_with_leaf_data_keeps_internal_parent_in_edge_derived_sibling_projection_dimensions() -> (
    None
):
    tree, annotations_df, leaf_data = _build_mixed_tree()

    bundle = run_gate_annotation_pipeline(tree, annotations_df.copy(), leaf_data=leaf_data)

    test_projection_dimensions_by_node = (
        bundle.edge_gate_result.spectral_context.test_projection_dimensions_by_node
    )
    assert test_projection_dimensions_by_node["L1"] == 0
    assert test_projection_dimensions_by_node["L2"] == 0
    assert test_projection_dimensions_by_node["L3"] == 0
    assert test_projection_dimensions_by_node["C"] == 0
    assert test_projection_dimensions_by_node["D"] == 0
    assert test_projection_dimensions_by_node["I"] > 0
    assert test_projection_dimensions_by_node["root"] > 0
    assert test_projection_dimensions_by_node["cal"] > 0
    assert test_projection_dimensions_by_node["top"] > 0

    sibling_projection_dimensions_from_edge_comparisons = (
        derive_sibling_projection_dimensions_from_child_edge_comparisons(
            tree,
            spectral_context=bundle.edge_gate_result.spectral_context,
        )
    )
    assert set(sibling_projection_dimensions_from_edge_comparisons) == {
        "root",
        "I",
        "cal",
        "top",
    }
    assert sibling_projection_dimensions_from_edge_comparisons["root"] > 0
    assert (
        sibling_projection_dimensions_from_edge_comparisons["I"]
        == (test_projection_dimensions_by_node["I"])
    )
    assert (
        sibling_projection_dimensions_from_edge_comparisons["cal"]
        == (test_projection_dimensions_by_node["cal"])
    )
    assert (
        0
        < sibling_projection_dimensions_from_edge_comparisons["top"]
        <= (test_projection_dimensions_by_node["top"])
    )


def test_decompose_without_leaf_data_raises_for_missing_spectral_contract() -> None:
    tree, annotations_df, _leaf_data = _build_cherry_tree()

    with pytest.raises(ValueError, match="require leaf_data"):
        run_gate_annotation_pipeline(
            tree,
            annotations_df.copy(),
            leaf_data=None,
        )


def test_edge_derived_sibling_dimensions_require_all_child_spectral_dimensions() -> None:
    tree, _annotations_df, _leaf_data = _build_cherry_tree()
    spectral_context = SpectralDecompositionResult(
        test_projection_dimensions_by_node={"A": 1},
        raw_mp_signal_counts_by_node={"A": 0},
        effective_independent_rows_by_node={"A": 1},
        mp_threshold_rows_by_node={"A": 1},
        principal_component_projections_by_node={},
        principal_component_eigenvalues_by_node={},
        descendant_leaf_row_counts_by_node={},
        internal_distribution_row_counts_by_node={},
        spectral_matrix_row_counts_by_node={},
        active_feature_counts_by_node={},
        stage_timings={},
    )

    with pytest.raises(ValueError, match="missing 'root'"):
        derive_sibling_projection_dimensions_from_child_edge_comparisons(
            tree,
            spectral_context=spectral_context,
        )


def test_only_the_projected_wald_sibling_test_uses_branch_length(
    expand_leaves_into_rows,
) -> None:
    def annotate(root_pair_scale: float, **gate_kwargs: str) -> tuple[pd.DataFrame, float]:
        tree, annotations_df, leaf_data = expand_leaves_into_rows(*_build_cherry_tree())
        # "top" differs from its children only in its stored distribution, which the row-based
        # edge test ignores, so start the tree at "root" to keep the root pair tested.
        tree.remove_node("top")
        tree.graph["root"] = "root"
        annotations_df = annotations_df.drop(index="top")
        for child in ("A", "B"):
            tree.edges["root", child]["branch_length"] *= root_pair_scale
        lengths = [length for *_, length in tree.edges(data="branch_length")]
        multiplier = 1.0 + (
            tree.edges["root", "A"]["branch_length"] + tree.edges["root", "B"]["branch_length"]
        ) / (2.0 * float(np.mean(lengths)))
        bundle = run_gate_annotation_pipeline(
            tree, annotations_df.copy(), leaf_data=leaf_data, **gate_kwargs
        )
        assert bundle.metadata.config.edge_branch_length_variance_policy == "none"
        return bundle.annotated_df, multiplier

    base_df, base_multiplier = annotate(1.0)
    long_df, long_multiplier = annotate(4.0)

    pd.testing.assert_series_equal(
        base_df["Child_Parent_Divergence_Test_Statistic"],
        long_df["Child_Parent_Divergence_Test_Statistic"],
    )
    assert (
        long_df.loc["root", "Sibling_Test_Statistic"]
        < base_df.loc["root", "Sibling_Test_Statistic"]
    )
    assert long_df.loc["root", "Sibling_Test_Statistic"] * long_multiplier == pytest.approx(
        base_df.loc["root", "Sibling_Test_Statistic"] * base_multiplier
    )

    fixed_base_df, _ = annotate(1.0, sibling_gate_method="fixed_coordinate_bh")
    fixed_long_df, _ = annotate(4.0, sibling_gate_method="fixed_coordinate_bh")
    for column in ("Sibling_Test_Statistic", "Sibling_Dense_Evidence_Test_Statistic"):
        assert fixed_long_df.loc["root", column] == pytest.approx(fixed_base_df.loc["root", column])
