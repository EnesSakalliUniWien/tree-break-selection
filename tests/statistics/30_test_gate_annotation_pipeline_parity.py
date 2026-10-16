from __future__ import annotations

import networkx as nx
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist

from benchmarks.shared.cases import get_default_test_cases
from benchmarks.shared.generators import generate_case_data
from tree_break_selection.hierarchy_analysis.decomposition.gates.orchestrator import (
    run_gate_annotation_pipeline,
)
from tree_break_selection.hierarchy_analysis.statistics.child_parent_divergence import (
    annotate_child_parent_divergence,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.inflated_projected_wald_annotation.pipeline import (
    annotate_sibling_divergence,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.projection.gate_inputs.parent_principal_component_inputs import (
    collect_parent_principal_component_inputs_for_sibling_tests,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.projection.gate_inputs.projection_dimensions import (
    derive_sibling_projection_dimensions_from_child_edge_comparisons,
)
from tree_break_selection.tree.construction.build import build_tree


def _build_small_binary_tree() -> tuple[nx.DiGraph, pd.DataFrame, pd.DataFrame]:
    tree = nx.DiGraph()
    tree.add_edge("root", "A", branch_length=0.25)
    tree.add_edge("root", "B", branch_length=0.20)
    tree.add_edge("cal", "C", branch_length=0.10)
    tree.add_edge("cal", "D", branch_length=0.10)

    root_dist = np.array([0.50, 0.50, 0.50, 0.50, 0.50, 0.50], dtype=np.float64)
    a_dist = np.array([0.12, 0.12, 0.12, 0.88, 0.88, 0.88], dtype=np.float64)
    b_dist = np.array([0.88, 0.88, 0.88, 0.12, 0.12, 0.12], dtype=np.float64)
    cal_dist = np.array([0.50, 0.50, 0.50, 0.50, 0.50, 0.50], dtype=np.float64)
    c_dist = np.array([0.45, 0.45, 0.45, 0.55, 0.55, 0.55], dtype=np.float64)
    d_dist = np.array([0.55, 0.55, 0.55, 0.45, 0.45, 0.45], dtype=np.float64)

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

    base_df = pd.DataFrame(
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
        [a_dist, b_dist, c_dist, d_dist],
        index=["A", "B", "C", "D"],
        columns=[f"feature_{feature_index}" for feature_index in range(len(a_dist))],
    )
    return tree, base_df, leaf_data


def test_gate_annotation_pipeline_matches_sequential_gate_annotations() -> None:
    tree, base_df, leaf_data = _build_small_binary_tree()

    edge_df, spectral_context = annotate_child_parent_divergence(
        tree,
        base_df.copy(),
        significance_level_alpha=0.01,
        leaf_data=leaf_data,
    )
    sibling_projection_dimensions = (
        derive_sibling_projection_dimensions_from_child_edge_comparisons(
            tree,
            spectral_context=spectral_context,
        )
    )
    (
        parent_principal_component_projections,
        parent_principal_component_eigenvalues,
    ) = collect_parent_principal_component_inputs_for_sibling_tests(
        sibling_projection_dimensions,
        spectral_context=spectral_context,
    )
    sequential_df = annotate_sibling_divergence(
        tree,
        edge_df,
        significance_level_alpha=0.01,
        sibling_projection_dimensions_from_edge_comparisons=sibling_projection_dimensions,
        parent_principal_component_projections=parent_principal_component_projections,
        parent_principal_component_eigenvalues=parent_principal_component_eigenvalues,
    )
    bundled = run_gate_annotation_pipeline(
        tree,
        base_df.copy(),
        edge_alpha=0.01,
        sibling_alpha=0.01,
        leaf_data=leaf_data,
    )

    fail_closed_nodes = ["root", "cal"]
    assert (
        sequential_df.loc[fail_closed_nodes, "Sibling_Test_Method"]
        .eq("empirical_null_no_internal_support")
        .all()
    )
    assert (
        sequential_df.loc[fail_closed_nodes, "Sibling_Gate_P_Value_Calibration"]
        .eq("undefined_no_internal_support")
        .all()
    )
    assert (
        sequential_df.loc[fail_closed_nodes, "Sibling_Gate_P_Value_Role"]
        .eq("fail_closed_sibling_gate")
        .all()
    )
    assert sequential_df.loc[fail_closed_nodes, "Sibling_Divergence_Skipped"].eq(True).all()
    assert sequential_df.loc[fail_closed_nodes, "Sibling_Divergence_Invalid"].eq(True).all()
    assert sequential_df.loc[fail_closed_nodes, "Sibling_BH_Different"].eq(False).all()
    _assert_matches_outside_evidence_channels(sequential_df, bundled.annotated_df)


def _sequential_gate_annotations(tree, base_df, leaf_data) -> pd.DataFrame:
    edge_df, spectral_context = annotate_child_parent_divergence(
        tree,
        base_df.copy(),
        significance_level_alpha=0.01,
        leaf_data=leaf_data,
    )
    sibling_projection_dimensions = (
        derive_sibling_projection_dimensions_from_child_edge_comparisons(
            tree,
            spectral_context=spectral_context,
        )
    )
    projections, eigenvalues = collect_parent_principal_component_inputs_for_sibling_tests(
        sibling_projection_dimensions,
        spectral_context=spectral_context,
    )
    return annotate_sibling_divergence(
        tree,
        edge_df,
        significance_level_alpha=0.01,
        sibling_projection_dimensions_from_edge_comparisons=sibling_projection_dimensions,
        parent_principal_component_projections=projections,
        parent_principal_component_eigenvalues=eigenvalues,
    )


def _assert_matches_outside_evidence_channels(
    sequential_df: pd.DataFrame, bundled_df: pd.DataFrame
) -> None:
    # Only the bundled run adds the fixed-subspace evidence channels; every other column
    # must match, so a wrong alpha or dimension forwarded by the pipeline fails here.
    evidence_prefixes = ("Sibling_Fixed_", "Sibling_Sparse_Evidence_", "Sibling_Dense_Evidence_")
    compared = [
        column for column in sequential_df.columns if not column.startswith(evidence_prefixes)
    ]
    pd.testing.assert_frame_equal(
        sequential_df[compared],
        bundled_df.loc[sequential_df.index, compared],
    )


def test_gate_annotation_pipeline_matches_sequential_gates_on_benchmark_case() -> None:
    # Unlike the fail-closed fixture above, this case has edge decisions that move with alpha.
    case = next(c for c in get_default_test_cases() if c["name"] == "binary_balanced_low_noise")
    data, _labels, _original, _metadata = generate_case_data(case)

    def fresh_tree():
        tree = build_tree(
            data, pdist(data.to_numpy(float), "hamming"), linkage_method="average"
        ).tree
        tree.populate_node_divergences(data)
        return tree

    sequential_tree = fresh_tree()
    sequential_df = _sequential_gate_annotations(
        sequential_tree, sequential_tree.annotations_df, data
    )
    bundled_tree = fresh_tree()
    bundled = run_gate_annotation_pipeline(
        bundled_tree,
        bundled_tree.annotations_df.copy(),
        edge_alpha=0.01,
        sibling_alpha=0.01,
        leaf_data=data,
    )

    assert sequential_df["Child_Parent_Divergence_Significant"].any()
    _assert_matches_outside_evidence_channels(sequential_df, bundled.annotated_df)
