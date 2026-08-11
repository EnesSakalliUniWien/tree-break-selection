"""
Production integration tests for fail-closed calibration on complex synthetic data.

Tests the full pipeline with:
- Balanced binary feature matrices with low entropy
- Unbalanced binary feature matrices with high entropy
"""

import pandas as pd
from benchmarks.shared.generators import generate_random_feature_matrix
from benchmarks.shared.util.decomposition import _labels_from_decomposition
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist
from tree_break_selection.hierarchy_analysis.decomposition.gates.orchestrator import (
    run_gate_annotation_pipeline,
)
from tree_break_selection.hierarchy_analysis.statistics.alpha_contract import (
    DEFAULT_EDGE_ALPHA,
)
from tree_break_selection.hierarchy_analysis.tree_decomposition import TreeDecomposition
from tree_break_selection.tree.construction import tree_from_linkage


def _run_pipeline_on_dataframe(data_df, significance_level=0.05, **kwargs):
    """Minimal pipeline helper for integration tests."""
    Z = linkage(pdist(data_df.values, metric="hamming"), method="complete")
    tree = tree_from_linkage(Z, leaf_names=data_df.index.tolist())
    tree.populate_node_divergences(data_df)
    gate_bundle = run_gate_annotation_pipeline(
        tree,
        tree.annotations_df.copy(),
        leaf_data=data_df,
        edge_alpha=DEFAULT_EDGE_ALPHA,
        sibling_alpha=significance_level,
        **kwargs,
    )
    decomposition = TreeDecomposition(
        tree=tree,
        gate_annotation_bundle=gate_bundle,
    ).decompose_tree()
    return decomposition, tree, gate_bundle


def _assert_selected_hierarchy_calibration_fails_closed(
    decomposition: dict,
    annotations: pd.DataFrame,
    sample_names: list[str],
) -> None:
    predicted = _labels_from_decomposition(decomposition, sample_names)
    assert all(label != -1 for label in predicted)
    assert decomposition["num_clusters"] == 1

    unresolved = annotations[
        annotations["Sibling_Gate_P_Value_Calibration"].eq(
            "undefined_unvalidated_reference_law"
        )
    ]
    assert not unresolved.empty
    assert unresolved["Sibling_Gate_P_Value_Role"].eq(
        "fail_closed_sibling_gate"
    ).all()
    assert unresolved["Sibling_Divergence_P_Value"].between(0.0, 1.0).all()
    assert unresolved["Sibling_Divergence_Skipped"].eq(True).all()
    assert unresolved["Sibling_Divergence_Invalid"].eq(True).all()
    assert not unresolved["Sibling_BH_Different"].any()


def test_complex_random_feature_matrix_balanced_clusters():
    """Low-entropy selected-hierarchy calibration must fail closed."""
    data_dict, _true_clusters = generate_random_feature_matrix(
        n_rows=72,
        n_cols=40,
        entropy_param=0.1,
        n_clusters=4,
        random_seed=314,
        balanced_clusters=True,
    )
    data_df = pd.DataFrame.from_dict(data_dict, orient="index").astype(int)

    decomposition, _, gate_bundle = _run_pipeline_on_dataframe(
        data_df,
        significance_level=0.05,
    )
    _assert_selected_hierarchy_calibration_fails_closed(
        decomposition,
        gate_bundle.annotated_df,
        data_df.index.tolist(),
    )


def test_complex_random_feature_matrix_unbalanced_clusters():
    """Unbalanced selected-hierarchy calibration must also fail closed."""
    data_dict, _true_clusters = generate_random_feature_matrix(
        n_rows=150,
        n_cols=36,
        entropy_param=0.25,
        n_clusters=4,
        random_seed=2024,
        balanced_clusters=False,
    )
    data_df = pd.DataFrame.from_dict(data_dict, orient="index").astype(int)

    decomposition, _, gate_bundle = _run_pipeline_on_dataframe(
        data_df,
        significance_level=0.05,
    )
    _assert_selected_hierarchy_calibration_fails_closed(
        decomposition,
        gate_bundle.annotated_df,
        data_df.index.tolist(),
    )
