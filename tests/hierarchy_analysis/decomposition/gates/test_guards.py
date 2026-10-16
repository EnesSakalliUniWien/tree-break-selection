"""Contracts for gate guard implementations."""

from __future__ import annotations

import networkx as nx
import numpy as np
import pandas as pd
import pytest

import tree_break_selection.hierarchy_analysis.decomposition.gates.guards as guard_module
from tree_break_selection.hierarchy_analysis.decomposition.gates.guards import (
    apply_root_selective_permutation_guard,
    apply_root_stability_guard,
    selected_global_sibling_min_permutation_p_value,
)
from tree_break_selection.hierarchy_analysis.decomposition.gates.orchestrator import (
    run_gate_annotation_pipeline,
)
from tree_break_selection.tree.feature_space import (
    FeatureSpace,
    infer_feature_space_from_columns,
)

from .gate_annotation_support import (
    _build_small_tree_with_leaf_data,
    _selected_root,
)


def _build_passthrough_tree_with_leaf_data() -> tuple[
    nx.DiGraph,
    pd.DataFrame,
    pd.DataFrame,
]:
    tree = nx.DiGraph()
    tree.graph["root"] = "root"
    tree.add_edges_from(
        [
            ("root", "A"),
            ("root", "B"),
            ("B", "C"),
            ("B", "D"),
        ]
    )
    for node in ("root", "A", "B", "C", "D"):
        tree.nodes[node]["is_leaf"] = node in {"A", "C", "D"}
        tree.nodes[node]["label"] = node
        tree.nodes[node]["leaf_count"] = 1 if tree.nodes[node]["is_leaf"] else 3
        tree.nodes[node]["distribution"] = np.array(
            [0.5, 0.5, 0.5, 0.5],
            dtype=np.float64,
        )

    annotations_df = pd.DataFrame(index=["root", "A", "B", "C", "D"])
    leaf_data = pd.DataFrame(
        [
            [0, 0, 1, 1],
            [1, 1, 0, 0],
            [0, 1, 0, 1],
        ],
        index=["A", "C", "D"],
        dtype=np.float64,
    )
    return tree, annotations_df, leaf_data


def _passthrough_guard_inputs(
    *,
    root_open: bool = False,
) -> tuple[nx.DiGraph, pd.DataFrame, pd.DataFrame, FeatureSpace]:
    tree, annotations_df, leaf_data = _build_passthrough_tree_with_leaf_data()
    feature_space = infer_feature_space_from_columns(tuple(leaf_data.columns))
    assert feature_space is not None
    annotated = annotations_df.copy()
    annotated["Child_Parent_Divergence_Significant"] = True
    annotated["Sibling_BH_Different"] = False
    annotated["Sibling_BH_Same"] = True
    annotated["Sibling_Divergence_Skipped"] = False
    annotated["Sibling_Divergence_P_Value"] = 1.0
    for node in ("root", "B") if root_open else ("B",):
        annotated.loc[node, "Sibling_BH_Different"] = True
        annotated.loc[node, "Sibling_BH_Same"] = False
        annotated.loc[node, "Sibling_Divergence_P_Value"] = 0.001
    return tree, annotated, leaf_data, feature_space


def _recording_selected_root(calls: list[tuple[object, ...]], selected_p_value: float):
    def runner(data, *_args, **_kwargs):
        calls.append(tuple(data.index))
        return {
            "root_observed_p_value": 0.001,
            "root_selective_p_value": selected_p_value,
            "root_selective_null_min_p_value": 0.01,
            "root_selective_null_q05_p_value": 0.05,
        }

    return runner


def test_root_stability_guard_closes_only_unstable_open_root() -> None:
    tree, annotations_df, leaf_data = _build_small_tree_with_leaf_data()
    bundle = run_gate_annotation_pipeline(
        tree,
        annotations_df.copy(),
        edge_alpha=0.01,
        sibling_alpha=0.01,
        leaf_data=leaf_data,
        sibling_gate_method="fixed_coordinate_bh",
    )
    root = "root"
    annotated = bundle.annotated_df.copy()
    annotated.loc[root, "Sibling_BH_Different"] = True
    annotated.loc[root, "Sibling_BH_Same"] = False
    annotated.loc["A", "Sibling_BH_Different"] = True

    guarded = apply_root_stability_guard(
        tree,
        annotated,
        {
            "root_stability_subsample_mean_ari": 0.20,
            "root_stability_subsample_median_ari": 0.25,
            "root_stability_subsample_q10_ari": 0.10,
        },
        threshold=0.50,
    )

    assert bool(guarded.loc[root, "Sibling_BH_Different"]) is False
    assert bool(guarded.loc[root, "Sibling_BH_Same"]) is True
    assert bool(guarded.loc[root, "Root_Stability_Guard_Blocked"]) is True
    assert bool(guarded.loc["A", "Sibling_BH_Different"]) is True
    assert guarded.loc[root, "Root_Stability_Subsample_Mean_ARI"] == 0.20


def test_root_selective_permutation_guard_closes_unselected_open_root(
    monkeypatch,
) -> None:
    tree, annotations_df, leaf_data = _build_small_tree_with_leaf_data()
    feature_space = infer_feature_space_from_columns(tuple(leaf_data.columns))
    bundle = run_gate_annotation_pipeline(
        tree,
        annotations_df.copy(),
        edge_alpha=0.01,
        sibling_alpha=0.01,
        leaf_data=leaf_data,
        feature_space=feature_space,
        sibling_gate_method="fixed_coordinate_bh",
    )
    root = "root"
    annotated = bundle.annotated_df.copy()
    annotated.loc[root, "Sibling_BH_Different"] = True
    annotated.loc[root, "Sibling_BH_Same"] = False
    annotated.loc["A", "Sibling_BH_Different"] = True

    fake_selected_root = _selected_root(0.50)

    monkeypatch.setattr(
        guard_module,
        "selected_root_permutation_p_value",
        fake_selected_root,
    )

    guarded = apply_root_selective_permutation_guard(
        tree,
        annotated,
        leaf_data,
        feature_space,
        method="fixed_coordinate_bh",
        bootstrap_replicates=5,
        seed=123,
        alpha=0.01,
    )

    assert bool(guarded.loc[root, "Sibling_BH_Different"]) is False
    assert bool(guarded.loc[root, "Sibling_BH_Same"]) is True
    assert bool(guarded.loc[root, "Root_Selective_Permutation_Guard_Blocked"]) is True
    assert bool(guarded.loc["A", "Sibling_BH_Different"]) is True
    assert guarded.loc[root, "Root_Selective_Permutation_P_Value"] == 0.50
    assert guarded.loc[root, "Root_Selective_Permutation_Guard_Alpha"] == 0.01


def test_selective_permutation_guard_closes_open_internal_context(
    monkeypatch,
) -> None:
    tree, annotations_df, leaf_data = _build_small_tree_with_leaf_data()
    feature_space = infer_feature_space_from_columns(tuple(leaf_data.columns))
    bundle = run_gate_annotation_pipeline(
        tree,
        annotations_df.copy(),
        edge_alpha=0.01,
        sibling_alpha=0.01,
        leaf_data=leaf_data,
        feature_space=feature_space,
        sibling_gate_method="fixed_coordinate_bh",
    )
    annotated = bundle.annotated_df.copy()
    annotated.loc["root", "Sibling_BH_Different"] = False
    annotated.loc["root", "Sibling_BH_Same"] = True
    annotated.loc["cal", "Sibling_BH_Different"] = True
    annotated.loc["cal", "Sibling_BH_Same"] = False

    fake_selected_root = _selected_root(0.50)

    monkeypatch.setattr(
        guard_module,
        "selected_root_permutation_p_value",
        fake_selected_root,
    )

    guarded = apply_root_selective_permutation_guard(
        tree,
        annotated,
        leaf_data,
        feature_space,
        method="fixed_coordinate_bh",
        bootstrap_replicates=5,
        seed=123,
        alpha=0.01,
        scope="open_internal",
    )

    assert bool(guarded.loc["root", "Sibling_BH_Different"]) is False
    assert bool(guarded.loc["cal", "Sibling_BH_Different"]) is False
    assert bool(guarded.loc["cal", "Sibling_BH_Same"]) is True
    assert bool(guarded.loc["cal", "Selective_Permutation_Guard_Blocked"]) is True
    assert guarded.loc["cal", "Selective_Permutation_P_Value"] == 0.50
    assert guarded.loc["cal", "Selective_Permutation_Guard_Scope"] == "open_internal"


def test_selective_permutation_guard_closes_only_passthrough_descendant(
    monkeypatch,
) -> None:
    tree, annotated, leaf_data, feature_space = _passthrough_guard_inputs()

    calls: list[tuple[object, ...]] = []

    fake_selected_root = _recording_selected_root(calls, 0.50)

    monkeypatch.setattr(
        guard_module,
        "selected_root_permutation_p_value",
        fake_selected_root,
    )

    guarded = apply_root_selective_permutation_guard(
        tree,
        annotated,
        leaf_data,
        feature_space,
        method="fixed_coordinate_bh",
        bootstrap_replicates=5,
        seed=123,
        alpha=0.01,
        scope="passthrough_descendant",
    )

    assert calls == [("C", "D")]
    assert bool(guarded.loc["root", "Sibling_BH_Different"]) is False
    assert bool(guarded.loc["B", "Sibling_BH_Different"]) is False
    assert bool(guarded.loc["B", "Selective_Permutation_Guard_Blocked"]) is True
    assert guarded.loc["B", "Selective_Permutation_Guard_Scope"] == ("passthrough_descendant")


def test_global_selected_family_guard_closes_passthrough_descendant(
    monkeypatch,
) -> None:
    tree, annotated, leaf_data, feature_space = _passthrough_guard_inputs()

    seen: dict[str, object] = {}

    def fake_selected_family(data, *_args, **kwargs):
        seen["index"] = tuple(data.index)
        seen["observed"] = kwargs["observed_p_value"]
        return {
            "root_observed_p_value": kwargs["observed_p_value"],
            "root_selective_p_value": 0.50,
            "root_selective_null_min_p_value": 0.0001,
            "root_selective_null_q05_p_value": 0.005,
        }

    monkeypatch.setattr(
        guard_module,
        "selected_global_sibling_min_permutation_p_value",
        fake_selected_family,
    )

    guarded = apply_root_selective_permutation_guard(
        tree,
        annotated,
        leaf_data,
        feature_space,
        method="fixed_coordinate_bh",
        bootstrap_replicates=5,
        seed=123,
        alpha=0.01,
        scope="global_sibling_min_passthrough_descendant",
    )

    assert seen == {"index": ("A", "C", "D"), "observed": 0.001}
    assert bool(guarded.loc["root", "Sibling_BH_Different"]) is False
    assert bool(guarded.loc["B", "Sibling_BH_Different"]) is False
    assert bool(guarded.loc["B", "Selective_Permutation_Guard_Blocked"]) is True
    assert guarded.loc["B", "Selective_Permutation_P_Value"] == 0.50
    assert guarded.loc["B", "Selective_Permutation_Guard_Scope"] == (
        "global_sibling_min_passthrough_descendant"
    )


def test_global_selected_family_guard_keeps_significant_passthrough_descendant(
    monkeypatch,
) -> None:
    tree, annotated, leaf_data, feature_space = _passthrough_guard_inputs()

    def fake_selected_family(*_args, **kwargs):
        return {
            "root_observed_p_value": kwargs["observed_p_value"],
            "root_selective_p_value": 0.01,
            "root_selective_null_min_p_value": 0.0001,
            "root_selective_null_q05_p_value": 0.005,
        }

    monkeypatch.setattr(
        guard_module,
        "selected_global_sibling_min_permutation_p_value",
        fake_selected_family,
    )

    guarded = apply_root_selective_permutation_guard(
        tree,
        annotated,
        leaf_data,
        feature_space,
        method="fixed_coordinate_bh",
        bootstrap_replicates=5,
        seed=123,
        alpha=0.01,
        scope="global_sibling_min_passthrough_descendant",
    )

    assert bool(guarded.loc["B", "Sibling_BH_Different"]) is True
    assert bool(guarded.loc["B", "Selective_Permutation_Guard_Blocked"]) is False
    assert bool(guarded.loc["B", "Selective_Permutation_Guard_Would_Block"]) is False


def test_refined_global_selected_family_guard_refines_floor_p_value(
    monkeypatch,
) -> None:
    tree, annotated, leaf_data, feature_space = _passthrough_guard_inputs()

    calls: list[int] = []

    def fake_selected_family(*_args, **kwargs):
        replicates = int(kwargs["bootstrap_replicates"])
        calls.append(replicates)
        selected_p = 0.01 if replicates == 99 else 0.02
        return {
            "root_observed_p_value": kwargs["observed_p_value"],
            "root_selective_p_value": selected_p,
            "root_selective_null_min_p_value": 0.0001,
            "root_selective_null_q05_p_value": 0.005,
        }

    monkeypatch.setattr(
        guard_module,
        "selected_global_sibling_min_permutation_p_value",
        fake_selected_family,
    )

    guarded = apply_root_selective_permutation_guard(
        tree,
        annotated,
        leaf_data,
        feature_space,
        method="fixed_coordinate_bh",
        bootstrap_replicates=99,
        seed=123,
        alpha=0.01,
        scope="global_sibling_min_passthrough_descendant_refined",
    )

    assert calls == [99, 999]
    assert bool(guarded.loc["B", "Sibling_BH_Different"]) is False
    assert bool(guarded.loc["B", "Selective_Permutation_Guard_Blocked"]) is True
    assert bool(guarded.loc["B", "Selective_Permutation_Guard_Refined"]) is True
    assert guarded.loc["B", "Selective_Permutation_Base_P_Value"] == 0.01
    assert guarded.loc["B", "Selective_Permutation_P_Value"] == 0.02
    assert guarded.loc["B", "Selective_Permutation_Guard_Replicates"] == 999


def test_block_permutation_null_sample_keeps_fractional_values() -> None:
    leaf_data = pd.DataFrame(
        [[0.7, 0.2, 0.8], [0.9, 0.5, 0.5], [0.2, 1.0, 0.0]],
        index=["a", "b", "c"],
        columns=["F0", "F1_c0", "F1_c1"],
    )
    feature_space = infer_feature_space_from_columns(tuple(leaf_data.columns))
    assert [block.family for block in feature_space.blocks] == ["bernoulli", "categorical"]

    null_sample = guard_module._block_permutation_null_sample(
        leaf_data,
        feature_space,
        np.random.default_rng(0),
    ).to_numpy()

    source = leaf_data.to_numpy()
    assert sorted(null_sample[:, 0]) == sorted(source[:, 0])
    assert sorted(map(tuple, null_sample[:, 1:])) == sorted(map(tuple, source[:, 1:]))


def test_global_selected_family_p_value_uses_add_one_monte_carlo(
    monkeypatch,
) -> None:
    leaf_data = pd.DataFrame([[0], [1], [0]], index=["a", "b", "c"])
    feature_space = infer_feature_space_from_columns(tuple(leaf_data.columns))
    null_values = iter([0.01, 0.20, 0.05])

    def fake_null_sample(data, *_args, **_kwargs):
        return data

    def fake_min_p_value(*_args, **_kwargs):
        return next(null_values)

    monkeypatch.setattr(
        guard_module,
        "_block_permutation_null_sample",
        fake_null_sample,
    )
    monkeypatch.setattr(
        guard_module,
        "_selected_tree_fixed_sibling_min_p_value",
        fake_min_p_value,
    )

    result = selected_global_sibling_min_permutation_p_value(
        leaf_data,
        feature_space,
        method="fixed_coordinate_bh",
        bootstrap_replicates=3,
        seed=123,
        observed_p_value=0.05,
    )

    assert result["root_observed_p_value"] == 0.05
    assert result["root_selective_p_value"] == 0.75
    assert result["root_selective_null_min_p_value"] == 0.01
    assert result["root_selective_null_q05_p_value"] == pytest.approx(0.014)


@pytest.mark.parametrize(
    ("selected_p_value", "expected_root_open"),
    [
        pytest.param(0.01, True, id="open-root-keeps-descendant"),
        pytest.param(0.50, False, id="blocked-root-skips-descendant"),
    ],
)
def test_selective_permutation_passthrough_scope_respects_root_guard_outcome(
    monkeypatch,
    selected_p_value: float,
    expected_root_open: bool,
) -> None:
    tree, annotated, leaf_data, feature_space = _passthrough_guard_inputs(root_open=True)

    calls: list[tuple[object, ...]] = []
    monkeypatch.setattr(
        guard_module,
        "selected_root_permutation_p_value",
        _recording_selected_root(calls, selected_p_value),
    )

    guarded = apply_root_selective_permutation_guard(
        tree,
        annotated,
        leaf_data,
        feature_space,
        method="fixed_coordinate_bh",
        bootstrap_replicates=5,
        seed=123,
        alpha=0.01,
        scope="passthrough_descendant",
    )

    assert calls == [("A", "C", "D")]
    assert bool(guarded.loc["root", "Sibling_BH_Different"]) is expected_root_open
    assert bool(guarded.loc["root", "Root_Selective_Permutation_Guard_Blocked"]) is (
        not expected_root_open
    )
    assert bool(guarded.loc["B", "Sibling_BH_Different"]) is True
    assert guarded.loc["B", "Selective_Permutation_Guard_Scope"] == ""


def test_non_root_guard_scope_requires_replicates() -> None:
    with pytest.raises(ValueError, match="requires root_selective_permutation_guard_replicates"):
        guard_module.validate_root_selective_permutation_guard_config(
            replicates=0,
            alpha=0.05,
            scope="passthrough_descendant",
        )


def test_passthrough_guard_candidates_skip_nodes_below_a_boundary() -> None:
    tree = nx.DiGraph()
    tree.add_edges_from(
        [
            ("R", "A"),
            ("R", "B"),
            ("A", "A1"),
            ("A", "A2"),
            ("B", "C"),
            ("B", "b1"),
            ("C", "D"),
            ("C", "c1"),
            ("D", "d1"),
            ("D", "d2"),
        ]
    )
    edge_divergent = {node: True for node in tree.nodes}
    edge_divergent["C"] = False
    edge_divergent["b1"] = False
    sibling_different = {node: node in {"A", "D"} for node in tree.nodes}
    annotations = pd.DataFrame(
        {
            "Child_Parent_Divergence_Significant": edge_divergent,
            "Sibling_BH_Different": sibling_different,
            "Sibling_Divergence_Skipped": {node: False for node in tree.nodes},
        }
    )

    # B's edges are closed, so traversal never reaches D.
    assert guard_module._passthrough_descendant_guard_candidates(tree, annotations, root="R") == [
        "A"
    ]
