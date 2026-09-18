"""Contracts for the gate annotation orchestrator."""

from __future__ import annotations

import json
from dataclasses import FrozenInstanceError, asdict

import numpy as np
import pytest
import tree_break_selection.hierarchy_analysis.decomposition.gates.guards as guard_module
from tree_break_selection.hierarchy_analysis.decomposition.gates.column_contracts import (
    EDGE_GATE_COLUMNS,
    SIBLING_GATE_COLUMNS,
)
from tree_break_selection.hierarchy_analysis.decomposition.gates.orchestrator import (
    build_gate_annotation_config_metadata,
    resolve_effective_sibling_alpha,
    run_gate_annotation_pipeline,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.inflation_correction.types.inflation_model import (
    CalibrationSupportThresholds,
)
from tree_break_selection.tree.feature_space import infer_feature_space_from_columns

from .gate_annotation_support import (
    _build_small_tree_with_leaf_data,
    _selected_root,
)


def _stable_root(*_args, **_kwargs):
    return {
        "root_stability_subsample_mean_ari": 1.0,
        "root_stability_subsample_median_ari": 1.0,
        "root_stability_subsample_q10_ari": 1.0,
    }


def test_pipeline_supports_current_gate_annotation_contract() -> None:
    tree, annotations_df, leaf_data = _build_small_tree_with_leaf_data()

    bundle = run_gate_annotation_pipeline(
        tree,
        annotations_df.copy(),
        edge_alpha=0.01,
        sibling_alpha=0.01,
        leaf_data=leaf_data,
    )
    out = bundle.annotated_df
    for col in EDGE_GATE_COLUMNS:
        assert col in out.columns
    for col in SIBLING_GATE_COLUMNS:
        assert col in out.columns
    assert bundle.metadata.pipeline == "gate_annotation"
    assert bundle.metadata.edge.gate == "edge"
    assert bundle.metadata.sibling.gate == "sibling"
    assert bundle.metadata.config.sibling_gate_method == "projected_wald_inflation"
    assert bundle.metadata.config.sibling_gate_alpha_penalty == 1.0
    tested = out[out["Sibling_Sparse_Evidence_Method"].eq("fixed_coordinate_bh")]
    assert not tested.empty
    assert tested["Sibling_Sparse_Evidence_P_Value"].between(0.0, 1.0).all()
    active_tested = out[out["Sibling_Gate_P_Value_Role"].eq("active_traversal_sibling_gate")]
    assert not active_tested.empty
    assert active_tested["Sibling_Gate_P_Value_Calibration"].eq("empirical_null_inflation").all()
    assert active_tested["Sibling_Divergence_P_Value"].between(0.0, 1.0).all()
    assert tested["Sibling_Dense_Evidence_Method"].eq("fixed_global_chi_square").all()
    assert tested["Sibling_Dense_Evidence_Calibration"].eq("fixed_subspace_chi_square").all()


def test_gate_metadata_owns_typed_versioned_support_policy_snapshot() -> None:
    thresholds = CalibrationSupportThresholds(
        min_supported_groups=3,
        min_family_supported_groups=2,
        min_family_effective_sample_size=1.5,
        min_local_effective_sample_size=1.25,
        max_weight_share=0.75,
        max_leave_one_group_delta_log_c=0.5,
    )

    metadata = build_gate_annotation_config_metadata(
        internal_support_thresholds=thresholds,
    )
    policy = metadata.internal_support_policy

    assert policy.support_contract_version == 2
    assert policy.thresholds == thresholds
    with pytest.raises(FrozenInstanceError):
        policy.support_contract_version = 3  # type: ignore[misc]

    payload = json.loads(json.dumps(asdict(metadata)))
    assert payload["internal_support_policy"] == {
        "support_contract_version": 2,
        "thresholds": {
            "min_supported_groups": 3,
            "min_family_supported_groups": 2,
            "min_family_effective_sample_size": 1.5,
            "min_local_effective_sample_size": 1.25,
            "max_weight_share": 0.75,
            "max_leave_one_group_delta_log_c": 0.5,
        },
    }


def test_pipeline_supports_opt_in_fixed_coordinate_sibling_gate() -> None:
    tree, annotations_df, leaf_data = _build_small_tree_with_leaf_data()

    bundle = run_gate_annotation_pipeline(
        tree,
        annotations_df.copy(),
        edge_alpha=0.01,
        sibling_alpha=0.01,
        leaf_data=leaf_data,
        sibling_gate_method="fixed_coordinate_bh",
    )
    out = bundle.annotated_df

    for col in EDGE_GATE_COLUMNS:
        assert col in out.columns
    for col in SIBLING_GATE_COLUMNS:
        assert col in out.columns
    assert bundle.metadata.sibling.gate == "sibling"
    assert bundle.metadata.config.sibling_gate_method == "fixed_coordinate_bh"
    assert bundle.metadata.config.sibling_gate_alpha_penalty == 1.0
    tested = out[out["Sibling_Test_Method"].eq("fixed_coordinate_bh")]
    assert not tested.empty
    assert tested["Sibling_Divergence_P_Value"].between(0.0, 1.0).all()
    assert np.allclose(
        tested["Sibling_Sparse_Evidence_P_Value"].astype(float),
        tested["Sibling_Divergence_P_Value"].astype(float),
    )
    assert tested["Sibling_Sparse_Evidence_Method"].eq("fixed_coordinate_bh").all()
    assert tested["Sibling_Gate_P_Value_Calibration"].eq("fixed_subspace_bh").all()
    assert tested["Sibling_Gate_P_Value_Role"].eq("active_traversal_sibling_gate").all()
    assert tested["Sibling_Dense_Evidence_Method"].eq("fixed_global_chi_square").all()
    assert tested["Sibling_Fixed_Block_BH_P_Value"].between(0.0, 1.0).all()
    assert tested["Sibling_Fixed_Global_P_Value"].between(0.0, 1.0).all()
    assert (tested["Sibling_Projection_Dimension"] == "").all() or tested[
        "Sibling_Projection_Dimension"
    ].isna().all()


def test_pipeline_supports_opt_in_fixed_global_sibling_gate() -> None:
    tree, annotations_df, leaf_data = _build_small_tree_with_leaf_data()

    bundle = run_gate_annotation_pipeline(
        tree,
        annotations_df.copy(),
        edge_alpha=0.01,
        sibling_alpha=0.01,
        leaf_data=leaf_data,
        sibling_gate_method="fixed_global_chi_square",
    )
    out = bundle.annotated_df

    assert bundle.metadata.config.sibling_gate_method == "fixed_global_chi_square"
    tested = out[out["Sibling_Test_Method"].eq("fixed_global_chi_square")]
    assert not tested.empty
    assert tested["Sibling_Divergence_P_Value"].between(0.0, 1.0).all()
    assert np.allclose(
        tested["Sibling_Dense_Evidence_P_Value"].astype(float),
        tested["Sibling_Divergence_P_Value"].astype(float),
    )
    assert tested["Sibling_Dense_Evidence_Method"].eq("fixed_global_chi_square").all()
    assert tested["Sibling_Gate_P_Value_Calibration"].eq("fixed_subspace_chi_square").all()
    assert tested["Sibling_Degrees_of_Freedom"].gt(0.0).all()


def test_pipeline_profile_avoids_adaptive_sibling_pca(monkeypatch) -> None:
    import tree_break_selection.hierarchy_analysis.decomposition.gates.orchestrator as orchestrator

    tree, annotations_df, leaf_data = _build_small_tree_with_leaf_data()

    def fail_if_called(*_args, **_kwargs):
        raise AssertionError("fixed profile must not resolve parent PCA inputs")

    stable_root = _stable_root

    monkeypatch.setattr(orchestrator, "_resolve_sibling_gate_inputs", fail_if_called)
    monkeypatch.setattr(
        orchestrator,
        "compute_root_feature_subsample_stability",
        stable_root,
    )

    bundle = orchestrator.run_gate_annotation_pipeline(
        tree,
        annotations_df.copy(),
        edge_alpha=0.01,
        sibling_alpha=0.01,
        leaf_data=leaf_data,
        sibling_gate_profile="fixed_coordinate_guarded_v1",
    )

    assert bundle.metadata.config.sibling_gate_profile_id == ("fixed_coordinate_guarded_v1")
    assert bundle.metadata.config.sibling_gate_method == "fixed_coordinate_bh"
    assert bundle.metadata.config.sibling_gate_alpha_penalty == 50.0
    assert bundle.metadata.config.root_stability_guard_threshold == 0.24
    assert bundle.metadata.config.root_stability_subsample_replicates == 12
    assert bundle.annotated_df["Sibling_Test_Method"].eq("fixed_coordinate_bh").any()


def test_fixed_sibling_gate_does_not_resolve_parent_pca_inputs(monkeypatch) -> None:
    import tree_break_selection.hierarchy_analysis.decomposition.gates.orchestrator as orchestrator

    tree, annotations_df, leaf_data = _build_small_tree_with_leaf_data()

    def fail_if_called(*_args, **_kwargs):
        raise AssertionError("fixed sibling gate must not resolve parent PCA inputs")

    monkeypatch.setattr(orchestrator, "_resolve_sibling_gate_inputs", fail_if_called)

    bundle = orchestrator.run_gate_annotation_pipeline(
        tree,
        annotations_df.copy(),
        edge_alpha=0.01,
        sibling_alpha=0.01,
        leaf_data=leaf_data,
        sibling_gate_method="fixed_coordinate_bh",
    )

    assert bundle.metadata.config.sibling_gate_method == "fixed_coordinate_bh"
    assert bundle.annotated_df["Sibling_Test_Method"].eq("fixed_coordinate_bh").any()


def test_sibling_gate_alpha_penalty_is_passed_as_effective_alpha(monkeypatch) -> None:
    import tree_break_selection.hierarchy_analysis.decomposition.gates.orchestrator as orchestrator

    tree, annotations_df, leaf_data = _build_small_tree_with_leaf_data()
    seen_alpha = None
    original = orchestrator.annotate_fixed_subspace_sibling_divergence

    def wrapped(*args, **kwargs):
        nonlocal seen_alpha
        seen_alpha = kwargs["significance_level_alpha"]
        return original(*args, **kwargs)

    monkeypatch.setattr(
        orchestrator,
        "annotate_fixed_subspace_sibling_divergence",
        wrapped,
    )

    bundle = orchestrator.run_gate_annotation_pipeline(
        tree,
        annotations_df.copy(),
        edge_alpha=0.01,
        sibling_alpha=0.01,
        sibling_gate_alpha_penalty=50.0,
        leaf_data=leaf_data,
        sibling_gate_method="fixed_coordinate_bh",
    )

    assert seen_alpha == 0.0002
    assert bundle.metadata.sibling.alpha == 0.01
    assert bundle.metadata.config.sibling_gate_alpha_penalty == 50.0


def test_pipeline_runs_opt_in_root_stability_guard(monkeypatch) -> None:
    import tree_break_selection.hierarchy_analysis.decomposition.gates.orchestrator as orchestrator

    tree, annotations_df, leaf_data = _build_small_tree_with_leaf_data()

    def unstable_root(*_args, **_kwargs):
        return {
            "root_stability_subsample_mean_ari": 0.0,
            "root_stability_subsample_median_ari": 0.0,
            "root_stability_subsample_q10_ari": 0.0,
        }

    monkeypatch.setattr(
        orchestrator,
        "compute_root_feature_subsample_stability",
        unstable_root,
    )

    bundle = orchestrator.run_gate_annotation_pipeline(
        tree,
        annotations_df.copy(),
        edge_alpha=0.01,
        sibling_alpha=0.01,
        leaf_data=leaf_data,
        sibling_gate_method="fixed_coordinate_bh",
        root_stability_guard_threshold=0.50,
        root_stability_subsample_replicates=4,
        root_stability_feature_fraction=0.8,
        root_stability_seed=7,
    )

    assert bundle.metadata.config.root_stability_guard_threshold == 0.50
    assert bundle.metadata.config.root_stability_subsample_replicates == 4
    assert bundle.metadata.config.root_stability_feature_fraction == 0.8
    assert bundle.metadata.config.root_stability_seed == 7
    assert bundle.metadata.config.root_stability_tree_distance_metric == "hamming"
    assert bundle.metadata.config.root_stability_tree_linkage_method == "average"
    assert "Root_Stability_Guard_Blocked" in bundle.annotated_df.columns


def test_pipeline_runs_opt_in_root_selective_permutation_guard(monkeypatch) -> None:
    import tree_break_selection.hierarchy_analysis.decomposition.gates.orchestrator as orchestrator

    tree, annotations_df, leaf_data = _build_small_tree_with_leaf_data()
    feature_space = infer_feature_space_from_columns(tuple(leaf_data.columns))

    fake_selected_root = _selected_root(0.50)

    monkeypatch.setattr(
        guard_module,
        "selected_root_permutation_p_value",
        fake_selected_root,
    )

    bundle = orchestrator.run_gate_annotation_pipeline(
        tree,
        annotations_df.copy(),
        edge_alpha=0.01,
        sibling_alpha=0.01,
        leaf_data=leaf_data,
        feature_space=feature_space,
        sibling_gate_method="fixed_coordinate_bh",
        root_selective_permutation_guard_replicates=5,
        root_selective_permutation_guard_seed=123,
        root_selective_permutation_guard_alpha=0.01,
    )

    assert bundle.metadata.config.root_selective_permutation_guard_replicates == 5
    assert bundle.metadata.config.root_selective_permutation_guard_seed == 123
    assert bundle.metadata.config.root_selective_permutation_guard_alpha == 0.01
    assert "Root_Selective_Permutation_Guard_Blocked" in bundle.annotated_df.columns


def test_pipeline_rejects_root_selective_guard_with_adaptive_sibling_gate() -> None:
    tree, annotations_df, leaf_data = _build_small_tree_with_leaf_data()

    with pytest.raises(ValueError, match="fixed-subspace sibling gate"):
        run_gate_annotation_pipeline(
            tree,
            annotations_df.copy(),
            edge_alpha=0.01,
            sibling_alpha=0.01,
            leaf_data=leaf_data,
            root_selective_permutation_guard_replicates=5,
        )


def test_pipeline_rejects_inert_root_stability_guard_config() -> None:
    tree, annotations_df, leaf_data = _build_small_tree_with_leaf_data()

    with pytest.raises(ValueError, match="root_stability_subsample_replicates"):
        run_gate_annotation_pipeline(
            tree,
            annotations_df.copy(),
            edge_alpha=0.01,
            sibling_alpha=0.01,
            leaf_data=leaf_data,
            sibling_gate_method="fixed_coordinate_bh",
            root_stability_guard_threshold=0.50,
            root_stability_subsample_replicates=0,
        )


def test_resolve_effective_sibling_alpha_rejects_invalid_penalty() -> None:
    assert resolve_effective_sibling_alpha(0.01, 50.0) == 0.0002
    with pytest.raises(ValueError, match="sibling_gate_alpha_penalty"):
        resolve_effective_sibling_alpha(0.01, 0.0)


def test_pipeline_rejects_unknown_sibling_gate_method() -> None:
    tree, annotations_df, leaf_data = _build_small_tree_with_leaf_data()

    with pytest.raises(ValueError, match="Unknown sibling_gate_method"):
        run_gate_annotation_pipeline(
            tree,
            annotations_df.copy(),
            edge_alpha=0.01,
            sibling_alpha=0.01,
            leaf_data=leaf_data,
            sibling_gate_method="adaptive_magic",
        )


def test_pipeline_requires_leaf_data_for_spectral_gate_context() -> None:
    tree, annotations_df, _leaf_data = _build_small_tree_with_leaf_data()

    with pytest.raises(ValueError, match="require leaf_data"):
        run_gate_annotation_pipeline(
            tree,
            annotations_df.copy(),
            edge_alpha=0.01,
            sibling_alpha=0.01,
            leaf_data=None,
        )
