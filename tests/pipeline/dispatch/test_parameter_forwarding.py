from __future__ import annotations

from benchmarks.shared.runners.dispatch import run_clustering_result
from benchmarks.shared.runners.method_registry import METHOD_SPECS
from benchmarks.shared.types.method_spec import MethodSpec
from scipy.spatial.distance import pdist
from tree_break_selection.tree.feature_space import continuous_feature_space_from_columns

from .helpers import _appending_runner, _capturing_runner, _toy_dataframe


def test_run_clustering_result_forwards_tbs_gate_profile_params(monkeypatch):
    captured = {}
    capture_runner = _capturing_runner(captured)

    monkeypatch.setitem(
        METHOD_SPECS,
        "tbs",
        MethodSpec(
            name="TBS Divergence",
            runner=capture_runner,
            param_grid=[
                {
                    "tree_distance_metric": "hamming",
                    "tree_linkage_method": "average",
                }
            ],
        ),
    )

    df = _toy_dataframe()
    run_clustering_result(
        data_df=df,
        method_id="tbs",
        params={
            "tree_distance_metric": "euclidean",
            "tree_linkage_method": "average",
            "sibling_gate_profile": "fixed_coordinate_global_passthrough_refined_v1",
            "sibling_gate_method": "fixed_coordinate_bh",
            "sibling_gate_alpha_penalty": 50.0,
            "root_stability_guard_threshold": 0.24,
            "root_stability_subsample_replicates": 12,
            "root_stability_feature_fraction": 0.8,
            "root_stability_seed": 7,
            "root_stability_tree_distance_metric": "jaccard",
            "root_stability_tree_linkage_method": "complete",
            "root_selective_permutation_guard_replicates": 99,
            "root_selective_permutation_guard_seed": 11,
            "root_selective_permutation_guard_alpha": 0.01,
            "root_selective_permutation_guard_scope": (
                "global_sibling_min_passthrough_descendant_refined"
            ),
            "root_selective_permutation_guard_tree_distance_metric": ("rogerstanimoto"),
            "root_selective_permutation_guard_tree_linkage_method": "weighted",
            "spectral_transport_passthrough_guard": True,
            "spectral_transport_max_cost": 0.75,
            "spectral_transport_require_mp_blocks": False,
            "spectral_transport_block_log_tolerance": 0.02,
            "spectral_transport_unmatched_mode_penalty": 1.5,
            "spectral_include_internal_barycenters": True,
            "spectral_internal_distribution_mode": "branch_length_state",
            "continuous_covariance_policy": "guarded_within_child",
            "continuous_covariance_min_child_leaf_count": 8,
            "neighborhood_bandwidth_profile": ("regional_tau_branch_length_support_only_v1"),
            "enforce_internal_support_thresholds": True,
            "passthrough": True,
        },
        seed=42,
        distance_condensed=pdist(df.values, metric="euclidean"),
    )

    assert captured["kwargs"]["sibling_gate_profile"] == (
        "fixed_coordinate_global_passthrough_refined_v1"
    )
    assert captured["kwargs"]["sibling_gate_method"] == "fixed_coordinate_bh"
    assert captured["kwargs"]["sibling_gate_alpha_penalty"] == 50.0
    assert captured["kwargs"]["root_stability_guard_threshold"] == 0.24
    assert captured["kwargs"]["root_stability_subsample_replicates"] == 12
    assert captured["kwargs"]["root_stability_feature_fraction"] == 0.8
    assert captured["kwargs"]["root_stability_seed"] == 7
    assert captured["kwargs"]["root_stability_tree_distance_metric"] == "jaccard"
    assert captured["kwargs"]["root_stability_tree_linkage_method"] == "complete"
    assert captured["kwargs"]["root_selective_permutation_guard_replicates"] == 99
    assert captured["kwargs"]["root_selective_permutation_guard_seed"] == 11
    assert captured["kwargs"]["root_selective_permutation_guard_alpha"] == 0.01
    assert captured["kwargs"]["root_selective_permutation_guard_scope"] == (
        "global_sibling_min_passthrough_descendant_refined"
    )
    assert (
        captured["kwargs"]["root_selective_permutation_guard_tree_distance_metric"]
        == "rogerstanimoto"
    )
    assert captured["kwargs"]["root_selective_permutation_guard_tree_linkage_method"] == "weighted"
    assert captured["kwargs"]["spectral_transport_passthrough_guard"] is True
    assert captured["kwargs"]["spectral_transport_max_cost"] == 0.75
    assert captured["kwargs"]["spectral_transport_require_mp_blocks"] is False
    assert captured["kwargs"]["spectral_transport_block_log_tolerance"] == 0.02
    assert captured["kwargs"]["spectral_transport_unmatched_mode_penalty"] == 1.5
    assert captured["kwargs"]["spectral_include_internal_barycenters"] is True
    assert captured["kwargs"]["spectral_internal_distribution_mode"] == ("branch_length_state")
    assert captured["kwargs"]["continuous_covariance_policy"] == "guarded_within_child"
    assert captured["kwargs"]["continuous_covariance_min_child_leaf_count"] == 8
    assert captured["kwargs"]["neighborhood_bandwidth_profile"] == (
        "regional_tau_branch_length_support_only_v1"
    )
    assert captured["kwargs"]["enforce_internal_support_thresholds"] is True
    assert captured["kwargs"]["passthrough"] is True


def test_run_clustering_result_uses_continuous_sibling_gate_only_for_continuous_blocks(
    monkeypatch,
):
    captured_calls = []
    capture_runner = _appending_runner(captured_calls)

    monkeypatch.setitem(
        METHOD_SPECS,
        "tbs",
        MethodSpec(
            name="TBS Adaptive Gate Test",
            runner=capture_runner,
            param_grid=[
                {
                    "tree_distance_metric": "euclidean",
                    "tree_linkage_method": "average",
                }
            ],
        ),
    )

    df = _toy_dataframe()
    params = {
        "tree_distance_metric": "euclidean",
        "tree_linkage_method": "average",
        "continuous_sibling_gate_method": "fixed_coordinate_bh",
    }
    distance_condensed = pdist(df.values, metric="euclidean")
    run_clustering_result(
        data_df=df,
        method_id="tbs",
        params=params,
        seed=42,
        distance_condensed=distance_condensed,
        feature_space=None,
    )
    run_clustering_result(
        data_df=df,
        method_id="tbs",
        params=params,
        seed=42,
        distance_condensed=distance_condensed,
        feature_space=continuous_feature_space_from_columns(tuple(df.columns)),
    )

    assert captured_calls[0]["kwargs"]["sibling_gate_method"] == ("projected_wald_inflation")
    assert captured_calls[1]["kwargs"]["sibling_gate_method"] == "fixed_coordinate_bh"


def test_run_clustering_result_forwards_adaptive_nnls_branch_time_params(monkeypatch):
    captured = {}
    original_spec = METHOD_SPECS["tbs_diffusion_adaptive_nnls"]
    capture_runner = _capturing_runner(captured)

    monkeypatch.setitem(
        METHOD_SPECS,
        "tbs_diffusion_adaptive_nnls",
        MethodSpec(
            name=original_spec.name,
            runner=capture_runner,
            param_grid=original_spec.param_grid,
        ),
    )

    df = _toy_dataframe()
    result = run_clustering_result(
        data_df=df,
        method_id="tbs_diffusion_adaptive_nnls",
        params=original_spec.param_grid[0],
        seed=42,
    )

    assert result.status == "ok"
    assert captured["kwargs"]["edge_branch_length_variance_policy"] == ("normalized_branch_length")
    assert captured["kwargs"]["branch_length_optimization_method"] == "fixed_topology_nnls"
    assert captured["kwargs"]["branch_length_optimization_target_metric"] == (
        "squared_standardized_euclidean"
    )
    assert captured["kwargs"]["branch_length_optimization_pair_sample_size"] == 50_000
    assert captured["kwargs"]["branch_length_optimization_random_state"] == 0
    assert captured["kwargs"]["branch_length_optimization_solver_tolerance"] == 1e-5
    assert captured["kwargs"]["branch_length_optimization_max_iterations"] == 1000


def test_run_clustering_result_forwards_graphtools_adaptive_k_params(monkeypatch):
    captured = {}
    original_spec = METHOD_SPECS["tbs_diffusion_graphtools_adaptive_nnls"]
    capture_runner = _capturing_runner(captured)

    monkeypatch.setitem(
        METHOD_SPECS,
        "tbs_diffusion_graphtools_adaptive_nnls",
        MethodSpec(
            name=original_spec.name,
            runner=capture_runner,
            param_grid=original_spec.param_grid,
        ),
    )

    result = run_clustering_result(
        data_df=_toy_dataframe(),
        method_id="tbs_diffusion_graphtools_adaptive_nnls",
        params=original_spec.param_grid[0],
        seed=42,
    )

    assert result.status == "ok"
    assert captured["kwargs"]["adaptive_neighbor_profile"] == "fragmentation_guard"
    assert captured["kwargs"]["adaptive_neighbor_grid"] == (5, 10, 15, 25, 40, 80, 160)
    assert captured["kwargs"]["tree_builder"] == "linkage"
    assert captured["kwargs"]["tree_rooting"] == "linkage_root"
    assert captured["kwargs"]["tree_linkage_method"] == "average"
    assert captured["kwargs"]["edge_branch_length_variance_policy"] == ("normalized_branch_length")
    assert captured["kwargs"]["branch_length_optimization_method"] == ("fixed_topology_nnls")


def test_run_clustering_result_forwards_graphtools_neighbor_joining_tree_params(
    monkeypatch,
):
    captured = {}
    original_spec = METHOD_SPECS["tbs_diffusion_graphtools_adaptive_nnls"]
    neighbor_joining_params = next(
        params
        for params in original_spec.param_grid
        if params["tree_builder"] == "neighbor_joining"
    )
    capture_runner = _capturing_runner(captured)

    monkeypatch.setitem(
        METHOD_SPECS,
        "tbs_diffusion_graphtools_adaptive_nnls",
        MethodSpec(
            name=original_spec.name,
            runner=capture_runner,
            param_grid=original_spec.param_grid,
        ),
    )

    result = run_clustering_result(
        data_df=_toy_dataframe(),
        method_id="tbs_diffusion_graphtools_adaptive_nnls",
        params=neighbor_joining_params,
        seed=42,
    )

    assert result.status == "ok"
    assert captured["kwargs"]["tree_builder"] == "neighbor_joining"
    assert captured["kwargs"]["tree_rooting"] == "mad"
    assert captured["kwargs"]["tree_linkage_method"] == "average"
    assert captured["kwargs"]["adaptive_neighbor_profile"] == "fragmentation_guard"
