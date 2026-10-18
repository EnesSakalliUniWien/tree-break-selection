from __future__ import annotations

import inspect
import re

import pytest
from scipy.spatial.distance import pdist

from benchmarks.shared.runners.dispatch import run_clustering_result
from benchmarks.shared.runners.method_registry import METHOD_SPECS
from benchmarks.shared.types.method_spec import MethodSpec
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
            "spectral_include_internal_barycenters": True,
            "spectral_internal_distribution_mode": "branch_length_state",
            "continuous_covariance_policy": "guarded_within_child",
            "continuous_covariance_min_child_leaf_count": 8,
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
    assert captured["kwargs"]["spectral_include_internal_barycenters"] is True
    assert captured["kwargs"]["spectral_internal_distribution_mode"] == ("branch_length_state")
    assert captured["kwargs"]["continuous_covariance_policy"] == "guarded_within_child"
    assert captured["kwargs"]["continuous_covariance_min_child_leaf_count"] == 8
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

    # Unset: the gate pipeline resolves it to the profile's method or projected_wald_inflation.
    assert captured_calls[0]["kwargs"]["sibling_gate_method"] is None
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
        "mean_squared_standardized_euclidean"
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
    assert captured["kwargs"]["tree_linkage_method"] == "average"
    assert captured["kwargs"]["edge_branch_length_variance_policy"] == ("normalized_branch_length")
    assert captured["kwargs"]["branch_length_optimization_method"] == ("fixed_topology_nnls")


def test_run_clustering_result_forwards_graphtools_tree_strategy_params(
    monkeypatch,
):
    captured = {}
    original_spec = METHOD_SPECS["tbs_diffusion_graphtools_adaptive_nnls"]
    ward_params = next(
        params for params in original_spec.param_grid if params.get("tree_linkage_method") == "ward"
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
        params=ward_params,
        seed=42,
    )

    assert result.status == "ok"
    assert captured["kwargs"]["tree_linkage_method"] == "ward"
    assert captured["kwargs"]["adaptive_neighbor_profile"] == "fragmentation_guard"


def test_run_clustering_result_rejects_unknown_tbs_params():
    with pytest.raises(ValueError, match="sibling_gate_alpha_penaly"):
        run_clustering_result(
            data_df=_toy_dataframe(),
            method_id="tbs",
            params={
                "tree_distance_metric": "hamming",
                "tree_linkage_method": "average",
                "sibling_gate_alpha_penaly": 50.0,
            },
        )


def test_tbs_param_allowlist_matches_keys_read_by_dispatch():
    from benchmarks.shared.runners import dispatch

    source = "".join(
        inspect.getsource(function)
        for function in (
            dispatch._run_tbs_distance_tree_method,
            dispatch._tbs_branch_length_optimization_kwargs,
            dispatch._resolve_tbs_sibling_gate_method,
        )
    )
    read_keys = set(re.findall(r'params(?:\.get\(|\[)\s*"([a-z0-9_]+)"', source))
    assert read_keys == dispatch._TBS_DISTANCE_TREE_PARAM_KEYS


def test_every_registered_method_allows_its_registry_params():
    from benchmarks.shared.benchmark_grid import BENCHMARK_METADATA_KEYS
    from benchmarks.shared.runners import dispatch

    for method_id, spec in METHOD_SPECS.items():
        allowed = dispatch._PARAM_KEYS_BY_METHOD[method_id] | BENCHMARK_METADATA_KEYS
        for params in spec.param_grid:
            assert set(params) <= allowed, (method_id, sorted(set(params) - allowed))


def test_run_clustering_result_rejects_unknown_diffusion_params():
    with pytest.raises(ValueError, match="edge_branch_length_variance_polcy"):
        run_clustering_result(
            data_df=_toy_dataframe(),
            method_id="tbs_diffusion_adaptive",
            params={
                **METHOD_SPECS["tbs_diffusion_adaptive"].param_grid[0],
                "edge_branch_length_variance_polcy": "normalized",
            },
        )


def _row_capturing_runner(captured: dict[str, object]):
    """Capture runner arguments and return one label per input row."""
    import numpy as np

    from benchmarks.shared.types import MethodRunResult

    def runner(*args, **kwargs):
        captured["args"], captured["kwargs"] = args, kwargs
        return MethodRunResult(
            labels=np.zeros(len(args[0]), dtype=int),
            found_clusters=1,
            report_df=None,
            status="ok",
            skip_reason=None,
            extra={},
        )

    return runner


def _two_cluster_continuous_frame(n_samples: int = 40, n_features: int = 6):
    import numpy as np
    import pandas as pd

    values = np.random.default_rng(0).normal(size=(n_samples, n_features))
    values[n_samples // 2 :] += 3.0
    return pd.DataFrame(
        values,
        index=[f"S{i}" for i in range(n_samples)],
        columns=[f"F{j}" for j in range(n_features)],
    )


def test_family_diffusion_tbs_uses_diffusion_over_the_case_distance(monkeypatch):
    import numpy as np

    captured = {}
    original_spec = METHOD_SPECS["tbs_family_diffusion"]
    monkeypatch.setitem(
        METHOD_SPECS,
        "tbs_family_diffusion",
        MethodSpec(
            name=original_spec.name,
            runner=_row_capturing_runner(captured),
            param_grid=original_spec.param_grid,
        ),
    )
    data = _two_cluster_continuous_frame()
    case_distance = np.ones(len(data) * (len(data) - 1) // 2)

    result = run_clustering_result(
        data_df=data,
        method_id="tbs_family_diffusion",
        params=original_spec.param_grid[0],
        seed=0,
        distance_condensed=case_distance,
        feature_space=continuous_feature_space_from_columns(tuple(data.columns)),
    )

    assert result.status == "ok"
    assert not np.allclose(captured["args"][1], case_distance)
    assert captured["kwargs"]["adaptive_projection_dimension_energy_fraction"] == 0.9
    assert captured["kwargs"]["branch_length_data_df"] is data
    assert captured["kwargs"]["extra"]["feature_family_diffusion_fallback_reason"] is None


def test_family_diffusion_tbs_falls_back_to_the_case_distance_on_duplicate_rows(monkeypatch):
    import numpy as np
    import pandas as pd

    captured = {}
    original_spec = METHOD_SPECS["tbs_family_diffusion"]
    monkeypatch.setitem(
        METHOD_SPECS,
        "tbs_family_diffusion",
        MethodSpec(
            name=original_spec.name,
            runner=_row_capturing_runner(captured),
            param_grid=original_spec.param_grid,
        ),
    )
    values = (np.random.default_rng(1).random((40, 12)) < 0.5).astype(float)
    values[:8] = values[0]
    data = pd.DataFrame(values, index=[f"S{i}" for i in range(40)])
    case_distance = pdist(values, metric="hamming")

    run_clustering_result(
        data_df=data,
        method_id="tbs_family_diffusion",
        params=original_spec.param_grid[0],
        seed=0,
        distance_condensed=case_distance,
    )

    np.testing.assert_allclose(captured["args"][1], case_distance)
    assert captured["kwargs"]["adaptive_projection_dimension_energy_fraction"] is None
    assert (
        "identical rows" in captured["kwargs"]["extra"]["feature_family_diffusion_fallback_reason"]
    )
