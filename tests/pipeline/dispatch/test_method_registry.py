from __future__ import annotations

import pytest
from benchmarks.shared.runners.dispatch import run_clustering_result
from benchmarks.shared.runners.method_registry import METHOD_SPECS
from benchmarks.shared.util.method_sets import TBS_RUNNER_METHODS
from scipy.spatial.distance import pdist

from .helpers import _toy_dataframe


@pytest.mark.parametrize(
    "expected_profiles",
    [
        pytest.param(
            {
                "tbs_global_passthrough_refined_diagnostic": (
                    "fixed_coordinate_global_passthrough_refined_v1"
                )
            },
            id="global-passthrough-refined",
        ),
    ],
)
def test_method_registry_exposes_sibling_gate_profiles(
    expected_profiles: dict[str, str],
) -> None:
    for method_id, profile_id in expected_profiles.items():
        params = METHOD_SPECS[method_id].param_grid[0]
        assert params["sibling_gate_profile"] == profile_id
        assert params["tree_distance_metric"] == "hamming"
        assert params["tree_linkage_method"] == "average"


def test_method_registry_exposes_fixed_fdr_benchmark_variants():
    expected_methods = {
        "tbs_fixed_coordinate_bh": "fixed_coordinate_bh",
        "tbs_fixed_coordinate_by": "fixed_coordinate_by",
        "tbs_fixed_coordinate_holm": "fixed_coordinate_holm",
        "tbs_fixed_coordinate_bonferroni": "fixed_coordinate_bonferroni",
        "tbs_fixed_block_bh": "fixed_block_bh",
        "tbs_fixed_block_simes_bh": "fixed_block_simes_bh",
    }

    for method_id, sibling_gate_method in expected_methods.items():
        spec = METHOD_SPECS[method_id]
        params = spec.param_grid[0]
        assert method_id in TBS_RUNNER_METHODS
        assert params["sibling_gate_method"] == sibling_gate_method
        assert params["tree_distance_metric"] == "hamming"
        assert params["tree_linkage_method"] == "average"


def test_method_registry_exposes_continuous_guarded_covariance_candidate():
    spec = METHOD_SPECS["tbs_continuous_guarded_within_covariance"]
    params = spec.param_grid[0]

    assert "tbs_continuous_guarded_within_covariance" in TBS_RUNNER_METHODS
    assert params["tree_distance_metric"] == "hamming"
    assert params["tree_linkage_method"] == "average"
    assert params["continuous_covariance_policy"] == "guarded_within_child"
    assert params["continuous_covariance_min_child_leaf_count"] == 8
    assert "sibling_gate_method" not in params
    assert params["continuous_sibling_gate_method"] == "fixed_coordinate_bh"


def test_continuous_guarded_covariance_candidate_runs_non_continuous_input():
    df = _toy_dataframe()
    params = METHOD_SPECS["tbs_continuous_guarded_within_covariance"].param_grid[0]

    result = run_clustering_result(
        data_df=df,
        method_id="tbs_continuous_guarded_within_covariance",
        params=params,
        seed=42,
        distance_condensed=pdist(df.values, metric="euclidean"),
        feature_space=None,
    )

    assert result.status == "ok"
    assert result.labels is not None


def test_method_registry_names_diffusion_methods_and_branch_lengths_explicitly():
    nn_diffusion = METHOD_SPECS["tbs_diffusion"]
    nn_params = nn_diffusion.param_grid[0]
    adaptive_diffusion = METHOD_SPECS["tbs_diffusion_adaptive"]
    adaptive_params = adaptive_diffusion.param_grid[0]
    adaptive_nnls_diffusion = METHOD_SPECS["tbs_diffusion_adaptive_nnls"]
    adaptive_nnls_params = adaptive_nnls_diffusion.param_grid[0]
    graphtools_diffusion = METHOD_SPECS["tbs_diffusion_graphtools"]
    graphtools_params = graphtools_diffusion.param_grid[0]
    graphtools_nnls_diffusion = METHOD_SPECS["tbs_diffusion_graphtools_nnls"]
    graphtools_nnls_params = graphtools_nnls_diffusion.param_grid[0]
    graphtools_adaptive_nnls_diffusion = METHOD_SPECS["tbs_diffusion_graphtools_adaptive_nnls"]
    graphtools_adaptive_nnls_params = graphtools_adaptive_nnls_diffusion.param_grid[0]

    assert nn_diffusion.name == "TBS (Hamming NN Diffusion)"
    assert nn_params["diffusion_method"] == "hamming_nn_diffusion"
    assert nn_params["k_neighbors"] == 15
    assert nn_params["tree_linkage_method"] == "average"
    assert nn_params["branch_length_optimization_method"] == "linkage_ultrametric"

    assert adaptive_diffusion.name == "TBS (Adaptive pydiffmap Diffusion)"
    assert adaptive_params["diffusion_method"] == "adaptive_pydiffmap_diffusion"
    assert adaptive_params["bandwidth_type"] == "-1/(d+2)"
    assert adaptive_params["epsilon"] == "median"
    assert adaptive_params["tree_linkage_method"] == "average"
    assert adaptive_params["branch_length_optimization_method"] == "linkage_ultrametric"

    assert adaptive_nnls_diffusion.name == ("TBS (Adaptive pydiffmap Diffusion, NNLS Branch-Time)")
    assert adaptive_nnls_params["diffusion_method"] == "adaptive_pydiffmap_diffusion"
    assert adaptive_nnls_params["bandwidth_type"] == "-1/(d+2)"
    assert adaptive_nnls_params["epsilon"] == "median"
    assert adaptive_nnls_params["edge_branch_length_variance_policy"] == (
        "normalized_branch_length"
    )
    assert adaptive_nnls_params["branch_length_optimization_method"] == "fixed_topology_nnls"
    assert adaptive_nnls_params["branch_length_optimization_target_metric"] == (
        "squared_standardized_euclidean"
    )
    assert adaptive_nnls_params["branch_length_optimization_pair_sample_size"] == 50_000
    assert adaptive_nnls_params["branch_length_optimization_random_state"] == 0
    assert adaptive_nnls_params["branch_length_optimization_solver_tolerance"] == 1e-5
    assert adaptive_nnls_params["branch_length_optimization_max_iterations"] == 1000

    assert graphtools_diffusion.name == "TBS (graphtools Kernel Diffusion)"
    assert graphtools_params["diffusion_method"] == "graphtools_kernel_diffusion"
    assert graphtools_params["k_neighbors"] == 10
    assert graphtools_params["decay"] == 40
    assert graphtools_params["kernel_symm"] == "+"
    assert graphtools_params["branch_length_optimization_method"] == "linkage_ultrametric"

    assert graphtools_nnls_diffusion.name == ("TBS (graphtools Kernel Diffusion, NNLS Branch-Time)")
    assert graphtools_nnls_params["diffusion_method"] == "graphtools_kernel_diffusion"
    assert graphtools_nnls_params["k_neighbors"] == 10
    assert graphtools_nnls_params["decay"] == 40
    assert graphtools_nnls_params["kernel_symm"] == "+"
    assert graphtools_nnls_params["edge_branch_length_variance_policy"] == (
        "normalized_branch_length"
    )
    assert graphtools_nnls_params["branch_length_optimization_method"] == ("fixed_topology_nnls")
    assert graphtools_nnls_params["branch_length_optimization_pair_sample_size"] == 50_000

    assert graphtools_adaptive_nnls_diffusion.name == (
        "TBS (graphtools Kernel Diffusion, Adaptive-K NNLS Branch-Time)"
    )
    assert graphtools_adaptive_nnls_params["diffusion_method"] == ("graphtools_kernel_diffusion")
    assert graphtools_adaptive_nnls_params["k_neighbors"] == 10
    assert graphtools_adaptive_nnls_params["adaptive_neighbor_profile"] == ("fragmentation_guard")
    assert graphtools_adaptive_nnls_params["adaptive_neighbor_grid"] == (
        5,
        10,
        15,
        25,
        40,
        80,
        160,
    )
    assert graphtools_adaptive_nnls_params["branch_length_optimization_method"] == (
        "fixed_topology_nnls"
    )
    assert graphtools_adaptive_nnls_params["tree_builder"] == "linkage"
    assert graphtools_adaptive_nnls_params["tree_linkage_method"] == "average"

    tree_strategy_grid = graphtools_adaptive_nnls_diffusion.param_grid
    assert len(tree_strategy_grid) == 8
    assert len({params["benchmark_run_id"] for params in tree_strategy_grid}) == 8
    assert {params["benchmark_class"] for params in tree_strategy_grid} == {"optional_gpl"}
    assert {params["benchmark_grid"] for params in tree_strategy_grid} == {
        "graphtools_adaptive_k_tree_strategy"
    }
    assert {params["benchmark_repeat"] for params in tree_strategy_grid} == {0}

    linkage_methods = {
        params["tree_linkage_method"]
        for params in tree_strategy_grid
        if params["tree_builder"] == "linkage"
    }
    assert linkage_methods == {
        "average",
        "complete",
        "weighted",
        "single",
        "centroid",
        "median",
        "ward",
    }
    for params in tree_strategy_grid:
        assert params["branch_length_optimization_method"] == "fixed_topology_nnls"

    graphtools_adaptive_nnls_neighbor_joining_params = next(
        params for params in tree_strategy_grid if params["tree_builder"] == "neighbor_joining"
    )
    assert graphtools_adaptive_nnls_neighbor_joining_params["tree_builder"] == ("neighbor_joining")
    assert graphtools_adaptive_nnls_neighbor_joining_params["tree_rooting"] == "mad"
    assert (
        graphtools_adaptive_nnls_neighbor_joining_params["branch_length_optimization_method"]
        == "fixed_topology_nnls"
    )


def test_method_registry_exposes_current_support_profiles():
    current_profile_ids = {
        "tbs_internal_filter_v1",
        "tbs_internal_filter_branch_length_v1",
        "tbs_bandwidth_context_v1",
    }
    internal = METHOD_SPECS["tbs_internal_filter_v1"].param_grid[0]
    branch_length = METHOD_SPECS["tbs_internal_filter_branch_length_v1"].param_grid[0]
    bandwidth = METHOD_SPECS["tbs_bandwidth_context_v1"].param_grid[0]

    assert current_profile_ids.issubset(METHOD_SPECS)
    assert current_profile_ids.issubset(TBS_RUNNER_METHODS)
    assert internal["spectral_include_internal_barycenters"] is True
    assert internal["spectral_internal_distribution_mode"] == "empirical_barycenter"
    assert internal["enforce_internal_support_thresholds"] is True
    assert branch_length["spectral_internal_distribution_mode"] == ("branch_length_state")
    assert branch_length["enforce_internal_support_thresholds"] is True
    assert bandwidth["neighborhood_bandwidth_profile"] == (
        "regional_tau_branch_length_support_only_v1"
    )
