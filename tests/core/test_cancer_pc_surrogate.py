"""Cancer-inspired singleton root and matched controls through native dispatch."""

import numpy as np
import pytest
from scipy.sparse.csgraph import connected_components

from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_replay import _geometry
from benchmarks.shared.cases import get_test_cases_by_suite
from benchmarks.shared.cases.geometry import case_recipe_geometry
from benchmarks.shared.generators.generate_case_data import generate_case_data
from benchmarks.shared.runners.dispatch import run_clustering_result
from benchmarks.shared.runners.method_registry import METHOD_SPECS
from benchmarks.shared.util.case_inputs import prepare_case_inputs


def cases():
    return {
        c["name"]: c
        for c in get_test_cases_by_suite("method_proof")
        if c["generator"] == "cancer_pc_surrogate"
    }


def test_cancer_pc_suite_selects_the_four_registered_cases():
    selected = get_test_cases_by_suite("cancer_pc")
    assert len(selected) == 4
    assert {case["name"]: case for case in selected} == cases()
    full = {case["name"]: case for case in get_test_cases_by_suite("full")}
    assert all(full[case["name"]] == case for case in selected)


def test_registered_surrogate_controls_preserve_random_draws_and_truth():
    recipes = cases()
    assert set(recipes) == {
        f"cancer_pc_{hub}_{graph}"
        for hub in ("hub", "nohub")
        for graph in ("disconnected", "connected")
    }
    generated = {name: generate_case_data(case) for name, case in recipes.items()}
    x, y, original, metadata = generated["cancer_pc_hub_disconnected"]
    control = generated["cancer_pc_nohub_disconnected"][0]
    assert case_recipe_geometry(recipes["cancer_pc_hub_disconnected"]) == (496, 30)
    np.testing.assert_array_equal(np.bincount(y), [256, 64, 32, 16, 64, 64])
    np.testing.assert_array_equal(x, original)
    np.testing.assert_array_equal(x.iloc[1:], control.iloc[1:])
    np.testing.assert_allclose(x.iloc[0], 0.03 * control.iloc[0])
    for name, (frame, truth, _, _) in generated.items():
        np.testing.assert_array_equal(truth, y)
        np.testing.assert_array_equal(frame, generate_case_data(recipes[name])[0])
    connected = generated["cancer_pc_hub_connected"][0]
    np.testing.assert_array_equal(x.iloc[:368], connected.iloc[:368])
    np.testing.assert_allclose(x.iloc[368:432, 1] - connected.iloc[368:432, 1], -17)
    np.testing.assert_allclose(x.iloc[432:, 1] - connected.iloc[432:, 1], 17)
    assert metadata["simulation_model"] == "gaussian_mixture_with_one_low_variance_row"
    assert metadata["hub_sample"] == "S0000"
    assert "counts" in metadata["scientific_caution"]
    prepared = prepare_case_inputs(recipes["cancer_pc_hub_disconnected"], ["tbs_family_diffusion"])
    assert prepared.data.equals(x)
    assert prepared.distance_condensed.shape == (496 * 495 // 2,)


@pytest.mark.parametrize("hub,root_sizes", [("hub", [1, 495]), ("nohub", [47, 449])])
def test_disconnected_surrogate_passes_through_the_singleton_root_like_its_control(hub, root_sizes):
    recipe = cases()[f"cancer_pc_{hub}_disconnected"]
    data, _, _, metadata = generate_case_data(recipe)
    signature, state = _geometry(data.to_numpy(), data.index.tolist())
    graph = state["scaled_distances"].maximum(state["scaled_distances"].T)
    graph.setdiag(0)
    graph.eliminate_zeros()
    _, components = connected_components(graph, directed=False)
    assert sorted(np.bincount(components).tolist()) == [64, 64, 368]
    assert np.count_nonzero(np.isclose(state["eigenvalues"], 1, atol=1e-10, rtol=0)) == 3
    result = run_clustering_result(
        data,
        "tbs_family_diffusion",
        dict(METHOD_SPECS["tbs_family_diffusion"].param_grid[0]),
        feature_space=metadata["feature_space"],
    )
    extra = result.extra
    z = extra["linkage_matrix"]
    sizes = np.r_[np.ones(len(data)), z[:, 3]]
    assert sorted(sizes[z[-1, :2].astype(int)].astype(int).tolist()) == root_sizes
    assert tuple(sorted(z[-1, :2].astype(int))) == signature[-1]
    assert extra["branch_length_optimization_applied_to_tree"]
    assert result.found_clusters == 6
    # Only the hub's 1/495 root is a singleton peel that traversal must pass through.
    assert extra["traversal_counters"]["live_pass_through_count"] == int(hub == "hub")
    if hub == "hub":
        assert 0 in signature[-1]
        assert graph[0].count_nonzero() > 200


def test_surrogate_rejects_nonfinite_noise_scale():
    recipe = cases()["cancer_pc_hub_disconnected"]
    with pytest.raises(ValueError, match="hub_scale"):
        generate_case_data(dict(recipe, hub_scale=float("nan")))
