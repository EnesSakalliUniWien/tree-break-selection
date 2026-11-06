"""Ensure the diagnostic retains the production runner's NNLS semantics."""

import numpy as np
import pytest
from sklearn.metrics import adjusted_rand_score

from benchmarks.diagnostics.calibration.diffusion_factorial import datasets, evaluate
from benchmarks.shared.runners.tbs_runner import run_tbs_on_distance
from tree_break_selection.space_separation.diffusion import feature_family_diffusion_geometry
from tree_break_selection.tree.construction import build_tree
from tree_break_selection.tree.optimized_branch_lengths import (
    fit_fixed_topology_nnls_branch_lengths,
)


@pytest.mark.parametrize("diffusion", [False, True])
@pytest.mark.parametrize("energy", [False, True])
def test_factorial_nnls_matches_existing_runner(diffusion, energy):
    _, _, _, _, data, labels, space, distances = next(datasets(0, 1))
    if diffusion:
        distances = feature_family_diffusion_geometry(data, space).distance_condensed
    tree = build_tree(data, distances, linkage_method="average").tree
    fit = fit_fixed_topology_nnls_branch_lengths(
        tree,
        data,
        pair_sample_size=50000,
        random_state=0,
        solver_tolerance=1e-5,
        max_iterations=1000,
    )
    assert fit.applied_to_tree
    result, annotations = evaluate(tree, data, labels, space, nnls=True, energy=energy)
    reference = run_tbs_on_distance(
        data,
        distances,
        0.01,
        feature_space=space,
        adaptive_projection_dimension_energy_fraction=0.9 if energy else None,
        branch_length_optimization_method="fixed_topology_nnls",
        branch_length_data_df=data,
        branch_length_optimization_pair_sample_size=50000,
        branch_length_optimization_solver_tolerance=1e-5,
        branch_length_optimization_max_iterations=1000,
        edge_branch_length_variance_policy="normalized_branch_length",
    )
    assert result["status"] == reference.status == "ok"
    assert result["ari"] == adjusted_rand_score(labels, reference.labels)
    expected = reference.extra["annotations"].loc[annotations.index]
    for column in (
        "Child_Parent_Divergence_Test_Statistic",
        "Sibling_Test_Statistic",
        "Sibling_Divergence_P_Value",
        "Sibling_Degrees_of_Freedom",
    ):
        np.testing.assert_allclose(annotations[column], expected[column], equal_nan=True)
