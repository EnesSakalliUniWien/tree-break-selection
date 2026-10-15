"""Discovery-only diffusion/NNLS integration and frozen inference contracts."""

from __future__ import annotations

import json
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from tree_break_selection.hierarchy_analysis.statistics import independent_gaussian as api


def paired_data():
    discovery_seed, inference_seed = np.random.SeedSequence(19382).spawn(2)
    labels = [f"L{i}" for i in range(24)]
    discovery = pd.DataFrame(
        np.random.default_rng(discovery_seed).normal(size=(24, 4)), index=labels
    )
    inference = pd.DataFrame(
        np.random.default_rng(inference_seed).normal(size=(24, 4)), index=labels
    )
    return discovery, inference


def fit(discovery, *, nnls=True):
    return api.fit_independent_gaussian_tree(
        discovery,
        known_covariance=np.eye(4),
        observation_ids=[f"d{i}" for i in range(24)],
        diffusion=api.IndependentGaussianDiffusion(fit_nnls=nnls),
    )


def test_diffusion_and_nnls_consume_only_discovery_and_are_not_refit():
    discovery, inference = paired_data()
    with (
        patch.object(
            api, "adaptive_diffusion_geometry", wraps=api.adaptive_diffusion_geometry
        ) as diffusion,
        patch.object(
            api,
            "fit_fixed_topology_nnls_branch_lengths",
            wraps=api.fit_fixed_topology_nnls_branch_lengths,
        ) as nnls,
    ):
        model = fit(discovery)
        pd.testing.assert_frame_equal(diffusion.call_args.args[0], discovery)
        pd.testing.assert_frame_equal(nnls.call_args.args[1], discovery)
        assert diffusion.call_count == nnls.call_count == 1
    assert not model.diffusion_coordinates.flags.writeable
    geometry = json.loads(model.geometry_metadata_json)
    branch = json.loads(model.branch_length_metadata_json)
    assert geometry["metric"] == "euclidean"
    assert branch["status"] == "ok" and branch["applied_to_tree"]
    assert branch["target_metric"] == "squared_standardized_euclidean"
    before = tuple(model.tree_branch_lengths)
    with (
        patch.object(
            api,
            "adaptive_diffusion_geometry",
            side_effect=AssertionError("test data reached diffusion"),
        ),
        patch.object(
            api,
            "fit_fixed_topology_nnls_branch_lengths",
            side_effect=AssertionError("test data reached NNLS"),
        ),
        patch.object(
            api, "eigendecompose_covariance", side_effect=AssertionError("test data reached PCA")
        ),
    ):
        result = api.test_independent_gaussian_tree(
            model, inference, observation_ids=[f"t{i}" for i in range(24)]
        )
    assert before == model.tree_branch_lengths
    assert result.tests.p_value.between(0, 1).all()
    assert len(result.decomposition["traversal_trace"]) > 0


def test_nnls_changes_geometry_but_not_known_covariance_pvalues():
    discovery, inference = paired_data()
    plain, optimized = fit(discovery, nnls=False), fit(discovery, nnls=True)
    assert plain.tree_edges == optimized.tree_edges
    assert plain.tree_branch_lengths != optimized.tree_branch_lengths
    a = api.test_independent_gaussian_tree(
        plain, inference, observation_ids=[f"t{i}" for i in range(24)]
    )
    b = api.test_independent_gaussian_tree(
        optimized, inference, observation_ids=[f"t{i}" for i in range(24)]
    )
    pd.testing.assert_frame_equal(a.tests, b.tests)
    assert a.decomposition["num_clusters"] == b.decomposition["num_clusters"]


def test_changed_fitted_branch_length_is_rejected():
    discovery, inference = paired_data()
    model = fit(discovery)
    u, v, length = model.tree_branch_lengths[0]
    model.tree.edges[u, v]["branch_length"] = length + 1
    with pytest.raises(ValueError, match="branch lengths changed"):
        api.test_independent_gaussian_tree(
            model, inference, observation_ids=[f"t{i}" for i in range(24)]
        )


def test_changed_diffusion_coordinates_are_rejected():
    discovery, inference = paired_data()
    model = fit(discovery)
    # A caller deliberately overrides array write protection; checksum still guards testing.
    model.diffusion_coordinates.setflags(write=True)
    model.diffusion_coordinates[0, 0] += 1
    with pytest.raises(ValueError, match="diffusion coordinates changed"):
        api.test_independent_gaussian_tree(
            model, inference, observation_ids=[f"t{i}" for i in range(24)]
        )


def test_nnls_nonconvergence_is_not_silently_accepted():
    discovery, _ = paired_data()
    failed = SimpleNamespace(
        status="solver_not_converged",
        applied_to_tree=False,
        to_dict=lambda: {"status": "solver_not_converged"},
    )
    with patch.object(api, "fit_fixed_topology_nnls_branch_lengths", return_value=failed):
        with pytest.raises(ValueError, match="NNLS did not converge"):
            fit(discovery)
