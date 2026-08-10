from __future__ import annotations

import numpy as np
from benchmarks.shared.runners.dispatch import run_clustering_result
from benchmarks.shared.runners.method_registry import METHOD_SPECS
from benchmarks.shared.types.method_spec import MethodSpec
from scipy.spatial.distance import pdist, squareform
from tree_break_selection.tree.continuous_distance import (
    continuous_time_distance_condensed,
    standardized_euclidean_distance_condensed,
)
from tree_break_selection.tree.feature_space import continuous_feature_space_from_columns

from .helpers import _capturing_runner, _toy_dataframe


def test_run_clustering_result_uses_provided_tbs_distance_condensed():
    df = _toy_dataframe()
    dist_condensed = pdist(df.values, metric="euclidean")
    result = run_clustering_result(
        data_df=df,
        method_id="tbs",
        params={"tree_distance_metric": "euclidean", "tree_linkage_method": "average"},
        seed=42,
        distance_condensed=dist_condensed,
    )

    assert result.status == "ok"
    assert result.labels is not None
    assert len(result.labels) == len(df)
    assert result.skip_reason is None


def test_run_clustering_result_builds_continuous_mahalanobis_time_distance(monkeypatch):
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
                    "tree_distance_metric": "mahalanobis_time",
                    "tree_linkage_method": "average",
                }
            ],
        ),
    )

    df = _toy_dataframe()
    feature_space = continuous_feature_space_from_columns(tuple(df.columns))
    run_clustering_result(
        data_df=df,
        method_id="tbs",
        params={"tree_distance_metric": "mahalanobis_time", "tree_linkage_method": "average"},
        seed=42,
        feature_space=feature_space,
    )

    np.testing.assert_allclose(
        captured["args"][1],
        continuous_time_distance_condensed(df.values, feature_space),
    )


def test_run_clustering_result_builds_continuous_standardized_euclidean_distance(
    monkeypatch,
):
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
                    "tree_distance_metric": "standardized_euclidean",
                    "tree_linkage_method": "average",
                }
            ],
        ),
    )

    df = _toy_dataframe()
    feature_space = continuous_feature_space_from_columns(tuple(df.columns))
    run_clustering_result(
        data_df=df,
        method_id="tbs",
        params={
            "tree_distance_metric": "standardized_euclidean",
            "tree_linkage_method": "average",
        },
        seed=42,
        feature_space=feature_space,
    )

    np.testing.assert_allclose(
        captured["args"][1],
        standardized_euclidean_distance_condensed(df.values, feature_space),
    )


def test_run_clustering_result_dispatches_iqtree_without_condensed_distance(monkeypatch):
    captured = {}
    capture_runner = _capturing_runner(captured)

    monkeypatch.setitem(
        METHOD_SPECS,
        "tbs_iqtree3",
        MethodSpec(
            name="TBS (IQ-TREE 3, MAD Root)",
            runner=capture_runner,
            param_grid=[
                {
                    "tree_distance_metric": "hamming",
                    "tree_linkage_method": "average",
                    "tree_builder": "iqtree3",
                    "tree_rooting": "mad",
                }
            ],
        ),
    )

    run_clustering_result(
        data_df=_toy_dataframe(),
        method_id="tbs_iqtree3",
        params={
            "tree_distance_metric": "hamming",
            "tree_linkage_method": "average",
            "tree_builder": "iqtree3",
            "tree_rooting": "mad",
        },
        seed=42,
        distance_condensed=pdist(_toy_dataframe().values, metric="euclidean"),
    )

    assert captured["args"][1] is None
    assert captured["kwargs"]["tree_builder"] == "iqtree3"
    assert captured["kwargs"]["tree_rooting"] == "mad"


def test_run_clustering_result_uses_provided_graph_distance_matrix():
    df = _toy_dataframe()
    dist_condensed = pdist(df.values, metric="euclidean")
    dist_matrix = squareform(dist_condensed)
    result = run_clustering_result(
        data_df=df,
        method_id="leiden",
        params={"n_neighbors": 2, "resolution": 1.0},
        seed=42,
        distance_matrix=dist_matrix,
    )

    assert result.status in {"ok", "skip"}
    if result.status == "ok":
        assert result.labels is not None
        assert len(result.labels) == len(df)
        assert result.skip_reason is None
    else:
        assert result.labels is None
        assert isinstance(result.skip_reason, str)
        assert result.skip_reason.strip()


def test_run_clustering_result_gives_distance_tree_nnls_methods_their_own_geometry(
    monkeypatch,
):
    """Distance-tree NNLS methods must get their own constructor that supplies
    branch_length_data_df explicitly, mirroring tbs_diffusion_graphtools_nnls,
    rather than falling through the generic TBS_RUNNER_METHODS dispatch branch
    that never wires branch_length_data_df through."""
    captured = {}
    capture_runner = _capturing_runner(captured, include_args=False)

    monkeypatch.setitem(
        METHOD_SPECS,
        "tbs_nnls",
        MethodSpec(
            name="TBS Divergence (NNLS branch lengths)",
            runner=capture_runner,
            param_grid=[{}],
        ),
    )

    df = _toy_dataframe()
    run_clustering_result(
        data_df=df,
        method_id="tbs_nnls",
        params={
            "tree_distance_metric": "hamming",
            "tree_linkage_method": "average",
            "branch_length_optimization_method": "fixed_topology_nnls",
        },
        seed=42,
    )

    assert captured["kwargs"]["branch_length_data_df"] is df
