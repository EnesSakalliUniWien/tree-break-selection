from __future__ import annotations

import numpy as np
from scipy.spatial.distance import pdist, squareform

from benchmarks.shared.runners.dispatch import run_clustering_result
from benchmarks.shared.runners.method_registry import METHOD_SPECS
from benchmarks.shared.types.method_spec import MethodSpec
from tree_break_selection.tree.continuous_distance import (
    continuous_time_distance_condensed,
    standardized_euclidean_distance_condensed,
)
from tree_break_selection.tree.feature_space import continuous_feature_space_from_columns

from .helpers import _capturing_runner, _toy_dataframe


def test_run_clustering_result_uses_provided_tbs_distance_condensed(monkeypatch):
    captured = {}
    params = {"tree_distance_metric": "euclidean", "tree_linkage_method": "average"}
    monkeypatch.setitem(
        METHOD_SPECS,
        "tbs",
        MethodSpec(name="TBS Divergence", runner=_capturing_runner(captured), param_grid=[params]),
    )
    df = _toy_dataframe()
    # Cityblock differs from the Euclidean distance dispatch would compute itself.
    supplied = pdist(df.values, metric="cityblock")

    run_clustering_result(
        data_df=df,
        method_id="tbs",
        params=params,
        seed=42,
        distance_condensed=supplied,
    )

    np.testing.assert_array_equal(captured["args"][1], supplied)


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


def test_run_clustering_result_uses_provided_graph_distance_matrix(monkeypatch):
    captured = {}
    params = {"n_neighbors": 2, "resolution": 1.0}
    monkeypatch.setitem(
        METHOD_SPECS,
        "leiden",
        MethodSpec(name="Leiden", runner=_capturing_runner(captured), param_grid=[params]),
    )
    df = _toy_dataframe()
    # Dispatch would otherwise build a Hamming matrix from the data.
    supplied = squareform(pdist(df.values, metric="cityblock"))

    run_clustering_result(
        data_df=df,
        method_id="leiden",
        params=params,
        seed=42,
        distance_matrix=supplied,
    )

    np.testing.assert_array_equal(captured["args"][0], supplied)


def test_run_clustering_result_gives_distance_tree_nnls_methods_their_own_geometry(
    monkeypatch,
):
    """Distance-tree NNLS methods must get their own constructor that supplies
    branch_length_data_df explicitly, mirroring tbs_diffusion_graphtools_nnls,
    rather than falling through the generic TBS_DISTANCE_TREE_METHODS dispatch branch
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
