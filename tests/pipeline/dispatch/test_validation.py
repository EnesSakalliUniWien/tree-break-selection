from __future__ import annotations

import pytest
from benchmarks.shared.runners.dispatch import (
    _tbs_branch_length_optimization_kwargs,
    run_clustering_result,
)
from benchmarks.shared.runners.method_registry import METHOD_SPECS
from benchmarks.shared.types.method_spec import MethodSpec

from .helpers import _capturing_runner, _toy_dataframe


def test_dispatch_result_rejects_invalid_runner_params():
    df = _toy_dataframe()
    with pytest.raises(ValueError):
        run_clustering_result(
            data_df=df,
            method_id="kmeans",
            params={"n_clusters": "bad"},
            seed=42,
        )


def test_dispatch_result_rejects_invalid_spectral_params():
    df = _toy_dataframe()
    with pytest.raises(ValueError):
        run_clustering_result(
            data_df=df,
            method_id="spectral",
            params={"n_clusters": "bad"},
            seed=42,
        )


def test_dispatch_result_reraises_unexpected_exception(monkeypatch):
    def _raise_runner(*_args, **_kwargs):
        raise RuntimeError("boom")

    monkeypatch.setitem(
        METHOD_SPECS,
        "kmeans",
        MethodSpec(
            name="K-Means",
            runner=_raise_runner,
            param_grid=[{"n_clusters": 2}],
        ),
    )

    df = _toy_dataframe()
    with pytest.raises(RuntimeError, match="boom"):
        run_clustering_result(
            data_df=df,
            method_id="kmeans",
            params={"n_clusters": 2},
            seed=42,
        )


def test_branch_length_optimization_apply_nonconverged_param_is_strict_bool():
    assert (
        _tbs_branch_length_optimization_kwargs(
            {"branch_length_optimization_apply_nonconverged": "false"}
        )["branch_length_optimization_apply_nonconverged"]
        is False
    )
    assert (
        _tbs_branch_length_optimization_kwargs(
            {"branch_length_optimization_apply_nonconverged": "true"}
        )["branch_length_optimization_apply_nonconverged"]
        is True
    )

    with pytest.raises(ValueError, match="branch_length_optimization_apply_nonconverged"):
        _tbs_branch_length_optimization_kwargs(
            {"branch_length_optimization_apply_nonconverged": "maybe"}
        )


def test_tbs_dispatch_bool_params_do_not_treat_false_strings_as_true(monkeypatch):
    captured = {}
    capture_runner = _capturing_runner(captured, include_args=False)

    monkeypatch.setitem(
        METHOD_SPECS,
        "tbs",
        MethodSpec(
            name="TBS",
            runner=capture_runner,
            param_grid=[{}],
        ),
    )
    df = _toy_dataframe()

    run_clustering_result(
        data_df=df,
        method_id="tbs",
        params={
            "tree_distance_metric": "euclidean",
            "tree_linkage_method": "average",
            "spectral_include_internal_barycenters": "false",
            "enforce_internal_support_thresholds": "false",
            "allow_linkage_ultrametric_branch_time": "false",
            "passthrough": "false",
        },
        seed=42,
    )

    for key in (
        "spectral_include_internal_barycenters",
        "enforce_internal_support_thresholds",
        "allow_linkage_ultrametric_branch_time",
        "passthrough",
    ):
        assert captured["kwargs"][key] is False
