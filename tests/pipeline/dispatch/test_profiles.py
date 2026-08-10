from __future__ import annotations

import pytest
from benchmarks.shared.runners.dispatch import run_clustering_result
from benchmarks.shared.runners.method_registry import METHOD_SPECS
from benchmarks.shared.types.method_spec import MethodSpec
from scipy.spatial.distance import pdist

from .helpers import _capturing_runner, _toy_dataframe


@pytest.mark.parametrize(
    ("method_id", "method_name", "profile_id"),
    [
        pytest.param(
            "tbs_conditional_topology_diagnostic",
            "TBS (Conditional Topology Diagnostic)",
            "fixed_coordinate_conditional_topology_diagnostic_v1",
            id="conditional-topology",
        ),
        pytest.param(
            "tbs_global_passthrough_refined_diagnostic",
            "TBS (Global Passthrough Refined Diagnostic)",
            "fixed_coordinate_global_passthrough_refined_v1",
            id="global-passthrough-refined",
        ),
        pytest.param(
            "tbs_spectral_transport_passthrough",
            "TBS (Spectral Transport Passthrough)",
            "fixed_coordinate_spectral_transport_passthrough_v1",
            id="spectral-transport",
        ),
    ],
)
def test_run_clustering_result_dispatches_profile_as_kl(
    monkeypatch,
    method_id: str,
    method_name: str,
    profile_id: str,
) -> None:
    captured: dict[str, object] = {}
    monkeypatch.setitem(
        METHOD_SPECS,
        method_id,
        MethodSpec(
            name=method_name,
            runner=_capturing_runner(captured),
            param_grid=[
                {
                    "tree_distance_metric": "hamming",
                    "tree_linkage_method": "average",
                    "sibling_gate_profile": profile_id,
                }
            ],
        ),
    )

    data = _toy_dataframe()
    run_clustering_result(
        data_df=data,
        method_id=method_id,
        params={
            "tree_distance_metric": "euclidean",
            "tree_linkage_method": "average",
            "sibling_gate_profile": profile_id,
        },
        seed=42,
        distance_condensed=pdist(data.values, metric="euclidean"),
    )

    assert captured["kwargs"]["sibling_gate_profile"] == profile_id
    assert captured["kwargs"]["tree_linkage_method"] == "average"
