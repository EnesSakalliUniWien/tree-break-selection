from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist

from benchmarks.shared.util.decomposition import (
    labels_and_report_from_decomposition,
    ok_result_from_labels,
    tbs_decompose,
)


def test_labels_and_report_from_decomposition_share_sample_assignment_table() -> None:
    decomposition = {
        "cluster_assignments": {
            1: {
                "root_node": "left",
                "leaves": ["S1", "S3"],
                "leaf_signature": ("S1", "S3"),
                "size": 2,
            },
            2: {"root_node": "right", "leaves": ["S2"], "leaf_signature": ("S2",), "size": 1},
        }
    }

    labels, report_df = labels_and_report_from_decomposition(
        decomposition,
        ["S3", "S0", "S2", "S1"],
    )

    assert np.array_equal(labels, np.array([1, -1, 2, 1], dtype=int))
    assert report_df.to_dict(orient="index") == {
        "S1": {"cluster_id": 1, "cluster_size": 2},
        "S3": {"cluster_id": 1, "cluster_size": 2},
        "S2": {"cluster_id": 2, "cluster_size": 1},
    }
    assert report_df.index.name == "sample_id"


def test_ok_result_from_labels_counts_non_noise_clusters() -> None:
    result = ok_result_from_labels(
        np.array([2, 2, -1, 4], dtype=int),
        pd.Index(["a", "b", "noise", "c"], name="sample_id"),
    )

    assert result.status == "ok"
    assert result.skip_reason is None
    assert result.found_clusters == 2
    assert np.array_equal(result.labels, np.array([2, 2, -1, 4], dtype=int))
    assert result.report_df is not None
    assert result.report_df.loc["noise", "cluster_size"] == 1


def test_tbs_decompose_returns_labels_in_data_row_order() -> None:
    rng = np.random.default_rng(0)
    centers = rng.uniform(size=(2, 30)) > 0.5
    rows = np.vstack(
        [np.where(rng.uniform(size=(20, 30)) < 0.05, ~center, center) for center in centers]
    ).astype(float)
    # Names out of natural order, so build_tree links samples in a different order.
    names = [f"cell_{index}" for index in rng.permutation(40)]
    data = pd.DataFrame(rows, index=names)

    labels, decomposition = tbs_decompose(
        data, pdist(rows, "hamming"), linkage_method="average", alpha=0.05
    )

    cluster_by_leaf = {
        leaf: cluster_id
        for cluster_id, info in decomposition["cluster_assignments"].items()
        for leaf in info["leaves"]
    }
    assert labels.tolist() == [cluster_by_leaf[name] for name in names]
