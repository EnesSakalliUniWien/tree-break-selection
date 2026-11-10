import numpy as np
import pandas as pd
import pytest

from benchmarks.shared.metrics import calculate_clustering_metrics, matched_group_recovery
from benchmarks.shared.statistical_profiles import build_statistical_profiles


def test_one_to_one_fragmentation_differs_from_legacy_majority_remap():
    truth = np.array(["A", "A", "A", "A", "B", "B"])
    report = pd.DataFrame({"cluster_id": [1, 1, 2, 2, 3, 3]}, index=list("abcdef"))
    metrics = calculate_clustering_metrics(report.iloc[::-1], report.index, truth)
    assert metrics.macro_f1 == 1
    assert metrics.matched_macro_f1 == pytest.approx(5 / 6)
    assert metrics.matched_rare_f1 == 1
    assert metrics.pair_split_share == pytest.approx(4 / 7)
    assert metrics.pair_merge_share == 0
    with pytest.raises(ValueError, match="Unmapped"):
        calculate_clustering_metrics(report.rename(index={"a": "z"}), report.index, truth)


def test_matching_aligns_series_and_leaves_unmatched_truth_groups_at_zero():
    truth = pd.Series(["A", "A", "B"], index=list("abc"))
    predicted = pd.Series([7, 7, 7], index=list("cba"))
    matched = matched_group_recovery(truth, predicted).set_index("group")
    assert matched.loc["A", "f1"] == pytest.approx(0.8)
    assert matched.loc["B", "f1"] == 0 and matched.loc["B", "rare"]
    assert pd.isna(matched.loc["B", "predicted"])
    metrics = calculate_clustering_metrics(
        predicted.to_frame("cluster_id"), truth.index, truth.to_numpy()
    )
    assert metrics.pair_split_share == 0
    assert metrics.pair_merge_share == pytest.approx(2 / 3)
    assert metrics.matched_rare_f1 == 0
    exact = pd.Series([8, 7, 7], index=list("cba"))
    assert matched_group_recovery(truth, exact).f1.tolist() == [1, 1]
    with pytest.raises(ValueError, match="index"):
        matched_group_recovery(truth, predicted.rename(index={"a": "z"}))


def test_profile_separation_and_spectra_are_rotation_and_translation_invariant():
    local = np.array([[-1, 0], [1, 0], [0, -2], [0, 2]])
    points = np.r_[local, local + [5, 0]]
    index = list("abcdefgh")
    labels = pd.Series(["a"] * 4 + ["b"] * 4, index=index).iloc[::-1]
    subtrees = {"s": ([0, 1, 2, 3], [4, 5, 6, 7])}
    profile = build_statistical_profiles(
        pd.DataFrame(points, index=index),
        true_labels=labels,
        subtrees=subtrees,
        geometry_disagreement={"s": 0.4},
    )
    rotated = build_statistical_profiles(
        points @ np.array([[0.6, -0.8], [0.8, 0.6]]) + [9, 4],
        true_labels=labels.reindex(index).to_numpy(),
        subtrees=subtrees,
    )
    group = profile.groups.set_index("group").loc["a"]
    assert group["size"] == 4
    np.testing.assert_allclose(group.covariance_eigenvalues, [8 / 3, 2 / 3])
    assert group.covariance_rank == 2 and 1 < group.effective_rank < 2
    pair = profile.pairs.iloc[0]
    assert pair.centroid_distance == 5
    assert pair.mahalanobis_distance == pytest.approx(5 / np.sqrt(2 / 3))
    assert pair.centroid_overlap_share == 0
    assert pair.effective_n == 2 and pair.smaller_share == 0.5
    for column in ["centroid_distance", "mahalanobis_distance", "centroid_overlap_share"]:
        assert profile.pairs[column].iloc[0] == pytest.approx(rotated.pairs[column].iloc[0])
        assert profile.pairs[column].iloc[0] == pytest.approx(profile.subtrees[column].iloc[0])
    assert profile.subtrees.geometry_disagreement.iloc[0] == 0.4
    assert np.isnan(rotated.subtrees.geometry_disagreement.iloc[0])


def test_profiles_keep_singular_and_unavailable_fields_explicit():
    points = np.array([[0, 0], [1, 0], [2, 1], [3, 1]])
    profile = build_statistical_profiles(points, subtrees={"s": ([0, 1], [2, 3])})
    assert profile.groups.empty
    node = profile.subtrees.iloc[0]
    assert node.pooled_covariance_rank == 1
    assert np.isnan(node.mahalanobis_distance)
    assert node.mahalanobis_status == "rank_deficient"
    with pytest.raises(ValueError, match="finite"):
        build_statistical_profiles([[0, 0], [0, np.nan]])
    with pytest.raises(ValueError, match="disjoint"):
        build_statistical_profiles(points, subtrees={"s": ([0, 1], [1, 2])})
    with pytest.raises(ValueError, match="geometry"):
        build_statistical_profiles(points, geometry_disagreement={"missing": 0.4})


def test_centroid_overlap_gives_equidistant_rows_half_weight_after_rotation():
    points = np.array([[-1, 0], [1, 0], [1, 0], [3, 0]])
    rotated = points @ np.array([[0.6, -0.8], [0.8, 0.6]]) + [10, 15]
    for x in (points, rotated):
        node = build_statistical_profiles(x, subtrees={"s": ([0, 1], [2, 3])}).subtrees.iloc[0]
        assert node.centroid_overlap_share == 0.25


@pytest.mark.parametrize(
    ("names", "truth", "report_index", "match"),
    [
        (list("abc"), [0, 0], list("abc"), "true_labels"),
        (list("abc"), [0, 0, 1], list("ab"), "report_df has"),
        (list("aab"), [0, 0, 1], list("aab"), "Duplicate"),
    ],
)
def test_calculate_clustering_metrics_rejects_misaligned_input(names, truth, report_index, match):
    report = pd.DataFrame({"cluster_id": range(len(report_index))}, index=report_index)
    with pytest.raises(ValueError, match=match):
        calculate_clustering_metrics(report, pd.Index(names), np.array(truth))
