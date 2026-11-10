"""Subtree-only candidate search, selected Ward events and retained sibling scores."""

import numpy as np
import pandas as pd
import pytest
from scipy.spatial.distance import pdist

from benchmarks.diagnostics.calibration.selected.hierarchy.subtree_gaussian_reference import (
    subtree_gaussian_chain_test,
)
from benchmarks.diagnostics.calibration.selected.hierarchy.subtree_topology_geometry import Cut
from benchmarks.diagnostics.calibration.selected.hierarchy.subtree_topology_tests import (
    native_sibling_statistic,
    sibling_cut_scores,
    subtree_alternative_cut_test,
    subtree_topology_search_test,
    topology_pair_reference,
    topology_sibling_score,
    ward_subtree_selection,
)
from benchmarks.shared.runners.tbs_runner import run_tbs_on_distance
from tree_break_selection.tree.feature_space import continuous_feature_space_from_columns


def test_sibling_score_matches_pooled_hotelling_and_parent_fallback():
    x = np.random.default_rng(12).normal(size=(20, 3))
    for count in (4, 10):
        left, right = x[:count], x[count:]
        if count < 8:
            covariance = np.cov(x, rowvar=False)
        else:
            residual = np.vstack([left - left.mean(0), right - right.mean(0)])
            covariance = residual.T @ residual / 18
        difference = left.mean(0) - right.mean(0)
        expected = (
            difference @ np.linalg.solve(covariance, difference) / (1 / count + 1 / (20 - count))
        )
        score = sibling_cut_scores(x, [Cut(tuple(range(count)), tuple(range(count, 20)))])
        assert score[0] == pytest.approx(expected, rel=1e-10)


def test_joint_search_adjusts_every_returned_cut_using_the_same_maximum_path():
    x = np.random.default_rng(83).normal(size=(5, 2))
    result = subtree_alternative_cut_test(x, subtree_indices=range(5), steps=19, seed=32, alpha=0.1)
    assert result["status"] == "complete"
    assert result["candidate_family"] == "full_span_strict_affine_cuts"
    assert result["cut_results"]
    assert max(row["statistic"] for row in result["cut_results"]) == result["observed_score"]
    for row in result["cut_results"]:
        assert row["global_null_adjusted_p_value"] >= result["p_value"]
        assert row["global_null_adjusted_p_value"] == np.mean(
            np.array(result["chain_scores"]) >= row["statistic"]
        )
    assert result == subtree_alternative_cut_test(
        x, subtree_indices=range(5), steps=19, seed=32, alpha=0.1
    )


def test_ward_descendant_event_freezes_ancestral_memberships():
    x = np.random.default_rng(3).normal(size=(12, 2))
    x[6:] += 30
    rows, event = ward_subtree_selection(x, path=(0,))
    assert np.array_equal(rows, ward_subtree_selection(x, path=iter([0]))[0])
    assert rows.tolist() == list(range(6)) and event(x)
    swapped = x.copy()
    swapped[[0, 8]] = swapped[[8, 0]]
    assert not event(swapped)

    def score(matrix):
        assert event(matrix)
        np.testing.assert_array_equal(matrix[6:], x[6:])
        return float(matrix[0, 0])

    result = subtree_gaussian_chain_test(
        x, subtree_indices=rows, selection_event=event, score=score, steps=19, seed=10, alpha=0.1
    )
    assert result["status"] == "complete" and result["accepted_moves"] > 0
    assert result["subtree_indices"] == list(range(6))


def test_preserving_node_split_is_stricter_than_its_ancestry():
    x = np.array([[0.0], [1.0], [8.0], [10.0], [100.0], [101.0], [110.0], [111.0]])
    rows, event = ward_subtree_selection(x, path=(0,))
    _, pair_event = ward_subtree_selection(x, path=(0,), keep_split=True)
    changed = x.copy()
    changed[[1, 2]] = changed[[2, 1]]
    assert rows.tolist() == [0, 1, 2, 3]
    assert event(changed) and not pair_event(changed)


def test_native_score_reproduces_the_actual_fitted_raw_root_statistic():
    x = np.random.default_rng(9).normal(size=(16, 3))
    x[8:, 0] += 5
    data = pd.DataFrame(x, index=[f"S{i:012d}" for i in range(len(x))])
    fitted = run_tbs_on_distance(
        data,
        pdist(x),
        0.05,
        tree_linkage_method="ward",
        feature_space=continuous_feature_space_from_columns(data.columns),
        adaptive_projection_dimension_energy_fraction=0.90,
        branch_length_optimization_method="fixed_topology_nnls",
        branch_length_data_df=data,
    )
    annotated = fitted.extra["annotations"].loc[fitted.extra["tree"].root()]
    # The displayed statistic has already been divided by empirical inflation.
    expected = annotated["Sibling_Test_Statistic"] * annotated["Sibling_Inflation_Factor"]
    assert np.isfinite(expected)
    assert native_sibling_statistic(x, range(16)) == pytest.approx(expected)


def test_paired_topology_reference_uses_a_difference_on_each_same_data_replay():
    x = np.random.default_rng(52).normal(size=(4, 1))
    first, second = ((0, 1), (2, 3)), ((0, 2), (1, 3))
    result = topology_pair_reference(
        x, subtree_indices=range(4), first=first, second=second, steps=19, seed=1, alpha=0.1
    )
    assert result["status"] == "complete"
    assert result["observed_score"] == pytest.approx(
        abs(topology_sibling_score(x, first) - topology_sibling_score(x, second))
    )
    same = topology_pair_reference(
        x, subtree_indices=range(4), first=first, second=first, steps=19, seed=1, alpha=0.1
    )
    assert same["p_value"] == 1 and same["observed_score"] == 0


def test_invalid_partition_and_duplicate_tree_leaf_cannot_be_scored():
    x = np.arange(8).reshape(4, 2)
    with pytest.raises(ValueError, match="partition"):
        sibling_cut_scores(x, [Cut((0, 1), (1, 2, 3))])
    with pytest.raises(ValueError, match="each local row"):
        topology_sibling_score(x, ((0, 1), (1, 3)))
    with pytest.raises(ValueError, match="child positions"):
        ward_subtree_selection(x, path=(2,))
    with pytest.raises(ValueError, match="child positions"):
        ward_subtree_selection(x, path=(0.0,))


def test_joint_topology_reference_counts_ties_and_replays_the_search():
    x = np.array([[-1.0, -1.0], [-1.0, 1.0], [1.0, -1.0], [1.0, 1.0]])
    result = subtree_topology_search_test(x, subtree_indices=range(4), steps=9, seed=40, alpha=0.1)
    assert result["status"] == "complete"
    assert len(result["topology_results"]) == 2
    assert result["observed_score"] == max(r["statistic"] for r in result["topology_results"])
    assert all(
        r["global_null_adjusted_p_value"] >= result["p_value"] for r in result["topology_results"]
    )


def test_observed_candidate_failure_remains_unresolved(monkeypatch):
    from benchmarks.diagnostics.calibration.selected.hierarchy import subtree_topology_tests as mod

    def fail(_):
        raise RuntimeError("LP failed")

    monkeypatch.setattr(mod, "_candidate_scores", fail)
    result = mod.subtree_alternative_cut_test(
        np.arange(8).reshape(4, 2), subtree_indices=range(4), steps=9, seed=1, alpha=0.1
    )
    assert result["status"] == "unresolved" and result["p_value"] is None
    assert result["error_phase"] == "observation" and result["cut_results"] == []
