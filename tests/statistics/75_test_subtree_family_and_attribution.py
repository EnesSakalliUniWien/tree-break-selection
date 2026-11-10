"""Common-event Ward family and native score attribution checks."""

import numpy as np
import pytest

from benchmarks.diagnostics.calibration.selected.hierarchy.subtree_family_reference import (
    ward_subtree_family_test,
)
from benchmarks.diagnostics.calibration.selected.hierarchy.subtree_statistic_attribution import (
    native_sibling_attribution,
)
from benchmarks.diagnostics.calibration.selected.hierarchy.subtree_topology_tests import (
    native_sibling_statistic,
)


def test_family_replays_ancestry_and_full_selected_subtree():
    x = np.array([0.0, 1.0, 7.0, 9.0, 60.0, 61.0, 68.0, 70.0])[:, None]

    def score(matrix, rows, left, right):
        outside = sorted(set(range(len(x))) - set(rows))
        np.testing.assert_array_equal(matrix[outside], x[outside])
        return abs(matrix[left, 0].mean() - matrix[right, 0].mean())

    r = ward_subtree_family_test(x, path=iter([0]), score=score, steps=39, seed=17)
    assert r["status"] == "complete"
    assert r["subtree_indices"] == list(range(4))
    assert r["hypotheses"] == 3
    assert {tuple(n["subtree_indices"]) for n in r["nodes"]} == {(0, 1), (2, 3), (0, 1, 2, 3)}
    assert all(n["selection_event"] == "caller_supplied" for n in r["nodes"])
    assert all(n["holm_p_value"] >= (n["p_value"] or 0) for n in r["nodes"])


def test_family_retains_failures_and_marks_insufficient_resolution():
    x = np.random.default_rng(4).normal(size=(8, 2))

    def score(matrix, rows, left, right):
        if len(rows) == 2:
            raise FloatingPointError("unsupported local score")
        return matrix[left, 0].mean() - matrix[right, 0].mean()

    r = ward_subtree_family_test(x, score=score, steps=99, seed=1)
    assert r["hypotheses"] == 7 and not r["resolution_sufficient"]
    assert any(n["p_value"] is None and n["holm_p_value"] == 1 for n in r["nodes"])
    assert all(n["family_decision"] == "unresolved" for n in r["nodes"])


def test_attribution_retains_native_score_and_changes_only_branch_factor():
    x = np.random.default_rng(9).normal(size=(16, 3))
    x[8:, 0] += 5
    a = native_sibling_attribution(x, range(16), oracle_covariance=np.eye(3), gram_rank=1)
    assert a["native_raw"] == pytest.approx(native_sibling_statistic(x, range(16)))
    assert a["unit_branch_raw"] == pytest.approx(a["branch_factor"] * a["native_raw"])
    assert a["native_df"] == a["unit_branch_df"]
    assert a["native_display"] == pytest.approx(a["native_raw"] / (a["scale"] * a["inflation"]))
    assert np.isfinite(a["oracle_covariance_raw"]) and np.isfinite(a["gram_rank_raw"])
    assert a["gram_rank_df"] == 1


def test_holm_keeps_failed_members_and_shares_one_event(monkeypatch):
    from benchmarks.diagnostics.calibration.selected.hierarchy import (
        subtree_family_reference as mod,
    )

    callbacks = []

    def chain(matrix, *, selection_event, seed, **kwargs):
        callbacks.append(selection_event)
        assert selection_event(matrix)
        p_value = (0.01, None, 0.02)[seed]
        return dict(p_value=p_value, decision="unresolved" if p_value is None else "reject")

    monkeypatch.setattr(mod, "subtree_gaussian_chain_test", chain)
    x = np.array([0.0, 1.0, 8.0, 10.0])[:, None]
    r = mod.ward_subtree_family_test(x)
    assert callbacks[0] is callbacks[1] is callbacks[2]
    assert [n["holm_p_value"] for n in r["nodes"]] == [0.03, 1, 0.04]
    assert [n["family_decision"] for n in r["nodes"]] == ["reject", "unresolved", "reject"]
    assert r["unresolved_nodes"] == 1 and r["decision"] == "reject"
    assert r["minimum_steps_for_holm"] == 59


def test_observed_ward_numerical_failure_is_unresolved():
    x = np.array([[1e200], [0.0], [-1e200]])
    r = ward_subtree_family_test(x)
    assert r["decision"] == r["status"] == "unresolved"
    assert r["error_phase"] == "observation" and not r["nodes"]
