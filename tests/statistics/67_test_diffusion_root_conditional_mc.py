"""Finite rank accounting and native diffusion membership controls."""

from types import SimpleNamespace

import numpy as np
import pytest

from benchmarks.diagnostics.calibration.selected.hierarchy import (
    diffusion_root_conditional_mc as mc,
)


def fake_path(statuses, observed_radius=1.0):
    states = iter(statuses)
    return SimpleNamespace(
        dimension=2,
        observed_radius=observed_radius,
        observed_signature=((0, 1),),
        source_hashes={},
        evaluate=lambda radius: SimpleNamespace(
            status=next(states),
            error="diagnostic",
            signature=((0, 1),),
            margins=np.array([0.0]),
            state={"merge_tie_tolerance": 1e-12},
        ),
    )


def test_rank_includes_observation_and_uses_only_selected_radii(monkeypatch):
    statuses = ["selected", "nonselected", "selected", "selected"]
    monkeypatch.setattr(mc, "make_diffusion_root_path", lambda *a, **kw: fake_path(statuses))
    result = mc.selected_diffusion_root_mc_test([[0], [1]], sigma=1, proposals=4, seed=7)
    radii = np.sqrt(np.random.default_rng(7).chisquare(2, size=4))[[0, 2, 3]]
    np.testing.assert_array_equal(result["accepted_radii"], radii)
    assert result["p_value"] == (1 + np.count_nonzero(radii >= 1)) / 4
    assert result["minimum_p_value"] == 0.25
    assert result["proposals_evaluated"] == 4


def test_zero_acceptance_is_conservative(monkeypatch):
    monkeypatch.setattr(
        mc, "make_diffusion_root_path", lambda *a, **kw: fake_path(["nonselected"] * 3)
    )
    result = mc.selected_diffusion_root_mc_test([[0], [1]], sigma=1, proposals=3)
    assert result["accepted_count"] == 0
    assert result["p_value"] == result["minimum_p_value"] == 1


def test_strict_history_violation_resolves_a_later_tie(monkeypatch):
    path = fake_path([])
    path.evaluate = lambda radius: SimpleNamespace(
        status="tie",
        error=None,
        signature=((0, 2),),
        margins=np.array([-1.0, 0.0]),
        state={"merge_tie_tolerance": 1e-12},
    )
    monkeypatch.setattr(mc, "make_diffusion_root_path", lambda *a, **kw: path)
    result = mc.selected_diffusion_root_mc_test([[0], [1]], sigma=1, proposals=3)
    assert result["status"] == "complete" and result["p_value"] == 1
    assert result["strict_exclusions_with_ties"] == 3


@pytest.mark.parametrize("status", ["tie", "unsupported"])
def test_unresolved_membership_cannot_be_discarded(monkeypatch, status):
    monkeypatch.setattr(
        mc, "make_diffusion_root_path", lambda *a, **kw: fake_path(["selected", status])
    )
    result = mc.selected_diffusion_root_mc_test([[0], [1]], sigma=1, proposals=3)
    assert result["status"] == "unresolved"
    assert result["p_value"] is None and result["proposals_evaluated"] == 2


@pytest.mark.parametrize("proposals", [0, -1, True, 1.5])
def test_invalid_proposal_budget(proposals):
    with pytest.raises(ValueError, match="positive integer"):
        mc.selected_diffusion_root_mc_test([[0], [1]], sigma=1, proposals=proposals)


def test_native_replay_is_reproducible_and_preserves_observed_data():
    values = np.random.default_rng(20261105).normal(size=(9, 1))
    original = values.copy()
    first = mc.selected_diffusion_root_mc_test(values, sigma=1, proposals=30, seed=8)
    second = mc.selected_diffusion_root_mc_test(values, sigma=1, proposals=30, seed=8)
    assert first == second
    assert first["status"] == "complete"
    assert 0 < first["p_value"] <= 1
    np.testing.assert_array_equal(values, original)


def test_native_first_null_case_excludes_history_despite_later_near_tie():
    values = np.random.default_rng(2026110401).normal(size=(9, 3))
    result = mc.selected_diffusion_root_mc_test(
        values,
        sigma=1,
        proposals=3,
        seed=2026120401,
    )
    assert result["status"] == "complete"
    assert result["strict_exclusions_with_ties"] >= 1
