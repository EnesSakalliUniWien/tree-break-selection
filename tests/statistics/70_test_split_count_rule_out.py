"""Contracts for the opt-in split-count rule-out sibling p-value."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.stats import chi2

from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.split_count_rule_out import (
    MAX_EQUIDISTANT_SHARE,
    equidistant_share,
    invariant_directions,
    is_separable,
    log_split_count,
    split_count_rule_out,
)

P = 30


def _groups(rng, centers, per):
    points = np.vstack([c + rng.standard_normal((per, P)) for c in centers])
    labels = np.repeat(np.arange(len(centers)), per)
    return points, labels


def test_log_split_count_counts_two_group_splits() -> None:
    assert np.isclose(np.exp(log_split_count(5)), 15)
    assert log_split_count(2) == 0.0
    assert np.isclose(log_split_count(400), 399 * np.log(2.0))


def test_invariant_directions_follow_the_noise_band() -> None:
    rng = np.random.default_rng(0)
    assert invariant_directions(rng.standard_normal((400, P))).shape == (0, P)
    centers = np.c_[np.eye(4) - 0.25, np.zeros((4, P - 4))] * 8.0
    points, _ = _groups(rng, centers, 250)
    assert invariant_directions(points).shape == (3, P)
    assert invariant_directions(rng.standard_normal((P + 1, P))).shape == (0, P)


def test_separability_and_equidistant_cells() -> None:
    rng = np.random.default_rng(1)
    shift = np.r_[6.0, np.zeros(P - 1)]
    points, labels = _groups(rng, [shift, -shift], 200)
    in_left = labels == 0
    assert is_separable(points, in_left)
    assert equidistant_share(points, in_left) == 0.0
    noise = rng.standard_normal((400, P))
    cut = noise[:, 0] > 0
    assert not is_separable(noise[:, 1:4], cut)
    assert equidistant_share(noise, cut) > MAX_EQUIDISTANT_SHARE


def test_rule_out_keeps_real_splits_and_rules_out_noise() -> None:
    rng = np.random.default_rng(2)
    shift = np.r_[6.0, np.zeros(P - 1)]
    points, labels = _groups(rng, [shift, -shift], 200)
    kept = split_count_rule_out(5000.0, 30.0, points, labels == 0)
    assert kept.reason == "kept"
    expected = np.exp(log_split_count(400) + chi2.logsf(5000.0, 30.0))
    assert np.isclose(kept.p_value, expected)

    noise = rng.standard_normal((400, P))
    ruled_out = split_count_rule_out(5000.0, 30.0, noise, noise[:, 0] > 0)
    assert ruled_out.reason == "no_invariant_directions"
    assert ruled_out.p_value == 1.0

    few = np.r_[0:16, 200:215]  # P + 1 cells from both groups
    small = split_count_rule_out(5000.0, 30.0, points[few], labels[few] == 0)
    assert small.reason == "no_invariant_directions"


def test_rule_out_rejects_one_sided_input() -> None:
    with pytest.raises(ValueError, match="both sides"):
        split_count_rule_out(1.0, 1.0, np.zeros((4, P)), np.ones(4, dtype=bool))
