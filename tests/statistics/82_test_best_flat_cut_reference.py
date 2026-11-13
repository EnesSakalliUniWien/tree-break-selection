"""Contracts for the opt-in Gaussian best-flat-cut reference sibling p-value."""

from __future__ import annotations

import numpy as np
import pytest

from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence import (
    best_flat_cut_reference as module,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.best_flat_cut_reference import (
    best_flat_cut_reference,
    load_table,
    reference_p_value,
    split_rho,
)

P = 30


def _direct_rho(z: np.ndarray, in_left: np.ndarray) -> float:
    left, right = z[in_left], z[~in_left]
    n, pi = len(z), in_left.mean()
    pooled = (
        (left - left.mean(0)).T @ (left - left.mean(0))
        + (right - right.mean(0)).T @ (right - right.mean(0))
    ) / (n - 2)
    delta = left.mean(0) - right.mean(0)
    return float(pi * (1 - pi) * delta @ np.linalg.solve(pooled, delta))


def test_split_rho_is_pi_one_minus_pi_times_pooled_mahalanobis_squared() -> None:
    rng = np.random.default_rng(0)
    z = rng.standard_normal((40, 3))
    in_left = rng.random(40) < 0.4
    assert np.isclose(split_rho(z, in_left), _direct_rho(z, in_left))


def test_split_rho_is_invariant_under_invertible_linear_map_and_shift() -> None:
    rng = np.random.default_rng(1)
    z = rng.standard_normal((60, 4))
    in_left = z[:, 0] + 0.5 * z[:, 1] > 0.2
    transform = rng.standard_normal((4, 4)) + 3 * np.eye(4)
    assert np.isclose(split_rho(z, in_left), split_rho(z @ transform.T + 7.0, in_left))


def test_perfect_separation_gives_infinite_rho() -> None:
    z = np.c_[np.r_[np.zeros(5), np.ones(5)], np.zeros(10)]
    z[:, 1] = np.tile([0.0, 1.0], 5)  # within-side covariance singular in coordinate 0
    assert np.isinf(split_rho(z, np.arange(10) < 5))


def test_table_cells_are_sorted_and_finite_except_small_n_above_10_directions() -> None:
    table = load_table()
    n_grid, k_grid, upper = table["n_grid"], table["k_grid"], table["upper"]
    assert upper.shape == (len(n_grid), len(k_grid), int(np.ceil(0.1 * int(table["B"]))))
    assert n_grid.min() >= module.MIN_CELLS
    assert k_grid.tolist() == list(range(1, module.MAX_DIRECTIONS + 1))
    not_built = (k_grid[None, :] > 10) & (n_grid[:, None] < module.MIN_CELLS_ABOVE_10_DIRECTIONS)
    assert np.isnan(upper[not_built]).all()
    assert np.isfinite(upper[~not_built]).all()
    assert (np.diff(upper[~not_built], axis=-1) >= 0).all()


def test_p_value_is_monte_carlo_rank_against_stored_values() -> None:
    table = load_table()
    n, k = 100, 3
    upper = table["upper"][list(table["n_grid"]).index(n), k - 1].astype(np.float64)
    bound = int(table["B"]) + 1
    for rho in (upper[0], upper[50], upper[-1], upper[-1] * 2):
        p_value, used = reference_p_value(float(rho), n, k)
        assert used == n
        assert p_value == (1 + np.sum(upper >= rho)) / bound
    assert reference_p_value(float(upper[0] * 0.99), n, k) == (1.0, n)
    assert reference_p_value(float("inf"), n, k)[0] == 1 / bound


def test_n_between_grid_values_uses_the_lower_grid_n() -> None:
    rho = float(load_table()["upper"][:, 2, 100].max())
    assert reference_p_value(rho, 120, 3)[1] == 100
    assert reference_p_value(rho, 149, 3)[1] == 100
    assert reference_p_value(rho, 150, 3)[1] == 150
    assert reference_p_value(rho, 15, 3)[1] == 15
    assert reference_p_value(rho, 5000, 3)[1] == 800


def test_fail_closed_reasons() -> None:
    noise = np.random.default_rng(0).standard_normal((400, P))
    result = best_flat_cut_reference(noise, noise[:, 0] > 0)
    assert (result.reason, result.p_value, result.invariant_directions) == (
        "no_invariant_directions",
        1.0,
        0,
    )

    rng = np.random.default_rng(2)
    wide = rng.standard_normal((600, 40))
    wide[:, :18] *= 6.0
    result = best_flat_cut_reference(wide, wide[:, 0] > 0)
    assert result.reason == "too_many_invariant_directions" and result.p_value == 1.0
    assert result.invariant_directions > module.MAX_DIRECTIONS

    small = rng.standard_normal((12, 3))
    small[:6, 0] += 8.0
    result = best_flat_cut_reference(small, np.arange(12) < 6)
    assert result.reason == "too_few_cells" and result.p_value == 1.0
    assert result.invariant_directions >= 1


def test_split_with_11_to_15_invariant_directions_is_tested() -> None:
    rng = np.random.default_rng(2)
    wide = rng.standard_normal((600, P))
    wide[:, :12] *= 6.0
    result = best_flat_cut_reference(wide, wide[:, 0] > 0)
    assert result.reason == "tested" and result.reference_n == 400
    assert 10 < result.invariant_directions <= module.MAX_DIRECTIONS
    assert 1 / (int(load_table()["B"]) + 1) <= result.p_value <= 1.0


def test_above_10_directions_and_fewer_than_50_cells_fails_closed() -> None:
    rng = np.random.default_rng(5)
    narrow = rng.standard_normal((40, P))
    narrow[:, :12] *= 6.0
    result = best_flat_cut_reference(narrow, narrow[:, 0] > 0)
    assert result.invariant_directions > 10
    assert (result.reason, result.p_value, result.reference_n) == ("too_few_cells", 1.0, 0)
    with pytest.raises(ValueError, match="No reference"):
        reference_p_value(1.0, 40, 12)


def test_tested_split_of_two_separated_groups_has_the_smallest_p_value() -> None:
    rng = np.random.default_rng(3)
    shift = np.r_[8.0, np.zeros(P - 1)]
    points = np.vstack(
        [rng.standard_normal((150, P)) + shift, rng.standard_normal((150, P)) - shift]
    )
    result = best_flat_cut_reference(points, np.arange(300) < 150)
    assert result.reason == "tested" and result.reference_n == 300
    assert result.p_value == 1 / (int(load_table()["B"]) + 1)


def test_one_sided_input_raises() -> None:
    with pytest.raises(ValueError, match="both sides"):
        best_flat_cut_reference(np.zeros((4, P)), np.ones(4, dtype=bool))
    with pytest.raises(ValueError, match="both sides"):
        split_rho(np.zeros((4, 2)), np.zeros(4, dtype=bool))
