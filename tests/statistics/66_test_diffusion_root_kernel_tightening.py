"""Avoid interval dependence on shared translations and irrelevant neighbor order."""

from decimal import Decimal
from fractions import Fraction

import numpy as np
import pytest

from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_kernel import (
    _i,
    _neighbor_support,
    _quadratic_range,
    _radial_squared_distances,
)
from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_replay import (
    make_diffusion_root_path,
)


def test_convex_quadratic_bounds_include_interior_minimum_and_endpoint_maximum():
    # (r-2)^2 has a zero at an interior point, not at either endpoint.
    lower, upper = _quadratic_range(Fraction(4), _i(-4), Fraction(1), 1, 4)
    assert lower == 0 and upper == 4
    lower, upper = _quadratic_range(Fraction(4), (Decimal(-5), Decimal(-3)), Fraction(1), 1, 4)
    assert lower <= Fraction(-9, 4) and upper >= 8


def test_neighbor_support_does_not_require_order_inside_selected_set():
    row = [_i(0), (Decimal(1), Decimal(3)), (Decimal(2), Decimal(4)), _i(8)]
    assert set(_neighbor_support(row, self_index=0, counts=(2,))[:2]) == {1, 2}
    with pytest.raises(ValueError, match="support"):
        _neighbor_support(row, self_index=0, counts=(1,))


def test_within_child_distances_remain_exact_across_wide_radial_interval():
    path = make_diffusion_root_path(np.random.default_rng(20261105).normal(size=(9, 1)))
    narrow, _ = _radial_squared_distances(path, path.observed_radius, path.observed_radius)
    wide, _ = _radial_squared_distances(path, 0, 10)
    groups = {i: [i] for i in range(len(path.values))}
    for step, (a, b) in enumerate(path.observed_signature):
        groups[len(path.values) + step] = groups[a] + groups[b]
    for child in path.observed_signature[-1]:
        members = groups[child]
        for i in members:
            for j in members:
                assert wide[i][j] == narrow[i][j]
                actual = sum(
                    (Fraction(float(a)) - Fraction(float(b))) ** 2
                    for a, b in zip(path.values[i], path.values[j])
                )
                assert Fraction(wide[i][j][0]) <= actual <= Fraction(wide[i][j][1])
    for radius in (0, 1, path.observed_radius, 10):
        observed = path.offset + radius * path.slope
        for i in range(len(observed)):
            for j in range(len(observed)):
                squared = float(np.sum((observed[i] - observed[j]) ** 2))
                assert float(wide[i][j][0]) - 1e-12 <= squared <= float(wide[i][j][1]) + 1e-12
