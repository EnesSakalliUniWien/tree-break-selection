"""Real-arithmetic kernel enclosure checks against native point observations."""

from fractions import Fraction

import numpy as np
import pytest

from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_kernel import (
    enclose_diffusion_operator,
    exact_path_radius_bounds,
)
from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_replay import (
    make_diffusion_root_path,
)


def test_actual_reentry_operator_has_positive_width_validated_enclosure():
    path = make_diffusion_root_path(np.random.default_rng(20261105).normal(size=(9, 1)))
    r = path.observed_radius
    enclosure = enclose_diffusion_operator(path, r - 1e-10, r + 1e-10)
    assert enclosure.retained_count_target == 1 and enclosure.time == 3
    assert all(isinstance(x, Fraction) for row in enclosure.lower for x in row)
    assert enclosure.proof_details["kde_dimension"] == 1
    assert (
        enclosure.proof_details["kernel_contract"]
        == "exact-real formulas; not bitwise NumPy enclosure"
    )
    assert float(enclosure.observed_radius_bounds[0]) <= r + 1e-14
    lo, hi = exact_path_radius_bounds(path)
    assert lo <= hi
    for fraction in [-0.9, 0, 0.9]:
        point = path.evaluate(r + fraction * 1e-10)
        # Recover native symmetric normalized weights from its captured actual kernel distances.
        d = point.state["scaled_distances"].toarray()
        mask = point.state["scaled_distances"].astype(bool).toarray()
        w = np.where(mask, np.exp(-d * d / (4 * point.state["outer_epsilon"])), 0.0)
        np.fill_diagonal(w, 0.0)
        w = np.maximum(w, w.T)
        degree = w.sum(axis=1)
        operator = w / np.sqrt(np.outer(degree, degree))
        for i in range(9):
            for j in range(9):
                # Positive-width enclosure dwarfs native point rounding; agreement isn't its proof.
                assert float(enclosure.lower[i][j]) <= operator[i, j] + 1e-14
                assert operator[i, j] <= float(enclosure.upper[i][j]) + 1e-14


def test_wide_branch_ambiguity_remains_unresolved():
    path = make_diffusion_root_path(np.random.default_rng(20261105).normal(size=(9, 1)))
    with pytest.raises(ValueError, match="unresolved"):
        enclose_diffusion_operator(path, 0.0, 10.0)


def test_signed_interval_arithmetic_and_exact_radius_enclose_rational_target():
    from decimal import Decimal, localcontext

    from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_kernel import (
        _div,
        _exp,
        _i,
        _ln,
        _mul,
    )

    assert _mul((_i(-3)[0], _i(-2)[1]), (_i(2)[0], _i(4)[1])) == (Decimal(-12), Decimal(-4))
    assert _div((_i(-3)[0], _i(-2)[1]), (_i(2)[0], _i(4)[1])) == (Decimal("-1.5"), Decimal("-.5"))
    for primitive, value, operation in [(_exp, Fraction(1, 3), "exp"), (_ln, Fraction(7, 3), "ln")]:
        interval = primitive(_i(value))
        with localcontext() as context:
            context.prec = 110
            reference = getattr(Decimal(value.numerator) / Decimal(value.denominator), operation)()
        assert interval[0] <= reference <= interval[1]
    path = make_diffusion_root_path(np.random.default_rng(20261105).normal(size=(9, 1)))
    groups = {i: [i] for i in range(9)}
    for step, (a, b) in enumerate(path.observed_signature):
        groups[9 + step] = groups[a] + groups[b]
    a, b = path.observed_signature[-1]
    eta = [
        Fraction(1, len(groups[a])) if i in groups[a] else -Fraction(1, len(groups[b]))
        for i in range(9)
    ]
    v = sum(eta[i] * Fraction(float(path.values[i, 0])) for i in range(9))
    squared = v * v / sum(e * e for e in eta)
    lo, hi = map(Fraction, exact_path_radius_bounds(path))
    assert lo * lo <= squared <= hi * hi


def test_actual_eighty_row_attribution_case_is_not_artificially_unsupported():
    from benchmarks.diagnostics.calibration.diffusion_attribution import gaussian_tree_pairs

    pair = next(gaussian_tree_pairs(replicates=1))
    path = make_diffusion_root_path(pair.selected_data)
    assert path.values.shape == (80, 24)
    assert path.observed_radius == pytest.approx(8.765613678802739)
    r = path.observed_radius
    enclosure = enclose_diffusion_operator(path, r - 1e-12, r + 1e-12)
    assert enclosure.retained_count_target == 24
    assert enclosure.proof_details["neighbor_k"] == 15
    assert enclosure.proof_details["kde_bgh_index"] == 37
    assert enclosure.proof_details["kde_dimension"] == 2
    assert len(enclosure.lower) == 80
    for lower, upper in zip(enclosure.lower, enclosure.upper):
        assert all(left <= right for left, right in zip(lower, upper))
    point = path.evaluate(r)
    scaled = point.state["scaled_distances"]
    distances = scaled.toarray()
    diagonal = scaled.diagonal()
    assert point.state["positive_self_scaled_distance_count"] == int(
        np.count_nonzero(diagonal > 0.0)
    )
    assert point.state["zero_self_scaled_distance_count"] == int(np.count_nonzero(diagonal == 0.0))
    assert point.state["positive_scaled_distance_count"] == int(np.count_nonzero(scaled.data > 0.0))
    assert point.state["actual_positive_scaled_sq_median"] == np.median(
        scaled.data[scaled.data > 0.0] ** 2
    )
    actual_weights = np.where(
        scaled.astype(bool).toarray(),
        np.exp(-(distances**2) / (4 * point.state["outer_epsilon"])),
        0.0,
    )
    np.fill_diagonal(actual_weights, 0.0)
    actual_weights = np.maximum(actual_weights, actual_weights.T)
    actual_degrees = actual_weights.sum(axis=1)
    actual_operator = actual_weights / np.sqrt(np.outer(actual_degrees, actual_degrees))
    # The repository boundary now reconciles native self identity with the exact-real map.
    np.fill_diagonal(distances, 0.0)
    exact_real_epsilon = np.median(distances[distances > 0.0] ** 2)
    assert point.state["off_diagonal_positive_scaled_sq_median"] == exact_real_epsilon
    weights = np.where(
        distances > 0.0,
        np.exp(-(distances**2) / (4 * exact_real_epsilon)),
        0.0,
    )
    np.fill_diagonal(weights, 0.0)
    weights = np.maximum(weights, weights.T)
    degree = weights.sum(axis=1)
    native_formula_with_zero_self = weights / np.sqrt(np.outer(degree, degree))
    lower = np.array(enclosure.lower, dtype=float)
    upper = np.array(enclosure.upper, dtype=float)
    # Numerical agreement only; exact-real containment comes from interval operations.
    assert np.all(native_formula_with_zero_self >= lower - 1e-14)
    assert np.all(native_formula_with_zero_self <= upper + 1e-14)

    assert point.state["positive_self_scaled_distance_count"] == 0
    assert point.state["zero_self_scaled_distance_count"] == 80
    assert point.state["outer_epsilon"] == exact_real_epsilon
    assert (
        point.state["metadata"]["self_distance_policy"]
        == "zero_identity_diagonal_before_epsilon_and_kernel"
    )
    assert np.all(actual_operator >= lower - 1e-14)
    assert np.all(actual_operator <= upper + 1e-14)
