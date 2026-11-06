"""Exact spectral bounds; matrix controls alone do not certify native KDE branches."""

from fractions import Fraction

import numpy as np
import pytest
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist

from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_cells import (
    _classify_linkage,
    _distance_enclosures,
    _spectral_reference,
    _sqrt_down,
    _sqrt_up,
)
from tree_break_selection.space_separation.diffusion import compute_diffusion_coordinates


@pytest.mark.parametrize(
    "value", [Fraction(0), Fraction(2), Fraction(1, 10**100), Fraction(10**70, 3)]
)
def test_square_root_bounds_use_exact_integer_arithmetic(value):
    lower, upper = _sqrt_down(value), _sqrt_up(value)
    assert lower * lower <= value <= upper * upper
    assert lower <= upper


def _matrix_control():
    rng = np.random.default_rng(20261109)
    weights = rng.uniform(0.01, 0.2, (9, 9))
    weights[:4, :4] += 0.8
    weights[4:, 4:] += 0.8
    weights = np.maximum(weights, weights.T)
    np.fill_diagonal(weights, 0.0)
    degree = weights.sum(1)
    operator = weights / np.sqrt(degree[:, None] * degree[None, :])
    coordinates = compute_diffusion_coordinates(weights, diffusion_time=3, n_components=1)
    signature = tuple(
        tuple(sorted(map(int, row[:2]))) for row in linkage(pdist(coordinates), method="average")
    )
    return operator, coordinates, signature


def test_exact_point_spectral_enclosure_certifies_a_full_history_control():
    operator, coordinates, signature = _matrix_control()
    exact = [[Fraction(float(v)) for v in row] for row in operator]
    reference = _spectral_reference(operator, target_dimension=1)
    lower, upper, details = _distance_enclosures(reference, exact, exact)
    status, margins = _classify_linkage(lower, upper, signature, 1e-12)
    assert status == "selected"
    assert Fraction(details["factor_error"]) > 0
    assert Fraction(margins["minimum_margin_lower"]) > 0
    from scipy.spatial.distance import squareform

    distances = squareform(pdist(coordinates))
    for i in range(9):
        for j in range(9):
            assert float(lower[i][j]) - 1e-14 <= distances[i, j] <= float(upper[i][j]) + 1e-14


def test_nonzero_operator_box_bounds_distances_and_a_competing_history():
    operator, _, signature = _matrix_control()
    reference = _spectral_reference(operator, target_dimension=1)
    width = Fraction(1, 10**12)
    lower_operator = [[Fraction(float(v)) - width for v in row] for row in operator]
    upper_operator = [[Fraction(float(v)) + width for v in row] for row in operator]
    lower, upper, _ = _distance_enclosures(reference, lower_operator, upper_operator)
    status, _ = _classify_linkage(lower, upper, signature, 1e-12)
    assert status == "selected"
    competing = list(signature)
    a, b = competing[0]
    other = next(i for i in range(9) if i not in (a, b))
    competing[0] = tuple(sorted((a, other)))
    status, _ = _classify_linkage(lower, upper, tuple(competing), 1e-12)
    assert status == "nonselected"


def test_spectral_branch_uncertainty_is_not_classified():
    operator, _, _ = _matrix_control()
    reference = _spectral_reference(operator, target_dimension=1)
    lower = [[Fraction(float(v)) - Fraction(1, 2) for v in row] for row in operator]
    upper = [[Fraction(float(v)) + Fraction(1, 2) for v in row] for row in operator]
    with pytest.raises(ValueError, match="separated"):
        _distance_enclosures(reference, lower, upper)


def test_linkage_tangency_and_ties_remain_unresolved():
    lower = [
        [Fraction(0), Fraction(1), Fraction(1)],
        [Fraction(1), Fraction(0), Fraction(2)],
        [Fraction(1), Fraction(2), Fraction(0)],
    ]
    status, _ = _classify_linkage(lower, lower, ((0, 1), (2, 3)), 1e-12)
    assert status == "unresolved"


def test_actual_exact_real_diffusion_has_certified_reentry_windows_and_unknown_tail():
    from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_cells import (
        analyze_diffusion_root_cells,
    )
    from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_mass import (
        bounded_root_gate,
    )
    from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_replay import (
        make_diffusion_root_path,
    )

    path = make_diffusion_root_path(np.random.default_rng(20261105).normal(size=(9, 1)))
    analysis = analyze_diffusion_root_cells(
        path, observed_half_width=1e-12, anchor_ratios=(1.0, 1.1, 1.2)
    )
    assert analysis.selection_certified
    assert analysis.cells[0].lo == 0.0
    assert analysis.cells[-1].hi == np.inf
    assert analysis.cells[-1].status == "unresolved"
    assert analysis.cells[-1].proof == "uncertified_unbounded_tail"
    for first, second in zip(analysis.cells[:-1], analysis.cells[1:]):
        assert first.hi == second.lo
    for ratio, status in [(1.0, "selected"), (1.1, "nonselected"), (1.2, "selected")]:
        radius = path.observed_radius * ratio
        cell = next(c for c in analysis.cells if c.lo < radius < c.hi)
        assert cell.status == status
        assert cell.proof == "outward_native_kernel_and_exact_spectral_linkage_bounds"
        assert cell.details["kernel"]["exact_gaussian_path"]
        assert Fraction(cell.details["spectral"]["spectral_gap_lower"]) > 0
        if status == "selected":
            assert Fraction(cell.details["linkage"]["minimum_margin_lower"]) > 0
    rlo, rhi = analysis.observed_radius_bounds
    observed_cell = next(c for c in analysis.cells if c.lo < path.observed_radius < c.hi)
    assert Fraction(observed_cell.lo) < rlo <= rhi < Fraction(observed_cell.hi)
    bounds, gate = bounded_root_gate(analysis, path, 0.01)
    assert bounds.selection_mass_lower > 0
    assert bounds.selection_validated
    assert gate.status == "unresolved"
    assert gate.clusters == 1


def test_exact_gaussian_radius_bounds_do_not_assume_float_path_is_unit():
    from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_cells import (
        _exact_observed_radius_bounds,
    )
    from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_kernel import (
        exact_path_radius_bounds,
    )
    from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_replay import (
        make_diffusion_root_path,
    )

    path = make_diffusion_root_path(np.random.default_rng(20261105).normal(size=(9, 1)))
    integer_lower, integer_upper = _exact_observed_radius_bounds(path)
    decimal_lower, decimal_upper = map(Fraction, exact_path_radius_bounds(path))
    assert integer_lower <= decimal_upper and decimal_lower <= integer_upper


def test_native_dyadic_tie_tolerance_is_not_rounded_down_to_decimal_rational():
    from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_cells import (
        _SpectralReference,
    )

    exact_decimal = Fraction(1, 10**10) + Fraction(1, 10**12)
    native_dyadic = Fraction(1e-10) + Fraction(1e-12)
    assert native_dyadic > exact_decimal
    gap = (exact_decimal + native_dyadic) / 2
    eigenvalues = (Fraction(1), 1 - gap, Fraction(1, 10))
    operator = tuple(
        tuple(value if i == j else Fraction(0) for j, value in enumerate(eigenvalues))
        for i in range(3)
    )
    gram = tuple(
        tuple(eigenvalues[1] ** 6 if i == j == 1 else Fraction(0) for j in range(3))
        for i in range(3)
    )
    reference = _SpectralReference(eigenvalues, operator, gram, Fraction(0), Fraction(1), 1, 2, 3)
    with pytest.raises(ValueError, match="uniformly separated"):
        _distance_enclosures(reference, operator, operator)


def test_arithmetic_enclosure_failure_stays_unresolved(monkeypatch):
    from decimal import Overflow

    from benchmarks.diagnostics.calibration.selected.hierarchy import diffusion_root_kernel
    from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_cells import (
        analyze_diffusion_root_cells,
    )
    from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_replay import (
        make_diffusion_root_path,
    )

    path = make_diffusion_root_path(np.random.default_rng(20261105).normal(size=(9, 1)))

    def overflow(*args, **kwargs):
        raise Overflow("kernel interval overflow")

    monkeypatch.setattr(diffusion_root_kernel, "enclose_diffusion_operator", overflow)
    analysis = analyze_diffusion_root_cells(path)
    assert not analysis.selection_certified
    assert all(cell.status == "unresolved" for cell in analysis.cells)
    assert all("overflow" in cell.details["reason"] for cell in analysis.cells[:-1])


def test_actual_nontrivial_width_uses_linear_factor_distance_certificate():
    from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_cells import (
        analyze_diffusion_root_cells,
    )
    from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_replay import (
        make_diffusion_root_path,
    )

    path = make_diffusion_root_path(np.random.default_rng(20261105).normal(size=(9, 1)))
    analysis = analyze_diffusion_root_cells(path, observed_half_width=1e-8)
    observed = next(c for c in analysis.cells if c.lo < path.observed_radius < c.hi)
    assert observed.status == "selected"
    assert observed.details["spectral"]["distance_bound_method"] == "diffusion_factor_linear_error"
    assert analysis.selection_certified


def test_bounded_expansion_increases_certified_mass_without_dropping_unknown_regions():
    from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_cells import (
        analyze_diffusion_root_cells,
    )
    from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_mass import (
        bounded_root_gate,
    )
    from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_replay import (
        make_diffusion_root_path,
    )

    path = make_diffusion_root_path(np.random.default_rng(20261105).normal(size=(9, 1)))
    base = analyze_diffusion_root_cells(path, observed_half_width=1e-8)
    grown = analyze_diffusion_root_cells(path, observed_half_width=1e-8, expansion_steps=3)
    original = next(c for c in base.cells if c.lo < path.observed_radius < c.hi)
    expanded = next(c for c in grown.cells if c.lo < path.observed_radius < c.hi)
    assert expanded.status == "selected"
    assert expanded.lo < original.lo and expanded.hi > original.hi
    assert 0 < expanded.details["validated_expansion_steps"] <= 3
    assert grown.cells[0].lo == 0 and grown.cells[-1].hi == np.inf
    assert grown.cells[-1].status == "unresolved"
    assert all(a.hi == b.lo for a, b in zip(grown.cells[:-1], grown.cells[1:]))
    before, _ = bounded_root_gate(base, path, 0.01)
    after, gate = bounded_root_gate(grown, path, 0.01)
    assert after.selection_mass_lower > before.selection_mass_lower > 0
    assert gate.status == "unresolved" and gate.clusters == 1


def test_selected_spectral_bound_certifies_wider_actual_window():
    from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_cells import (
        analyze_diffusion_root_cells,
    )
    from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_replay import (
        make_diffusion_root_path,
    )

    path = make_diffusion_root_path(np.random.default_rng(20261105).normal(size=(9, 1)))
    analysis = analyze_diffusion_root_cells(path, observed_half_width=1e-6)
    observed = next(c for c in analysis.cells if c.lo < path.observed_radius < c.hi)
    assert observed.status == "selected"
    assert observed.details["spectral"]["factor_bound_method"] == "selected_spectrum_telescope"
    assert analysis.selection_certified
