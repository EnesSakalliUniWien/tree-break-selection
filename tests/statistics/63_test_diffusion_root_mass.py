"""Certified probability arithmetic; synthetic cells are mathematical controls only."""

from dataclasses import replace
from decimal import Decimal, localcontext
from fractions import Fraction
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.stats import chi

from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_mass import (
    bounded_root_gate,
    chi_interval_mass,
    conditional_p_bounds,
    root_gate,
)


def cell(lo, hi, status):
    return SimpleNamespace(lo=lo, hi=hi, status=status)


def test_unrestricted_and_tiny_chi_two_mass():
    assert chi_interval_mass(0, np.inf, 2) == (Decimal(1), Decimal(1))
    lo, hi = chi_interval_mass(0, 1e-30, 2)
    with localcontext() as context:
        context.prec = 150
        x = Decimal.from_float(1e-30) ** 2 / 2
        reference = 1 - (-x).exp()
    assert 0 < lo <= reference <= hi


@pytest.mark.parametrize("dimension", [1, 2, 3, 24, 25])
def test_probability_mass_matches_independent_scipy_control(dimension):
    for a, b in ((0, 1), (1, 3), (3, 9), (20, np.inf)):
        lo, hi = chi_interval_mass(a, b, dimension)
        approximate = chi.sf(a, dimension) - chi.sf(b, dimension)
        assert float(lo) - 1e-12 <= approximate <= float(hi) + 1e-12
        assert 0 <= lo <= hi <= 1


def test_disconnected_selected_union_and_unknown_tail_bound():
    cells = [
        cell(0, 1, "selected"),
        cell(1, 2, "nonselected"),
        cell(2, 3, "selected"),
        cell(3, np.inf, "unresolved"),
    ]
    bounds = conditional_p_bounds(cells, 2.5, 2)
    with localcontext() as context:
        context.prec = 120

        def sf(r):
            return (-(Decimal(str(r)) ** 2) / 2).exp()

        below = 1 - sf(1) + sf(2) - sf(2.5)
        above = sf(2.5) - sf(3)
        expected_lower = above / (above + below)
        expected_upper = (above + sf(3)) / (above + sf(3) + below)
    assert bounds.lower <= expected_lower <= bounds.upper
    assert bounds.lower <= expected_upper <= bounds.upper
    assert bounds.selection_mass_lower > 0


def test_missing_or_unsupported_regions_are_not_nonselection():
    for cells in ([cell(0, 2, "unresolved")], [cell(0, np.inf, "unsupported")]):
        bounds = conditional_p_bounds(cells, 1, 2)
        assert (bounds.lower, bounds.upper) == (Decimal(0), Decimal(1))
        assert root_gate(bounds, 0.01).clusters == 1
        assert root_gate(bounds, 0.01).status == "unresolved"


def test_gate_requires_selection_certificate_and_does_not_resolve_straddle():
    bounds = conditional_p_bounds([cell(0, np.inf, "selected")], 5, 2)
    assert bounds.upper < Decimal("0.01")
    assert root_gate(bounds, 0.01).status == "unresolved"
    # This flag represents an audited classifier in the mathematical gate control.
    certified = replace(bounds, selection_validated=True)
    assert root_gate(certified, 0.01).status == "reject"
    assert root_gate(certified, 0.01).clusters == 2
    ordinary = replace(
        conditional_p_bounds([cell(0, np.inf, "selected")], 1, 2), selection_validated=True
    )
    assert root_gate(ordinary, 0.01).status == "nonreject"
    straddle = replace(certified, lower=Decimal(0), upper=Decimal(1))
    assert root_gate(straddle, 0.01).status == "unresolved"


def test_overlapping_cells_are_rejected():
    with pytest.raises(ValueError, match="overlap"):
        conditional_p_bounds([cell(0, 2, "selected"), cell(1, 3, "selected")], 1, 2)


def test_unknown_below_cannot_be_ignored_in_conditional_denominator():
    bounds = conditional_p_bounds([cell(0, 1, "unresolved"), cell(1, np.inf, "selected")], 2, 2)
    with localcontext() as context:
        context.prec = 120
        lower = Decimal(-2).exp()
        upper = Decimal("-1.5").exp()
    assert bounds.lower <= lower <= bounds.upper
    assert bounds.lower <= upper <= bounds.upper


def test_far_tail_stays_positive_below_float_underflow():
    lower, upper = chi_interval_mass(1000, np.inf, 24)
    assert 0 < lower <= upper < Decimal("1e-200000")
    assert (upper - lower) / lower < Decimal("1e-60")


def test_unverified_scanner_labels_cannot_create_a_pipeline_pvalue():
    analysis = SimpleNamespace(cells=[cell(0, np.inf, "selected")])
    path = SimpleNamespace(observed_radius=5, dimension=2)
    bounds, gate = bounded_root_gate(analysis, path, 0.01)
    assert (bounds.lower, bounds.upper) == (Decimal(0), Decimal(1))
    assert gate.status == "unresolved" and gate.clusters == 1 and not gate.passthrough


def test_gate_threshold_boundary():
    bounds = conditional_p_bounds([cell(0, np.inf, "selected")], 1, 2)
    at_threshold = replace(
        bounds, lower=Decimal("0.01"), upper=Decimal("0.01"), selection_validated=True
    )
    assert root_gate(at_threshold, 0.01).status == "reject"
    straddle = replace(at_threshold, upper=Decimal("0.02"))
    assert root_gate(straddle, 0.01).status == "unresolved"


def test_exact_observation_radius_enclosure_propagates_to_p_interval():
    cells = [cell(0, np.inf, "selected")]
    bounds = conditional_p_bounds(
        cells, 2, 2, observation_radius_bounds=(Fraction(199, 100), Fraction(201, 100))
    )
    with localcontext() as context:
        context.prec = 120
        exact_low = (-(Decimal("2.01") ** 2) / 2).exp()
        exact_high = (-(Decimal("1.99") ** 2) / 2).exp()
    assert bounds.lower <= exact_low <= exact_high <= bounds.upper
    assert bounds.observation_radius_lower <= Decimal("1.99")
    assert bounds.observation_radius_upper >= Decimal("2.01")


def test_selection_flag_without_exact_observation_bound_cannot_authorize_reject():
    analysis = SimpleNamespace(cells=[cell(0, np.inf, "selected")], selection_certified=True)
    bounds, gate = bounded_root_gate(
        analysis, SimpleNamespace(observed_radius=5, dimension=2), 0.01
    )
    assert gate.status == "unresolved" and not bounds.selection_validated
