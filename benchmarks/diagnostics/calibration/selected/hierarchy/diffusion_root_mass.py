"""Outward-rounded chi masses and a conservative diagnostic root-only gate.

Decimal exp/sqrt are correctly rounded to nearest; one adjacent Decimal on each
side encloses their exact values. Algebra uses directed rounding at precision80.
Gamma CDF uses a positive series with geometric remainder. Odd-dimensional
far tails use alternating Taylor remainder bounds for the erfc integral.
No quadrature error estimate or sampled selection label is a certificate.
"""

from dataclasses import dataclass, replace
from decimal import ROUND_CEILING, ROUND_FLOOR, Decimal, localcontext
from fractions import Fraction
from functools import lru_cache
from types import SimpleNamespace

ZERO, ONE = Decimal(0), Decimal(1)
PRECISION = 80


def _d(value):
    if isinstance(value, Decimal):
        return value
    if isinstance(value, (float,)):
        return Decimal.from_float(value)
    return Decimal(str(value))


def _rounded(operation, rounding):
    with localcontext() as context:
        context.prec = PRECISION
        context.rounding = rounding
        return operation()


def _add(a, b):
    return (
        _rounded(lambda: a[0] + b[0], ROUND_FLOOR),
        _rounded(lambda: a[1] + b[1], ROUND_CEILING),
    )


def _sub(a, b):
    return (
        _rounded(lambda: a[0] - b[1], ROUND_FLOOR),
        _rounded(lambda: a[1] - b[0], ROUND_CEILING),
    )


def _mul(a, b):
    # All callers multiply nonnegative intervals.
    return (
        _rounded(lambda: a[0] * b[0], ROUND_FLOOR),
        _rounded(lambda: a[1] * b[1], ROUND_CEILING),
    )


def _div(a, b):
    if b[0] <= 0:
        raise ValueError("Interval division requires positive denominator lower bound.")
    return (
        _rounded(lambda: a[0] / b[1], ROUND_FLOOR),
        _rounded(lambda: a[1] / b[0], ROUND_CEILING),
    )


def _power(a, exponent):
    result = (ONE, ONE)
    while exponent:
        if exponent % 2:
            result = _mul(result, a)
        exponent //= 2
        if exponent:
            a = _mul(a, a)
    return result


def _sqrt(a):
    with localcontext() as context:
        context.prec = PRECISION
        return max(ZERO, a[0].sqrt().next_minus()), a[1].sqrt().next_plus()


def _exp_negative(a):
    with localcontext() as context:
        context.prec = PRECISION
        return max(ZERO, (-a[1]).exp().next_minus()), (-a[0]).exp().next_plus()


def _clip(a):
    return max(ZERO, a[0]), min(ONE, a[1])


@lru_cache(maxsize=1)
def _pi():
    # Machin identity; exact-rational alternating arctan remainders certify pi.
    def atan(q):
        total = sum(
            (Fraction((-1) ** k, (2 * k + 1) * q ** (2 * k + 1)) for k in range(80)), Fraction(0)
        )
        next_term = Fraction(1, 161 * q**161)
        return total, total + next_term

    fifth, two_thirty_ninth = atan(5), atan(239)
    lower = 16 * fifth[0] - 4 * two_thirty_ninth[1]
    upper = 16 * fifth[1] - 4 * two_thirty_ninth[0]
    return (
        _rounded(lambda: Decimal(lower.numerator) / Decimal(lower.denominator), ROUND_FLOOR),
        _rounded(lambda: Decimal(upper.numerator) / Decimal(upper.denominator), ROUND_CEILING),
    )


@lru_cache(maxsize=None)
def _gamma_a_plus_one(dimension):
    if dimension % 2 == 0:
        result = (ONE, ONE)
        for j in range(1, dimension // 2 + 1):
            result = _mul(result, (Decimal(j), Decimal(j)))
        return result
    result = _sqrt(_pi())
    for j in range(dimension // 2 + 1):
        factor = Decimal(j) + Decimal("0.5")
        result = _mul(result, (factor, factor))
    return result


def _cdf_x(x, dimension):
    if x[1] == 0:
        return ZERO, ZERO
    if x[1] > 512:
        return ZERO, ONE
    a = Decimal(dimension) / 2
    x_power = _power(x, dimension // 2)
    if dimension % 2:
        x_power = _mul(x_power, _sqrt(x))
    factor = _div(_mul(_exp_negative(x), x_power), _gamma_a_plus_one(dimension))
    term = total = (ONE, ONE)
    for n in range(10000):
        denominator = a + n + 1
        ratio = _div(x, (denominator, denominator))
        next_term = _mul(term, ratio)
        if ratio[1] < 1:
            remaining = _div(next_term, _sub((ONE, ONE), ratio))
            if remaining[1] <= total[0] * Decimal("1e-65"):
                return _clip(_mul(factor, _add(total, (ZERO, remaining[1]))))
        total = _add(total, next_term)
        term = next_term
    return ZERO, ONE


def _sf_x(x, dimension):
    if x[1] == 0:
        return ONE, ONE
    exponential = _exp_negative(x)
    if dimension % 2 == 0:
        term = total = (ONE, ONE)
        for j in range(1, dimension // 2):
            term = _div(_mul(term, x), (Decimal(j), Decimal(j)))
            total = _add(total, term)
        return _clip(_mul(exponential, total))
    if x[0] < 32:
        return _clip(_sub((ONE, ONE), _cdf_x(x, dimension)))
    # erfc(sqrt(x)) = exp(-x)/(sqrt(pi*x)) * E[(1+E/x)^(-1/2)], E~Exp(1).
    # Taylor's integral remainder has the sign of the next term for E>=0.
    partial = term = (ONE, ONE)
    lower, upper = ZERO, ONE
    for k in range(1, 41):
        term = _div(_mul(term, (Decimal(k) - Decimal("0.5"),) * 2), x)
        partial = _sub(partial, term) if k % 2 else _add(partial, term)
        if k % 2:
            lower = max(lower, partial[0])
        else:
            upper = min(upper, partial[1])
    half = _mul(_div(exponential, _sqrt(_mul(_pi(), x))), (lower, upper))
    total = half
    x_power = _sqrt(x)
    gamma = _div(_sqrt(_pi()), (Decimal(2), Decimal(2)))
    for j in range(dimension // 2):
        total = _add(total, _div(_mul(exponential, x_power), gamma))
        x_power = _mul(x_power, x)
        factor = Decimal(j) + Decimal("1.5")
        gamma = _mul(gamma, (factor, factor))
    return _clip(total)


@lru_cache(maxsize=8192)
def _endpoint(radius, dimension):
    if radius.is_infinite():
        return (ZERO, ZERO), (ONE, ONE)
    x = _div(_power((radius, radius), 2), (Decimal(2), Decimal(2)))
    return _sf_x(x, dimension), x


def chi_interval_mass(lo, hi, dimension):
    """Enclose P(lo < chi_dimension < hi); float endpoints mean exact float values."""
    lo, hi = _d(lo), _d(hi)
    if not isinstance(dimension, int) or isinstance(dimension, bool) or dimension < 1:
        raise ValueError("dimension must be a positive integer.")
    if not lo.is_finite() or lo < 0 or hi.is_nan() or hi < lo:
        raise ValueError("Invalid nonnegative chi interval endpoints.")
    if lo == hi:
        return ZERO, ZERO
    if lo == 0 and hi.is_infinite():
        return ONE, ONE
    survival_lo, x_lo = _endpoint(lo, dimension)
    if hi.is_infinite():
        return survival_lo
    survival_hi, x_hi = _endpoint(hi, dimension)
    mass = _clip(_sub(survival_lo, survival_hi))
    if mass[0] > 0 and mass[1] - mass[0] < mass[0] * Decimal("1e-40"):
        return mass
    cdf_mass = _clip(_sub(_cdf_x(x_hi, dimension), _cdf_x(x_lo, dimension)))
    return max(mass[0], cdf_mass[0]), min(mass[1], cdf_mass[1])


@dataclass(frozen=True)
class ConditionalPBounds:
    lower: Decimal
    upper: Decimal
    selection_mass_lower: Decimal
    selection_mass_upper: Decimal
    selected_below: tuple
    selected_above: tuple
    unresolved_below: tuple
    unresolved_above: tuple
    selection_validated: bool = False
    observation_radius_lower: Decimal | None = None
    observation_radius_upper: Decimal | None = None


@dataclass(frozen=True)
class RootGateResult:
    status: str
    clusters: int
    passthrough: bool = False


def conditional_p_bounds(cells, observed_radius, dimension, *, observation_radius_bounds=None):
    """Bound a conditional tail for declared interval labels, not certify labels.

    Missing intervals, unsupported/tied cells and the missing infinite tail are
    unresolved mass. A scanner's independently audited selection certificates
    are required before the root gate may use a definite decision.
    For an enclosing observed-radius interval, monotonicity bounds the tail
    below at its upper endpoint and above at its lower endpoint. Component
    masses then describe the lower endpoint, while the p bounds cover both.
    """
    if observation_radius_bounds is not None:

        def endpoint(value, rounding):
            if isinstance(value, Fraction):
                return _rounded(
                    lambda: Decimal(value.numerator) / Decimal(value.denominator), rounding
                )
            return _d(value)

        lo = endpoint(observation_radius_bounds[0], ROUND_FLOOR)
        hi = endpoint(observation_radius_bounds[1], ROUND_CEILING)
        if not lo.is_finite() or not hi.is_finite() or not ZERO <= lo <= hi:
            raise ValueError("Observed radius bounds must be finite, nonnegative and ordered.")
        lower_threshold = conditional_p_bounds(cells, lo, dimension)
        upper_threshold = conditional_p_bounds(cells, hi, dimension)
        return replace(
            lower_threshold,
            lower=upper_threshold.lower,
            observation_radius_lower=lo,
            observation_radius_upper=hi,
        )
    observed = _d(observed_radius)
    if not observed.is_finite() or observed < 0:
        raise ValueError("Observed radius must be finite and nonnegative.")
    selected_below = selected_above = unknown_below = unknown_above = (ZERO, ZERO)
    partition, cursor = [], ZERO
    for cell in sorted(cells, key=lambda c: _d(c.lo)):
        lo, hi = _d(cell.lo), _d(cell.hi)
        if lo < cursor:
            raise ValueError("Selection cells overlap.")
        if lo > cursor:
            partition.append((cursor, lo, "unresolved"))
        if hi <= lo or lo < 0 or lo.is_nan() or hi.is_nan():
            raise ValueError("Invalid selection cell endpoints.")
        if cell.status not in {"selected", "nonselected", "unresolved", "unsupported", "tie"}:
            raise ValueError("Unknown selection cell status.")
        partition.append((lo, hi, cell.status))
        cursor = hi
    if cursor.is_finite():
        partition.append((cursor, Decimal("Infinity"), "unresolved"))
    for lo, hi, status in partition:
        if status == "nonselected":
            continue
        below = (
            chi_interval_mass(lo, min(hi, observed), dimension) if lo < observed else (ZERO, ZERO)
        )
        above = (
            chi_interval_mass(max(lo, observed), hi, dimension) if hi > observed else (ZERO, ZERO)
        )
        if status == "selected":
            selected_below, selected_above = (
                _add(selected_below, below),
                _add(selected_above, above),
            )
        else:
            unknown_below, unknown_above = _add(unknown_below, below), _add(unknown_above, above)
    selected = _add(selected_below, selected_above)
    possible = _add(selected, _add(unknown_below, unknown_above))
    lower, upper = ZERO, ONE
    if selected[0] > 0:
        lower_numerator = (selected_above[0], selected_above[0])
        other_below = _add(selected_below, unknown_below)
        low_denominator = _add(lower_numerator, (other_below[1], other_below[1]))
        lower = _div(
            (selected_above[0], selected_above[0]), (low_denominator[1], low_denominator[1])
        )[0]
        high_numerator = _add(selected_above, unknown_above)
        upper_numerator = (high_numerator[1], high_numerator[1])
        high_denominator = _add(upper_numerator, (selected_below[0], selected_below[0]))
        upper = (
            _div(upper_numerator, (high_denominator[0], high_denominator[0]))[1]
            if high_denominator[0] > 0
            else ONE
        )
    return ConditionalPBounds(
        max(ZERO, lower),
        min(ONE, upper),
        selected[0],
        min(ONE, possible[1]),
        selected_below,
        selected_above,
        unknown_below,
        unknown_above,
    )


def root_gate(bounds, alpha):
    """Root-only gate: unresolved/nonreject stop at one cluster, never pass through."""
    threshold = _d(alpha)
    if not threshold.is_finite() or not ZERO < threshold <= ONE:
        raise ValueError("alpha must be finite in (0,1].")
    if not ZERO <= bounds.lower <= bounds.upper <= ONE:
        raise ValueError("Conditional probability bounds must lie in [0,1].")
    if not bounds.selection_validated or bounds.selection_mass_lower <= 0:
        return RootGateResult("unresolved", 1)
    if bounds.upper <= threshold:
        return RootGateResult("reject", 2)
    if bounds.lower > threshold:
        return RootGateResult("nonreject", 1)
    return RootGateResult("unresolved", 1)


def bounded_root_gate(analysis, path, alpha):
    """Consume scanner output; replay-only cells cannot authorize root rejection."""
    import hashlib
    from pathlib import Path

    from .diffusion_root_cells import CellAnalysis

    observed_bounds = getattr(analysis, "observed_radius_bounds", None)
    validated = (
        isinstance(analysis, CellAnalysis)
        and analysis.selection_certified
        and observed_bounds is not None
    )
    if validated:
        hashes = analysis.source_hashes
        replay_hashes = getattr(path, "source_hashes", {})
        validated = bool(replay_hashes) and all(
            hashes.get(name) == value for name, value in replay_hashes.items()
        )
        validated = validated and all(
            Path(name).is_file() and hashlib.sha256(Path(name).read_bytes()).hexdigest() == value
            for name, value in hashes.items()
        )
        validated = validated and all(
            cell.proof == "outward_native_kernel_and_exact_spectral_linkage_bounds"
            for cell in analysis.cells
            if cell.status in {"selected", "nonselected"}
        )
    cells = (
        analysis.cells
        if validated
        else [
            SimpleNamespace(lo=cell.lo, hi=cell.hi, status="unresolved") for cell in analysis.cells
        ]
    )
    bounds = conditional_p_bounds(
        cells, path.observed_radius, path.dimension, observation_radius_bounds=observed_bounds
    )
    bounds = replace(bounds, selection_validated=validated)
    return bounds, root_gate(bounds, alpha)
