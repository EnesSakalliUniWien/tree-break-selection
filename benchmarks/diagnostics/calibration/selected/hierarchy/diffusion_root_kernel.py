"""Outward enclosures for exact-real continuous-preset diffusion kernel formulas.

Not bitwise NumPy certificates. Every neighbor/BGH/dimension branch must be
proved on the entire cell. Ambiguity raises ValueError, never a selected label.
The original dyadic observations define the exact Gaussian conditional path;
floating replay offsets/slopes only propose branches and are never trusted.
"""

import hashlib
from dataclasses import dataclass
from decimal import ROUND_CEILING, ROUND_FLOOR, Decimal, localcontext
from fractions import Fraction
from pathlib import Path

from tree_break_selection.space_separation.diffusion import DENSE_DIFFUSION_EIGENSOLVER_MAX_SAMPLES

from .diffusion_root_mass import _add, _rounded, _sqrt, _sub

PRECISION = 80
ZERO = Decimal(0)
ONE = Decimal(1)


def _i(value):
    q = Fraction(value)
    return (
        _rounded(lambda: Decimal(q.numerator) / Decimal(q.denominator), ROUND_FLOOR),
        _rounded(lambda: Decimal(q.numerator) / Decimal(q.denominator), ROUND_CEILING),
    )


def _mul(a, b):
    return (
        min(_rounded(lambda: x * y, ROUND_FLOOR) for x in a for y in b),
        max(_rounded(lambda: x * y, ROUND_CEILING) for x in a for y in b),
    )


def _div(a, b):
    if b[0] <= 0:
        raise ValueError("unresolved positive interval denominator")
    reciprocal = (
        _rounded(lambda: ONE / b[1], ROUND_FLOOR),
        _rounded(lambda: ONE / b[0], ROUND_CEILING),
    )
    return _mul(a, reciprocal)


def _exp(a):
    with localcontext() as context:
        context.prec = PRECISION
        return max(ZERO, a[0].exp().next_minus()), a[1].exp().next_plus()


def _ln(a):
    if a[0] <= 0:
        raise ValueError("unresolved logarithm domain")
    with localcontext() as context:
        context.prec = PRECISION
        return a[0].ln().next_minus(), a[1].ln().next_plus()


def _power(a, b):
    return _exp(_mul(_ln(a), b))


def _sum(items):
    total = (ZERO, ZERO)
    for item in items:
        total = _add(total, item)
    return total


def _exact_path(path):
    x = [[Fraction(float(v)) for v in row] for row in path.values]
    n, p = len(x), len(x[0])
    groups = {i: [i] for i in range(n)}
    for step, (a, b) in enumerate(path.observed_signature):
        groups[n + step] = groups[a] + groups[b]
    a, b = path.observed_signature[-1]
    eta = [
        Fraction(1, len(groups[a])) if i in groups[a] else -Fraction(1, len(groups[b]))
        for i in range(n)
    ]
    h = sum(v * v for v in eta)
    v = [sum(eta[i] * x[i][j] for i in range(n)) for j in range(p)]
    vsq = sum(t * t for t in v)
    if not vsq > 0:
        raise ValueError("unresolved zero exact Gaussian contrast")
    offset = [[x[i][j] - eta[i] * v[j] / h for j in range(p)] for i in range(n)]
    slope_numerator = [[Fraction(path.sigma) * eta[i] * v[j] for j in range(p)] for i in range(n)]
    radius = _sqrt(_i(vsq / h / Fraction(path.sigma) ** 2))
    return offset, slope_numerator, h * vsq, radius


def exact_path_radius_bounds(path):
    """Decimal enclosure of exact sqrt(||X^T eta||²/(sigma²||eta||²))."""
    return _exact_path(path)[3]


def _quadratic_range(constant, linear, quadratic, lo, hi):
    """Exact extrema with an enclosed linear coefficient and nonnegative radius."""
    lo, hi = Fraction(lo), Fraction(hi)
    low_b, high_b = map(Fraction, linear)

    def value(r, b):
        return constant + b * r + quadratic * r * r

    lower = min(value(lo, low_b), value(hi, low_b))
    if quadratic > 0:
        vertex = -low_b / (2 * quadratic)
        if lo < vertex < hi:
            lower = min(lower, value(vertex, low_b))
    upper = max(value(lo, high_b), value(hi, high_b))
    return _i(lower)[0], _i(upper)[1]


def _radial_squared_distances(path, lo, hi):
    """Enclose pair polynomials before interval evaluation to cancel shared motion."""
    offset, numerator, norm_squared, radius = _exact_path(path)
    norm = _sqrt(_i(norm_squared))
    n = len(offset)
    distances = [[(ZERO, ZERO) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(i):
            delta = [a - b for a, b in zip(offset[i], offset[j])]
            velocity = [a - b for a, b in zip(numerator[i], numerator[j])]
            constant = sum(v * v for v in delta)
            quadratic = sum(v * v for v in velocity) / norm_squared
            linear = _div(_i(2 * sum(a * b for a, b in zip(delta, velocity))), norm)
            if not quadratic:
                bounds = _i(constant)
            else:
                bounds = _quadratic_range(constant, linear, quadratic, lo, hi)
            distances[i][j] = distances[j][i] = (max(ZERO, bounds[0]), bounds[1])
    return distances, radius


def _neighbor_support(row, *, self_index, counts):
    """Prove the two required neighbor sets; ordering inside either set is irrelevant."""
    ordered = sorted(
        (j for j in range(len(row)) if j != self_index), key=lambda j: row[j][0] + row[j][1]
    )
    if any(row[j][0] <= 0 for j in ordered):
        raise ValueError("unresolved positive distances or duplicates")
    for count in set(counts):
        if count < len(ordered) and (
            max(row[j][1] for j in ordered[:count]) >= min(row[j][0] for j in ordered[count:])
        ):
            raise ValueError("unresolved nearest-neighbor support")
    return ordered


@dataclass(frozen=True)
class OperatorEnclosure:
    lower: tuple
    upper: tuple
    retained_count_target: int
    time: int
    proof_details: dict
    observed_radius_bounds: tuple


def enclose_diffusion_operator(path, lo, hi):
    """Enclose actual preset's exact-real normalized operator on a finite cell.

    Supported: 9<=n<=2000 (native dense spectral route), separated branches.
    Work grows with all pairwise coordinate distances and 81 sparse kernel sums;
    this is a computationally expensive diagnostic, not a production path.
    Global KDE density normalization cancels exactly through rho and median
    epsilon; all row-dependent factors and all refitting remain enclosed.
    """
    n, p = path.values.shape
    if not 9 <= n <= DENSE_DIFFUSION_EIGENSOLVER_MAX_SAMPLES:
        raise ValueError(
            f"unresolved kernel size: supported 9<=n<={DENSE_DIFFUSION_EIGENSOLVER_MAX_SAMPLES} native dense route"
        )
    if not 0 <= lo <= hi or not float(hi) < float("inf"):
        raise ValueError("unresolved invalid radial cell")
    k = min(15, n - 1)
    distances, radius = _radial_squared_distances(path, lo, hi)
    neighbors = [
        _neighbor_support(row, self_index=i, counts=(7, k - 1)) for i, row in enumerate(distances)
    ]
    pilot = [_sqrt(_div(_sum(distances[i][j] for j in neighbors[i][:7]), _i(7))) for i in range(n)]
    # Pilot graph excludes self and has k-1 neighbors. Query graph includes self.
    scaled = [
        _div(distances[i][j], _mul(pilot[i], pilot[j]))
        for i in range(n)
        for j in neighbors[i][: k - 1]
    ]
    epsilons = [Fraction(2) ** j for j in range(-40, 41)]
    logs = []
    for epsilon in epsilons:
        denominator = _i(4 * epsilon)
        total = _add(_i(n), _sum(_exp(_mul(_i(-1), _div(d, denominator))) for d in scaled))
        logs.append(_ln(total))
    ln2 = _ln(_i(2))
    slopes = [_div(_sub(b, a), ln2) for a, b in zip(logs[:-1], logs[1:])]
    winner = max(range(80), key=lambda i: slopes[i][0])
    if any(slopes[winner][0] <= s[1] for i, s in enumerate(slopes) if i != winner):
        raise ValueError("unresolved BGH maximizing epsilon")
    twice = _mul(_i(2), slopes[winner])
    dimension = int(round(float((twice[0] + twice[1]) / 2)))
    if twice[0] <= Decimal(dimension) - Decimal(".5") or twice[1] >= Decimal(dimension) + Decimal(
        ".5"
    ):
        raise ValueError("unresolved rounded intrinsic dimension")
    if dimension < 0:
        raise ValueError("unresolved negative intrinsic dimension")
    epsilon = epsilons[winner]
    # q_i omits shared n*(4*pi*epsilon)^(d/2); this cancels from final W.
    density = []
    for i in range(n):
        total = _add(
            _i(1),
            _sum(
                _exp(
                    _mul(
                        _i(-1),
                        _div(_div(distances[i][j], _mul(pilot[i], pilot[j])), _i(4 * epsilon)),
                    )
                )
                for j in neighbors[i][: k - 1]
            ),
        )
        density.append(_div(total, _power(pilot[i], _i(dimension))))
    rho = [_power(q, _i(Fraction(-1, dimension + 2))) for q in density]
    scaled_final = {
        (i, j): _div(distances[i][j], _mul(rho[i], rho[j]))
        for i in range(n)
        for j in neighbors[i][: k - 1]
    }
    # Order statistics are coordinatewise monotone, so branch ordering needn't freeze.
    lowers = sorted(d[0] for d in scaled_final.values())
    uppers = sorted(d[1] for d in scaled_final.values())
    m = len(lowers)
    if m % 2:
        outer = (lowers[m // 2], uppers[m // 2])
    else:
        outer = _div(
            _add((lowers[m // 2 - 1], uppers[m // 2 - 1]), (lowers[m // 2], uppers[m // 2])), _i(2)
        )
    directed = {
        (i, j): _exp(_mul(_i(-1), _div(d, _mul(_i(4), outer))))
        for (i, j), d in scaled_final.items()
    }
    weights = [[(ZERO, ZERO) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(i):
            a, b = directed.get((i, j), (ZERO, ZERO)), directed.get((j, i), (ZERO, ZERO))
            weights[i][j] = weights[j][i] = (max(a[0], b[0]), max(a[1], b[1]))
    degrees = [_sum(row) for row in weights]
    operator = [
        [_div(weights[i][j], _sqrt(_mul(degrees[i], degrees[j]))) for j in range(n)]
        for i in range(n)
    ]
    details = dict(
        kernel_contract="exact-real formulas; not bitwise NumPy enclosure",
        exact_gaussian_path=True,
        self_distance_policy="exact zero; native identity diagonal reset before epsilon",
        distance_bound_method="exact radial quadratics with enclosed linear coefficients",
        neighbor_bound_method="pilot and kernel support sets; no internal ordering restriction",
        exact_real_exponent="rho=q**(-1/(d+2)) with rational exponent; not NumExpr float beta",
        bgh_winner_gap_lower=_sub(
            slopes[winner], (max(s[1] for i, s in enumerate(slopes) if i != winner),) * 2
        )[0],
        rounded_dimension_enclosure=twice,
        source_sha256={
            str(Path(__file__).resolve()): hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        },
        neighbors=tuple(tuple(row[: k - 1]) for row in neighbors),
        neighbor_k=k,
        pilot_bandwidths=tuple(pilot),
        reduced_bandwidths=tuple(rho),
        kde_bgh_index=winner,
        kde_epsilon=epsilon,
        kde_dimension=dimension,
        outer_epsilon_reduced=outer,
        precision=PRECISION,
        global_kde_normalization="cancels from rho-scaled distances / their median",
        interval_source="directed Decimal arithmetic; correctly-rounded exp/ln/sqrt widened one ulp",
    )
    return OperatorEnclosure(
        tuple(tuple(Fraction(v[0]) for v in row) for row in operator),
        tuple(tuple(Fraction(v[1]) for v in row) for row in operator),
        min(50, n - 2, p),
        3,
        details,
        tuple(Fraction(v) for v in radius),
    )
