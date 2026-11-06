"""Validated spectral/linkage cells for the declared real-arithmetic diffusion map.

Floating-point eigendecomposition proposes a reference only. Exact rational
residual, orthogonality, spectral-gap and distance bounds certify any label.
The public scanner obtains operator enclosures from the actual kernel encloser;
replay samples never certify a cell. Unresolved cells, including infinity, remain.
"""

from dataclasses import dataclass
from decimal import ROUND_CEILING, ROUND_FLOOR, Decimal, localcontext
from fractions import Fraction
from itertools import combinations
from math import isqrt

import numpy as np


@dataclass(frozen=True)
class SelectionCell:
    lo: float
    hi: float
    status: str
    proof: str
    details: dict


@dataclass(frozen=True)
class CellAnalysis:
    cells: tuple[SelectionCell, ...]
    selection_certified: bool
    source_hashes: dict
    observed_radius_bounds: tuple | None = None
    arithmetic_contract: str = "real-arithmetic constructor with validated outward enclosures; not bitwise floating-point replay"


def _fraction(value):
    return value if isinstance(value, Fraction) else Fraction(value)


def _report_bound(value, *, upper=True):
    """Compact outward decimal evidence; internal calculations remain rational."""
    value = _fraction(value)
    with localcontext() as context:
        context.prec = 40
        context.rounding = ROUND_CEILING if upper else ROUND_FLOOR
        return str(Decimal(value.numerator) / Decimal(value.denominator))


def _sqrt_up(value, bits=100):
    """An exact rational upper bound, using integer arithmetic only."""
    value = _fraction(value)
    if value < 0:
        raise ValueError("negative square root")
    scale = 1 << bits
    numerator = value.numerator * value.denominator * scale * scale
    root = isqrt(numerator)
    if root * root < numerator:
        root += 1
    return Fraction(root, value.denominator * scale)


def _sqrt_down(value, bits=100):
    value = _fraction(value)
    if value < 0:
        raise ValueError("negative square root")
    scale = 1 << bits
    return Fraction(
        isqrt(value.numerator * value.denominator * scale * scale), value.denominator * scale
    )


def _frobenius_upper(matrix):
    return _sqrt_up(sum(value * value for row in matrix for value in row))


@dataclass(frozen=True)
class _SpectralReference:
    eigenvalues: tuple
    reconstructed: tuple
    factor: tuple
    orthogonality_error: Fraction
    maximum_absolute_eigenvalue: Fraction
    start: int
    end: int
    time: int


def _spectral_reference(midpoint, target_dimension, time=3):
    """Produce candidate eigenpairs and validate their exact dyadic residual data."""
    matrix = np.asarray(midpoint, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or not np.isfinite(matrix).all():
        raise ValueError("finite square midpoint required")
    n = len(matrix)
    if (
        not isinstance(target_dimension, int)
        or isinstance(target_dimension, bool)
        or not isinstance(time, int)
        or isinstance(time, bool)
        or not 1 <= target_dimension < n - 1
        or time < 1
    ):
        raise ValueError("restricted certificate requires 1<=target<n-1 and time>=1")
    values, vectors = np.linalg.eigh((matrix + matrix.T) / 2)
    values, vectors = values[::-1], vectors[:, ::-1]
    eigenvalues = tuple(Fraction(float(v)) for v in values)
    q = [[Fraction(float(v)) for v in row] for row in vectors]
    orthogonality = [
        [sum(q[k][i] * q[k][j] for k in range(n)) - int(i == j) for j in range(n)] for i in range(n)
    ]
    epsilon = _frobenius_upper(orthogonality)
    if epsilon >= Fraction(1, 10):
        raise ValueError("reference eigenvectors cannot certify a nearby orthogonal basis")
    reconstructed = tuple(
        tuple(sum(q[i][k] * eigenvalues[k] * q[j][k] for k in range(n)) for j in range(n))
        for i in range(n)
    )
    start, end = 1, 1 + int(target_dimension)
    factor = tuple(
        tuple(
            sum(q[i][k] * eigenvalues[k] ** time * q[j][k] for k in range(start, end))
            for j in range(n)
        )
        for i in range(n)
    )
    return _SpectralReference(
        eigenvalues, reconstructed, factor, epsilon, max(map(abs, eigenvalues)), start, end, time
    )


def _distance_enclosures(reference, lower, upper):
    """Bound the selected diffusion factor by Weyl and spectral projectors.

    This restricted route requires an isolated leading eigenvalue, no numerical
    tie at the retained cutoff and strictly positive retained eigenvalues.
    The orthogonal polar factor U of Q obeys ||Q-U||<=epsilon when
    epsilon=||Q.T Q-I||F<1. Thus A*=U diag(lambda) U.T is an exact symmetric
    reference with the stated eigenvalues. All scalar bounds below are rational.
    """
    n = len(reference.eigenvalues)
    if len(lower) != n or len(upper) != n:
        raise ValueError("operator/reference shape mismatch")
    lo = [[_fraction(v) for v in row] for row in lower]
    hi = [[_fraction(v) for v in row] for row in upper]
    if any(len(row) != n for row in lo + hi):
        raise ValueError("square interval matrices required")
    if any(lo[i][j] > hi[i][j] for i in range(n) for j in range(n)):
        raise ValueError("reversed operator interval")
    # Intersect transpose enclosures: the declared real operator is symmetric.
    lo = [[max(lo[i][j], lo[j][i]) for j in range(n)] for i in range(n)]
    hi = [[min(hi[i][j], hi[j][i]) for j in range(n)] for i in range(n)]
    if any(lo[i][j] > hi[i][j] for i in range(n) for j in range(n)):
        raise ValueError("operator intervals exclude symmetry")
    residual = [
        [
            max(
                abs(lo[i][j] - reference.reconstructed[i][j]),
                abs(hi[i][j] - reference.reconstructed[i][j]),
            )
            for j in range(n)
        ]
        for i in range(n)
    ]
    epsilon = reference.orthogonality_error
    largest = reference.maximum_absolute_eigenvalue
    polar_factor_error = epsilon * (2 + epsilon)
    delta = _frobenius_upper(residual) + largest * polar_factor_error
    lam = reference.eigenvalues
    start, end = reference.start, reference.end
    norm_bound = largest + delta
    # np.isclose policy: conservative upper bound on every tie threshold.
    tie_threshold = (Fraction(1e-10) + Fraction(1e-12)) * norm_bound
    gaps = [lam[0] - lam[1], lam[end - 1] - lam[end]]
    if any(gap <= 2 * delta + tie_threshold for gap in gaps):
        raise ValueError("leading or retained spectral group is not uniformly separated")
    if lam[end - 1] <= delta:
        raise ValueError("retained positive spectrum is not uniformly separated from zero")
    # Davis--Kahan operator-norm form: Yu/Wang/Samworth (2015), Theorem 1.
    # https://www.statslab.cam.ac.uk/~rjs57/Biometrika-2015-Yu-315-23.pdf
    projector_error = sum(2 * delta / (gap - 2 * delta) for gap in gaps)
    projector_error = min(Fraction(1), projector_error)
    power = reference.time
    global_factor_error = (
        power * norm_bound ** (power - 1) * delta
        + largest**power * projector_error
        + largest**power * polar_factor_error
    )
    selected_maximum = max(map(abs, lam[start:end]))
    selected_perturbed_maximum = selected_maximum + delta
    # Project the power telescope before bounding; the removed leading mode contributes no power.
    restricted_factor_error = (
        (selected_perturbed_maximum**power + selected_maximum**power) * projector_error
        + delta
        * sum(
            selected_perturbed_maximum ** (power - 1 - k) * selected_maximum**k
            for k in range(power)
        )
        + selected_maximum**power * polar_factor_error
    )
    factor_error = min(global_factor_error, restricted_factor_error)
    distance_lower = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    distance_upper = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    pair_error = _sqrt_up(Fraction(2)) * factor_error
    for i, j in combinations(range(n), 2):
        squared = sum((reference.factor[k][i] - reference.factor[k][j]) ** 2 for k in range(n))
        # Distances are ||H(e_i-e_j)||, so perturbation <=sqrt(2)||Delta H||.
        low = max(Fraction(0), _sqrt_down(squared) - pair_error)
        high = _sqrt_up(squared) + pair_error
        distance_lower[i][j] = distance_lower[j][i] = low
        distance_upper[i][j] = distance_upper[j][i] = high
    return (
        distance_lower,
        distance_upper,
        dict(
            operator_error=_report_bound(delta),
            factor_error=_report_bound(factor_error),
            global_factor_error=_report_bound(global_factor_error),
            restricted_factor_error=_report_bound(restricted_factor_error),
            factor_bound_method=(
                "selected_spectrum_telescope"
                if restricted_factor_error <= global_factor_error
                else "global_spectrum_telescope"
            ),
            pair_distance_error=_report_bound(pair_error),
            distance_bound_method="diffusion_factor_linear_error",
            spectral_gap_lower=_report_bound(min(gaps) - 2 * delta, upper=False),
            retained_indices=(start, end),
        ),
    )


def _classify_linkage(lower, upper, signature, tie_tolerance):
    """Uniform full-history classification from interval average distances."""
    n = len(lower)
    groups = {i: (i,) for i in range(n)}
    active = set(groups)
    smallest_margin = None
    tolerance = Fraction(float(tie_tolerance)) * max(
        Fraction(1), max(value for row in upper for value in row)
    )
    ambiguous = False
    for step, chosen in enumerate(signature):
        if len(chosen) != 2 or not set(chosen) <= active:
            raise ValueError("invalid observed merge sequence")
        a, b = chosen
        chosen_low = sum(lower[i][j] for i in groups[a] for j in groups[b]) / (
            len(groups[a]) * len(groups[b])
        )
        chosen_high = sum(upper[i][j] for i in groups[a] for j in groups[b]) / (
            len(groups[a]) * len(groups[b])
        )
        for c, d in combinations(sorted(active), 2):
            if (c, d) == chosen:
                continue
            count = len(groups[c]) * len(groups[d])
            competitor_low = sum(lower[i][j] for i in groups[c] for j in groups[d]) / count
            competitor_high = sum(upper[i][j] for i in groups[c] for j in groups[d]) / count
            margin_lower = competitor_low - chosen_high
            smallest_margin = (
                margin_lower if smallest_margin is None else min(smallest_margin, margin_lower)
            )
            if competitor_high < chosen_low:
                return "nonselected", dict(
                    step=step,
                    chosen=chosen,
                    competitor=(c, d),
                    margin_upper=_report_bound(competitor_high - chosen_low),
                )
            ambiguous |= margin_lower <= tolerance
        groups[n + step] = groups[a] + groups[b]
        active.difference_update(chosen)
        active.add(n + step)
    if len(active) != 1:
        raise ValueError("merge sequence does not finish at one root")
    return ("unresolved" if ambiguous else "selected"), dict(
        minimum_margin_lower=_report_bound(smallest_margin, upper=False),
        tie_tolerance_upper=_report_bound(tolerance),
    )


def _exact_observed_radius_bounds(path):
    """Outward radius for the exact Gaussian decomposition of dyadic input X."""
    n = len(path.values)
    groups = {i: (i,) for i in range(n)}
    for step, (a, b) in enumerate(path.observed_signature):
        groups[n + step] = groups[a] + groups[b]
    a, b = path.observed_signature[-1]
    eta = [
        Fraction(1, len(groups[a])) if i in groups[a] else -Fraction(1, len(groups[b]))
        for i in range(n)
    ]
    squared_norm = sum(value * value for value in eta)
    contrast = [
        sum(eta[i] * Fraction(float(path.values[i, j])) for i in range(n))
        for j in range(path.dimension)
    ]
    squared_radius = sum(v * v for v in contrast) / squared_norm / Fraction(float(path.sigma)) ** 2
    return _sqrt_down(squared_radius), _sqrt_up(squared_radius)


def analyze_diffusion_root_cells(
    path,
    *,
    radius_bound=None,
    observed_half_width=None,
    subdivision_depth=0,
    anchor_ratios=(1.0,),
    expansion_steps=0,
):
    """Cover [0,infinity] with audited interval cells; retain all unresolved mass.

    The default deliberately includes a small cell around the observed radius
    because broad cells often cross adaptive branches. A small cell is labeled
    only by upstream interval and downstream exact residual proofs, never by its
    sampled signature. Tail stability is not assumed. This is partial coverage,
    not a complete enumeration of every boundary.
    """
    import hashlib
    from pathlib import Path

    from . import diffusion_root_kernel, diffusion_root_mass

    if not isinstance(subdivision_depth, int) or subdivision_depth < 0:
        raise ValueError("subdivision_depth must be a nonnegative integer")
    if not isinstance(expansion_steps, int) or expansion_steps < 0:
        raise ValueError("expansion_steps must be a nonnegative integer")
    observed = float(path.observed_radius)
    bound = max(1.0, 2 * observed) if radius_bound is None else float(radius_bound)
    width = (
        max(1e-12, observed * 1e-13) if observed_half_width is None else float(observed_half_width)
    )
    if not np.isfinite(bound) or bound <= observed or not np.isfinite(width) or width <= 0:
        raise ValueError("finite radius_bound>observed_radius and positive finite width required")
    for name, expected in path.source_hashes.items():
        if hashlib.sha256(Path(name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"replay source changed since path creation: {name}")
    ratios = tuple(float(ratio) for ratio in anchor_ratios)
    if any(not np.isfinite(ratio) or ratio <= 0 for ratio in ratios):
        raise ValueError("anchor_ratios must be positive finite numbers")
    anchors = sorted(
        set(observed * ratio for ratio in (1.0,) + ratios if 0 < observed * ratio < bound)
    )
    cells, cache, references = [], {}, {}

    def certificate(lo, hi, reference=None):
        if (lo, hi) in cache:
            return cache[lo, hi]
        details = {}
        try:
            enclosure = diffusion_root_kernel.enclose_diffusion_operator(path, lo, hi)
            midpoint = [
                [float((_fraction(a) + _fraction(b)) / 2) for a, b in zip(lrow, urow)]
                for lrow, urow in zip(enclosure.lower, enclosure.upper)
            ]
            reused_reference = reference is not None
            if reference is None:
                reference = _spectral_reference(
                    midpoint, enclosure.retained_count_target, enclosure.time
                )
            references[lo, hi] = reference
            lower, upper, spectral_details = _distance_enclosures(
                reference, enclosure.lower, enclosure.upper
            )
            status, linkage_details = _classify_linkage(
                lower, upper, path.observed_signature, path.tie_tolerance
            )
            details = dict(
                kernel=enclosure.proof_details,
                spectral=spectral_details,
                linkage=linkage_details,
                spectral_reference_reused=reused_reference,
            )
            if status != "unresolved":
                result = SelectionCell(
                    lo,
                    hi,
                    status,
                    "outward_native_kernel_and_exact_spectral_linkage_bounds",
                    details,
                )
                cache[lo, hi] = result
                return result
            details["reason"] = "interval linkage comparisons do not separate selection and ties"
        except (ValueError, ArithmeticError) as error:
            details["reason"] = str(error)
        result = SelectionCell(lo, hi, "unresolved", "no_uniform_certificate", details)
        cache[lo, hi] = result
        return result

    def classify(lo, hi, depth):
        result = certificate(lo, hi)
        midpoint = lo + (hi - lo) / 2
        if result.status == "unresolved" and depth > 0 and lo < midpoint < hi:
            classify(lo, midpoint, depth - 1)
            classify(midpoint, hi, depth - 1)
        else:
            cells.append(result)

    splits = {0.0, bound}
    for index, anchor in enumerate(anchors):
        left_limit = 0.0 if index == 0 else (anchors[index - 1] + anchor) / 2
        right_limit = bound if index == len(anchors) - 1 else (anchor + anchors[index + 1]) / 2
        current_width = width
        lo, hi = max(left_limit, anchor - current_width), min(right_limit, anchor + current_width)
        current = certificate(lo, hi)
        anchor_reference = references.get((lo, hi))
        for step in range(expansion_steps):
            if current.status == "unresolved":
                break
            current_width *= 2
            trial_lo = max(left_limit, anchor - current_width)
            trial_hi = min(right_limit, anchor + current_width)
            if (trial_lo, trial_hi) == (lo, hi):
                break
            trial = certificate(trial_lo, trial_hi, anchor_reference)
            if trial.status != current.status:
                break
            lo, hi, current = trial_lo, trial_hi, trial
            current.details["validated_expansion_steps"] = step + 1
        splits.update([lo, hi])
    splits = sorted(splits)

    for lo, hi in zip(splits[:-1], splits[1:]):
        if lo < hi:
            classify(lo, hi, subdivision_depth)
    cells.append(SelectionCell(bound, np.inf, "unresolved", "uncertified_unbounded_tail", {}))
    hashes = dict(path.source_hashes)
    for filename in (
        Path(__file__),
        Path(diffusion_root_kernel.__file__),
        Path(diffusion_root_mass.__file__),
    ):
        hashes[str(filename.resolve())] = hashlib.sha256(filename.read_bytes()).hexdigest()
    observed_bounds = _exact_observed_radius_bounds(path)
    observed_event_validated = any(
        cell.status == "selected"
        and Fraction(cell.lo) <= observed_bounds[0]
        and Fraction(cell.hi) >= observed_bounds[1]
        for cell in cells
    )
    return CellAnalysis(
        tuple(cells),
        observed_event_validated,
        hashes,
        observed_radius_bounds=observed_bounds,
    )
