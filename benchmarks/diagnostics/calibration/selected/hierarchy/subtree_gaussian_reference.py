"""Gaussian common-mean references for a specified subtree."""

import math
from numbers import Integral, Real

import numpy as np
from scipy.linalg import helmert
from scipy.stats import beta


def _rows(indices, n, *, name, minimum):
    try:
        rows = tuple(indices)
    except TypeError as error:
        raise ValueError(f"{name} must be a sequence of row indices") from error
    if (
        len(rows) < minimum
        or any(isinstance(i, bool) or not isinstance(i, Integral) or not 0 <= i < n for i in rows)
        or len(set(rows)) != len(rows)
    ):
        raise ValueError(f"{name} must contain at least {minimum} distinct valid row indices")
    return np.asarray(rows, dtype=np.intp)


def _inputs(values, subtree_indices, alpha):
    if isinstance(alpha, bool) or not isinstance(alpha, Real) or not 0 < alpha < 1:
        raise ValueError("alpha must lie strictly between zero and one")
    x = np.asarray(values, dtype=float)
    if x.ndim != 2 or x.shape[1] < 1 or not np.isfinite(x).all():
        raise ValueError("values must be a finite matrix with at least one feature")
    rows = _rows(subtree_indices, len(x), name="subtree_indices", minimum=2)
    # Membership alone fixes the proposal basis, including for data-ordered input indices.
    return x, np.sort(rows)


def _rotated(state, rng):
    proposed = state.copy()
    if len(state) == 1:
        sign = 1 if rng.integers(2) else -1
        proposed[0] *= sign
        return proposed, (0, 0, 0.0, sign)
    i, j = sorted(map(int, rng.choice(len(state), size=2, replace=False)))
    angle = float(rng.uniform(-np.pi, np.pi))
    sign = 1 if rng.integers(2) else -1
    c, s = math.cos(angle), math.sin(angle)
    proposed[i] = c * state[i] - s * state[j]
    proposed[j] = sign * (s * state[i] + c * state[j])
    return proposed, (i, j, angle, sign)


def subtree_gaussian_chain_test(
    values, *, subtree_indices, score, selection_event=None, steps=999, seed=0, alpha=0.05
):
    """Rank a deterministic full-matrix score under the common-mean Gaussian null.

    Condition on the subtree mean, centered feature Gram matrix, outside rows
    and, if supplied, the caller's selection event. Subtree rows must be iid
    Gaussian with one common mean and covariance and independent of outside
    rows. Unknown covariance is allowed. The event must describe how the tested
    membership and any other conditioned selection were obtained. Omitting it
    requires prespecified membership. Callbacks receive read-only matrices.
    Any score-defining choices captured by callbacks must be prespecified or
    part of the event. A failed evaluation invalidates the rank. No descendant or traversal event
    is inferred here, and native floating-point event membership is uncertified.

    Chain moves store accepted transforms as (i, j, angle, sign), relative to
    the preceding state when moving outward from the observed anchor. An empty
    entry keeps the preceding state. Equal indices denote the one-axis case.
    """
    x, rows = _inputs(values, subtree_indices, alpha)
    for name, value, minimum in (("steps", steps, 1), ("seed", seed, 0)):
        if isinstance(value, bool) or not isinstance(value, Integral) or value < minimum:
            raise ValueError(f"{name} must be an integer >= {minimum}")
    if not callable(score) or (selection_event is not None and not callable(selection_event)):
        raise ValueError("score and any selection_event must be callable")

    def in_event(matrix):
        if selection_event is None:
            return True
        accepted = selection_event(matrix)
        if not isinstance(accepted, (bool, np.bool_)):
            raise ValueError("selection_event must return a boolean")
        return bool(accepted)

    def score_value(matrix):
        value = score(matrix)
        if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
            raise ValueError("score must return a finite real number")
        return float(value)

    rng = np.random.default_rng(seed)
    anchor = int(rng.integers(steps + 1))
    scores = [None] * (steps + 1)
    moves = [None] * (steps + 1)
    result = dict(
        status="unresolved",
        decision="unresolved",
        p_value=None,
        error=None,
        sampling_status="unresolved",
        observed_score=None,
        anchor_index=anchor,
        accepted_moves=0,
        completed_steps=0,
        event_evaluations=0,
        chain_scores=scores,
        chain_moves=moves,
        alpha=float(alpha),
        seed=int(seed),
        steps_requested=int(steps),
        resolution_sufficient=1 / (steps + 1) <= alpha,
        subtree_indices=rows.tolist(),
        numerical_event_certified=False,
        selection_event="fixed_membership" if selection_event is None else "caller_supplied",
        observation_model="independent_gaussian_common_mean_covariance",
    )
    observed = x.copy()
    observed.setflags(write=False)
    phase = "observation"
    try:
        if not in_event(observed):
            raise ValueError("observed matrix does not satisfy selection_event")
        scores[anchor] = result["observed_score"] = score_value(observed)
        with np.errstate(over="raise", invalid="raise"):
            basis = helmert(len(rows), full=False)
            mean = x[rows].mean(axis=0)
            original = basis @ (x[rows] - mean)
            phase = "sampling"
            for positions in (range(anchor - 1, -1, -1), range(anchor + 1, steps + 1)):
                current, current_score = original, scores[anchor]
                for index in positions:
                    proposed, move = _rotated(current, rng)
                    if not np.array_equal(proposed, current):
                        candidate = x.copy()
                        candidate[rows] = mean + basis.T @ proposed
                        candidate.setflags(write=False)
                        if selection_event is not None:
                            result["event_evaluations"] += 1
                        if in_event(candidate):
                            current_score = score_value(candidate)
                            current = proposed
                            moves[index] = move
                            result["accepted_moves"] += 1
                    scores[index] = current_score
                    result["completed_steps"] += 1
    except Exception as error:
        result.update(error_phase=phase, error=f"{type(error).__name__}: {error}")
        return result

    p_value = sum(value >= scores[anchor] for value in scores) / (steps + 1)
    moved = result["accepted_moves"] > 0
    result.update(
        status="complete", p_value=p_value, sampling_status="moved" if moved else "stalled"
    )
    if moved and result["resolution_sufficient"]:
        result["decision"] = "reject" if p_value <= alpha else "do_not_reject"
    return result


def fixed_subtree_contrast_test(
    values, *, subtree_indices, left_indices, rank=None, all_labels=False, alpha=0.05
):
    """Test a fixed binary row contrast on a prespecified common-mean subtree.

    The statistic is squared projection of the unit row contrast onto the
    leading rank-k left singular subspace of the centered subtree. The right
    directions and numerical rank are functions of its feature Gram matrix.
    For fixed membership and labels the conditional law is Beta(k/2,(m-1-k)/2).
    An optional Bonferroni factor counts all unordered nonempty binary splits.
    It permits label search on the specified subtree. Selected membership and
    contrast-selected rank are outside the reference.
    """
    x, rows = _inputs(values, subtree_indices, alpha)
    left = _rows(left_indices, len(x), name="left_indices", minimum=1)
    left_mask = np.isin(rows, left)
    if left_mask.sum() != len(left) or len(left) == len(rows):
        raise ValueError("left_indices must form a nonempty proper subset of subtree_indices")
    if rank is not None and (isinstance(rank, bool) or not isinstance(rank, Integral) or rank < 0):
        raise ValueError("rank must be a nonnegative integer or None")
    if not isinstance(all_labels, (bool, np.bool_)):
        raise ValueError("all_labels must be a boolean")

    m = len(rows)
    count = (1 << (m - 1)) - 1 if all_labels else 1
    result = dict(
        status="unresolved",
        decision="unresolved",
        error=None,
        p_value=None,
        unadjusted_p_value=None,
        statistic=None,
        rank=None,
        subtree_indices=rows.tolist(),
        left_indices=left.tolist(),
        number_of_labels=count,
        alpha=float(alpha),
        selection_event="fixed_membership",
        all_labels=bool(all_labels),
        observation_model="independent_gaussian_common_mean_covariance",
    )
    try:
        with np.errstate(over="raise", invalid="raise"):
            basis = helmert(m, full=False)
            centered = x[rows] - x[rows].mean(axis=0)
            u, singular, _ = np.linalg.svd(basis @ centered, full_matrices=False)
        tolerance = singular[0] * max(centered.shape) * np.finfo(float).eps
        full_rank = int(np.count_nonzero(singular > tolerance))
        k = full_rank if rank is None else int(rank)
        result["rank"] = k
        if not 0 < k < m - 1 or k > full_rank:
            raise ValueError(
                "rank is zero, exceeds the numerical rank, or spans the centered row space"
            )
        if k < full_rank and singular[k - 1] - singular[k] <= tolerance:
            raise ValueError("rank cuts a tied Gram eigenspace")
        contrast = np.where(left_mask, 1 / len(left), -1 / (m - len(left)))
        contrast /= np.linalg.norm(contrast)
        q = float(np.linalg.norm(u[:, :k].T @ (basis @ contrast)) ** 2)
        if not 0 <= q < 1:
            raise FloatingPointError("projection statistic is outside numerical Beta support")
        log_tail = float(beta.logsf(q, k / 2, (m - 1 - k) / 2))
        if not math.isfinite(log_tail):
            raise FloatingPointError("Beta tail underflowed or is nonfinite")
    except (ValueError, FloatingPointError, np.linalg.LinAlgError) as error:
        result["error"] = f"{type(error).__name__}: {error}"
        return result

    p_value = math.exp(min(0.0, log_tail + math.log(count)))
    return result | dict(
        status="complete",
        decision="reject" if p_value <= alpha else "do_not_reject",
        p_value=p_value,
        unadjusted_p_value=math.exp(log_tail),
        statistic=q,
    )
