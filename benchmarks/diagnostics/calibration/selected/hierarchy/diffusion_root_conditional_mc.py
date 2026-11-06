"""Known-covariance root selective inference through native diffusion replay.

Condition on the full average-linkage merge history, nuisance offset and
contrast direction. The proposal budget is fixed independently of the observed
radius. Conditional on the accepted count, the observation and accepted radii
are exchangeable under the declared Gaussian contrast null.

This diagnostic uses floating replay and has no interval-certification or
production-traversal guarantee. Run serially because the existing replay
temporarily instruments native functions. The chain tests extend from the root
to a pre-set root-plus-children family, each node conditional on its ancestor path.
"""

import hashlib
from dataclasses import replace
from pathlib import Path

import numpy as np
from joblib import parallel_config
from scipy.stats import chi

from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_replay import (
    _geometry,
    make_diffusion_root_path,
)


def selected_diffusion_root_mc_test(values, *, sigma, proposals=999, seed=0):
    """Return a randomized conditional rank p-value for a known Gaussian scale.

    A zero acceptance count gives p=1. Any unresolved membership gives no
    p-value, since silently excluding it would change the conditional law.
    """
    if isinstance(proposals, bool) or not isinstance(proposals, (int, np.integer)):
        raise ValueError("proposals must be a positive integer")
    if proposals < 1:
        raise ValueError("proposals must be a positive integer")
    rng = np.random.default_rng(seed)
    path = make_diffusion_root_path(values, sigma=sigma)
    draws = np.sqrt(rng.chisquare(path.dimension, size=proposals))
    accepted = []
    result = {
        "status": "complete",
        "selection_event": "full_merge_history",
        "sigma": float(sigma),
        "dimension": path.dimension,
        "observed_radius": path.observed_radius,
        "nominal_p_value": float(chi.sf(path.observed_radius, path.dimension)),
        "proposals_requested": int(proposals),
        "seed": seed,
        "observed_signature": path.observed_signature,
        "source_hashes": path.source_hashes,
        "accepted_radii": accepted,
        "p_value": None,
        "minimum_p_value": None,
        "strict_exclusions_with_ties": 0,
    }
    for index, radius in enumerate(draws):
        point = path.evaluate(radius)
        result["proposals_evaluated"] = index + 1
        if point.status == "selected":
            accepted.append(float(radius))
        elif (
            point.status == "tie"
            and point.signature != path.observed_signature
            and np.any(point.margins < -point.state["merge_tie_tolerance"])
        ):
            # One violated merge inequality excludes the complete observed history.
            result["strict_exclusions_with_ties"] += 1
        elif point.status != "nonselected":
            result.update(status="unresolved", error=point.error, replay_status=point.status)
            break
    result["accepted_count"] = len(accepted)
    if result["status"] == "complete":
        exceedances = sum(radius >= path.observed_radius for radius in accepted)
        result["p_value"] = (1 + exceedances) / (1 + len(accepted))
        result["minimum_p_value"] = 1 / (1 + len(accepted))
    return result


def _root_partition(signature, n):
    groups = {i: (i,) for i in range(n)}
    for step, (a, b) in enumerate(signature):
        groups[n + step] = tuple(sorted(groups[a] + groups[b]))
    a, b = signature[-1]
    return tuple(sorted((groups[a], groups[b])))


def _node_rows(signature, n, node):
    """Row indices under ``node`` in a merge signature."""
    groups = {i: (i,) for i in range(n)}
    for step, (a, b) in enumerate(signature):
        groups[n + step] = tuple(sorted(groups[a] + groups[b]))
    return groups[node]


def _merge_pairs(signature, n):
    """Unordered child row-set pairs of every merge in a signature."""
    groups = {i: (i,) for i in range(n)}
    pairs = set()
    for step, (a, b) in enumerate(signature):
        pairs.add(tuple(sorted((groups[a], groups[b]))))
        groups[n + step] = tuple(sorted(groups[a] + groups[b]))
    return pairs


def _ancestor_path_pairs(signature, n, node):
    """Child pairs of ``node`` and of each ancestor: the splits that select the node."""
    groups = {i: (i,) for i in range(n)}
    parent, pair = {}, {}
    for step, (a, b) in enumerate(signature):
        merged = n + step
        parent[a] = parent[b] = merged
        pair[merged] = tuple(sorted((groups[a], groups[b])))
        groups[merged] = tuple(sorted(groups[a] + groups[b]))
    required = {pair[node]}
    while node in parent:
        node = parent[node]
        required.add(pair[node])
    return required


def _node_path(path, node):
    """Radial path along one node's child contrast; rows outside the node stay fixed."""
    n = len(path.values)
    left, right = path.observed_signature[node - n]
    contrast = np.zeros(n)
    for child, sign in ((left, 1.0), (right, -1.0)):
        rows = list(_node_rows(path.observed_signature, n, child))
        contrast[rows] = sign / len(rows)
    e = contrast / np.linalg.norm(contrast)
    vector = e @ path.values
    norm = float(np.linalg.norm(vector))
    if not norm > 0:
        raise ValueError("node contrast direction undefined")
    u = vector / norm
    return replace(
        path,
        e=e,
        u=u,
        observed_radius=norm / path.sigma,
        offset=path.values - np.outer(e, vector),
        slope=path.sigma * np.outer(e, u),
    )


def _validated_chain_inputs(values, *, steps, seed, sigma, alpha):
    for name, value, lower in [("steps", steps, 1), ("seed", seed, 0)]:
        if (
            isinstance(value, (bool, np.bool_))
            or not isinstance(value, (int, np.integer))
            or value < lower
        ):
            raise ValueError(f"{name} must be an integer >= {lower}")
    for name, value in [("sigma", sigma), ("alpha", alpha)]:
        if (
            isinstance(value, (bool, np.bool_))
            or not isinstance(value, (int, float, np.integer, np.floating))
            or not np.isfinite(value)
            or value <= 0
        ):
            raise ValueError(f"{name} must be a finite positive scalar")
    if alpha >= 1:
        raise ValueError("alpha must be strictly between zero and one")
    x = np.asarray(values, dtype=float)
    if x.ndim != 2 or len(x) < 2 or x.shape[1] < 1 or not np.isfinite(x).all():
        raise ValueError("values must be a finite n-by-p matrix with n>=2 and p>=1")
    return x, int(steps), int(seed)


def _sample_chi_chain(observed, dimension, event, *, steps, seed):
    """Random-position reversible chi rank; failed event evaluations invalidate the run."""
    rng = np.random.default_rng(seed)
    anchor = int(rng.integers(steps + 1))
    states = np.full(steps + 1, np.nan)
    states[anchor] = observed
    result = dict(
        status="complete",
        p_value=None,
        anchor_index=anchor,
        accepted_moves=0,
        geometry_evaluations=0,
        completed_steps=0,
        error=None,
    )
    try:
        for indices, source_delta in [
            (range(anchor - 1, -1, -1), 1),
            (range(anchor + 1, steps + 1), -1),
        ]:
            for index in indices:
                current = states[index + source_delta]
                proposed = current + rng.normal()
                states[index] = current
                if not np.isfinite(proposed):
                    raise FloatingPointError("nonfinite radius proposal")
                if proposed > 0:
                    log_ratio = (dimension - 1) * np.log(proposed / current) - (
                        proposed**2 - current**2
                    ) / 2
                    if not np.isfinite(log_ratio):
                        raise FloatingPointError("nonfinite chi density ratio")
                    if np.log(rng.random()) <= min(0.0, log_ratio):
                        result["geometry_evaluations"] += 1
                        if event(proposed):
                            states[index] = proposed
                            result["accepted_moves"] += int(proposed != current)
                result["completed_steps"] += 1
    except (
        ValueError,
        RuntimeError,
        np.linalg.LinAlgError,
        FloatingPointError,
        OverflowError,
    ) as error:
        result.update(
            status="unresolved",
            error=f"{type(error).__name__}: {error}",
            failed_radius=float(proposed),
        )
    if result["status"] == "complete":
        result["p_value"] = float(np.mean(states >= observed))
    result.update(
        chain_radii=[float(x) if np.isfinite(x) else None for x in states],
        unique_states=int(len(np.unique(states[np.isfinite(states)]))),
        sampling_status=(
            "unresolved"
            if result["status"] != "complete"
            else "moved"
            if result["accepted_moves"]
            else "stalled"
        ),
    )
    return result


def selected_diffusion_root_chain_test(values, *, sigma, steps=999, seed=0, alpha=0.05):
    """Test one selected diffusion root with a fixed random-position chi chain.

    Requires independent Gaussian rows with externally known common covariance
    ``sigma**2 * I`` in the supplied features. Fitted PCA, fitted covariance and
    dependent rows are outside this contract. The full feature dimension is used.

    The event preserves the unordered root children, not every merge. The observed
    path uses the existing strict constructor; proposals use native deterministic
    floating-point tie breaking. No exact-real event certificate is claimed.

    The unit proposal scale, step budget and event must be fixed before inspecting
    the result. The observation occupies an independent uniform chain position;
    all upper-tail ties count. No burn-in, thinning, retries or convergence-based
    stopping is used. Poor movement limits power, not a license to rerun until
    significance. A stalled run retains its conservative p=1 but gives no decision.

    Invalid numerical inputs raise ValueError. Unsupported observed geometry or
    a required proposal failure returns an unresolved result with no p-value.
    Run serially because the shared native replay instruments module functions.
    """
    x, steps, seed = _validated_chain_inputs(
        values, steps=steps, seed=seed, sigma=sigma, alpha=alpha
    )
    result = dict(
        status="unresolved",
        decision="unresolved",
        p_value=None,
        selection_event="root_split",
        observation_model="independent_gaussian_known_spherical_covariance",
        sigma=float(sigma),
        alpha=float(alpha),
        dimension=x.shape[1],
        steps_requested=steps,
        seed=seed,
        proposal_scale=1.0,
        minimum_p_value=1 / (steps + 1),
        resolution_sufficient=1 / (steps + 1) <= alpha,
        numerical_event_certified=False,
        error=None,
    )
    try:
        path = make_diffusion_root_path(values, sigma=float(sigma))
    except (ValueError, RuntimeError, np.linalg.LinAlgError) as error:
        result.update(error_phase="observation", error=f"{type(error).__name__}: {error}")
        return result
    root = _root_partition(path.observed_signature, len(path.values))
    result.update(
        observed_radius=path.observed_radius,
        observed_signature=path.observed_signature,
        root_children=tuple(tuple(path.labels[i] for i in group) for group in root),
        source_hashes={
            **path.source_hashes,
            str(Path(__file__).resolve()): hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        },
    )

    def event(radius):
        signature, _ = _geometry(path.offset + radius * path.slope, path.labels)
        return _root_partition(signature, len(path.values)) == root

    sampled = _sample_chi_chain(path.observed_radius, path.dimension, event, steps=steps, seed=seed)
    result.update(sampled)
    if sampled["status"] != "complete":
        result["error_phase"] = "sampling"
    elif sampled["sampling_status"] == "moved" and result["resolution_sufficient"]:
        result["decision"] = "reject" if sampled["p_value"] <= alpha else "do_not_reject"
    return result


def _node_chain_test(path, node, *, steps, seed, alpha):
    n = len(path.values)
    node_path = path if node == 2 * n - 2 else _node_path(path, node)
    required = _ancestor_path_pairs(path.observed_signature, n, node)

    def event(radius):
        signature, _ = _geometry(node_path.offset + radius * node_path.slope, path.labels)
        return required <= _merge_pairs(signature, n)

    # Sequential joblib leaves merge signatures unchanged but avoids about 70 ms of worker
    # start-up per low-dimensional neighbor search, which dominates small replays.
    with parallel_config(backend="sequential"):
        result = _sample_chi_chain(
            node_path.observed_radius, node_path.dimension, event, steps=steps, seed=seed
        )
    labels = lambda node_id: tuple(  # noqa: E731
        path.labels[i] for i in _node_rows(path.observed_signature, n, node_id)
    )
    result.update(
        node=int(node),
        rows=labels(node),
        children=tuple(labels(child) for child in path.observed_signature[node - n]),
        observed_radius=node_path.observed_radius,
        seed=int(seed),
        decision="unresolved",
    )
    if result["sampling_status"] == "moved" and 1 / (steps + 1) <= alpha:
        result["decision"] = "reject" if result["p_value"] <= alpha else "do_not_reject"
    return result


def selected_diffusion_family_chain_tests(values, *, sigma, steps=999, seed=0, alpha=0.05):
    """Test the root and each internal root child, each conditional on its ancestor path.

    The family is pre-set: the root, then the root's internal children in signature
    order. Each node tests its two children's mean contrast with the random-position
    chi chain on its own contrast radius; only the node's rows move. A proposal is in
    the event when every split from the root down to the node recurs, which is what
    places the node in the family; at the root this is the root-split event of
    ``selected_diffusion_root_chain_test``. Decisions use the fixed Bonferroni level
    ``alpha / 3`` whatever the number of internal children; position ``j`` uses chain
    seed ``seed + j``. Same model and limits as the root test: independent Gaussian
    rows with known covariance ``sigma**2 * I``; no production-traversal guarantee.
    """
    x, steps, seed = _validated_chain_inputs(
        values, steps=steps, seed=seed, sigma=sigma, alpha=alpha
    )
    result = dict(
        status="unresolved",
        selection_event="ancestor_path",
        observation_model="independent_gaussian_known_spherical_covariance",
        sigma=float(sigma),
        alpha=float(alpha),
        family_alpha=float(alpha) / 3,
        steps_requested=steps,
        seed=seed,
        nodes=[],
        numerical_event_certified=False,
        error=None,
    )
    try:
        path = make_diffusion_root_path(values, sigma=float(sigma))
    except (ValueError, RuntimeError, np.linalg.LinAlgError) as error:
        result.update(error_phase="observation", error=f"{type(error).__name__}: {error}")
        return result
    n = len(path.values)
    positions = [2 * n - 2, *(child for child in path.observed_signature[-1] if child >= n)]
    for position, node in enumerate(positions):
        try:
            node_result = _node_chain_test(
                path, node, steps=steps, seed=seed + position, alpha=result["family_alpha"]
            )
        except ValueError as error:
            node_result = dict(
                node=int(node),
                status="unresolved",
                decision="unresolved",
                p_value=None,
                error=f"{type(error).__name__}: {error}",
            )
        node_result["role"] = "root" if position == 0 else "child"
        result["nodes"].append(node_result)
    result.update(
        status="complete",
        source_hashes={
            **path.source_hashes,
            str(Path(__file__).resolve()): hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        },
    )
    return result
