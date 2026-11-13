"""Conditional familywise homogeneity reference for a selected Ward subtree."""

import math
from numbers import Integral

import numpy as np

from .subtree_gaussian_reference import subtree_gaussian_chain_test
from .subtree_topology_geometry import Cut
from .subtree_topology_tests import _validated_ward_path, _ward_description, sibling_cut_scores


def _fixed_split_score(matrix, rows, left, right):
    # Under two-row Gram conditioning the only orientations swap the rows.
    if len(rows) == 2:
        return 0.0
    count = len(left)
    cut = Cut(tuple(range(count)), tuple(range(count, len(rows))))
    return float(sibling_cut_scores(matrix[np.r_[left, right]], (cut,))[0])


def ward_subtree_family_test(values, *, path=(), score=None, steps=999, seed=0, alpha=0.05):
    """Test every internal node under one shared Ward selection event.

    The predetermined path selects the subtree. The common event preserves
    its ancestral child pairs and the whole labeled topology inside it. For
    each node, freeze all outside rows and rotate only that node. The true
    null at a tested node requires independent rows with one common Gaussian
    mean and covariance, independent of outside rows. A custom deterministic
    score receives (matrix, node rows, left rows, right rows); all row arrays
    use input IDs. Holm controls false homogeneity rejections under that
    model and event, including configurations with false nulls.
    """
    path = _validated_ward_path(path)
    if score is not None and not callable(score):
        raise ValueError("score must be callable")
    score = _fixed_split_score if score is None else score
    if isinstance(seed, bool) or not isinstance(seed, Integral) or seed < 0:
        raise ValueError("seed must be a nonnegative integer")
    observed = _ward_description(values, path, keep_subtree=True)
    if observed is None:
        return dict(
            status="unresolved",
            decision="unresolved",
            nodes=[],
            error="The specified path must reach an internal node",
        )
    subtree, _, family = observed

    def event(matrix):
        return _ward_description(matrix, path, keep_subtree=True) == observed

    nodes = []
    for i, (node, left, right) in enumerate(family):
        rows, a, b = map(np.asarray, (node, left, right))
        for indices in (rows, a, b):
            indices.setflags(write=False)
        result = subtree_gaussian_chain_test(
            values,
            subtree_indices=rows,
            score=lambda x: score(x, rows, a, b),
            selection_event=event,
            steps=steps,
            seed=seed + i,
            alpha=alpha,
        )
        result.update(left=list(left), right=list(right))
        nodes.append(result)

    m = len(nodes)
    resolution = 1 / (steps + 1) <= alpha / m
    order = sorted(
        range(m), key=lambda j: nodes[j]["p_value"] if nodes[j]["decision"] != "unresolved" else 1
    )
    running = 0.0
    for i, j in enumerate(order):
        node = nodes[j]
        p_value = node["p_value"] if node["decision"] != "unresolved" else 1.0
        running = max(running, min(1.0, (m - i) * p_value))
        node["holm_p_value"] = running
        node["family_decision"] = (
            "unresolved"
            if not resolution or node["decision"] == "unresolved"
            else "reject"
            if running <= alpha
            else "do_not_reject"
        )
    decisions = [node["family_decision"] for node in nodes]
    return dict(
        status="complete",
        decision="reject"
        if "reject" in decisions
        else "unresolved"
        if "unresolved" in decisions
        else "do_not_reject",
        unresolved_nodes=decisions.count("unresolved"),
        subtree_indices=list(subtree),
        hypotheses=m,
        nodes=nodes,
        path=list(path),
        selection_event="common_ward_ancestry_and_subtree_topology",
        alpha=float(alpha),
        resolution_sufficient=resolution,
        holm_first_cutoff=alpha / m,
        minimum_p_value=1 / (steps + 1),
        minimum_steps_for_holm=math.ceil(m / alpha) - 1,
    )
