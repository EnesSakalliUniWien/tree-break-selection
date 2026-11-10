"""Opt-in subtree score references; these do not change TBS gates or traversal.

All inferential outputs test one common-mean Gaussian population in the tested
subtree. Topology differences do not test equality of reconstruction accuracy.
"""

from __future__ import annotations

from numbers import Integral

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist

from benchmarks.diagnostics.calibration.selected.hierarchy.subtree_gaussian_reference import (
    _inputs,
    subtree_gaussian_chain_test,
)
from benchmarks.diagnostics.calibration.selected.hierarchy.subtree_topology_geometry import (
    Cut,
    Tree,
    spectral_subspaces,
    unordered_cuts,
    ward_topology_family,
)
from benchmarks.shared.runners.tbs_runner import run_tbs_on_distance
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.pair_testing.collection.record_collection import (
    collect_sibling_pair_records,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.pair_testing.wald_statistic.sibling_divergence_test import (
    sibling_divergence_test,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.projection.gate_inputs.projection_dimensions import (
    derive_sibling_projection_dimensions_from_child_edge_comparisons,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.split_count_rule_out import (
    is_separable,
)
from tree_break_selection.tree.construction import build_tree
from tree_break_selection.tree.distributions import pooled_within_child_covariance_by_block
from tree_break_selection.tree.feature_space import continuous_feature_space_from_columns


def _frame(values):
    x = np.asarray(values, dtype=float)
    if x.ndim != 2 or len(x) < 2 or x.shape[1] < 1 or not np.isfinite(x).all():
        raise ValueError("Expected at least two finite rows and one feature")
    return pd.DataFrame(x, index=[f"S{i:012d}" for i in range(len(x))])


def sibling_cut_scores(values, cuts):
    """Existing full-space raw sibling Wald scores, zero branches, guarded covariance.

    This matches the earlier full-space sibling diagnostic, including parent
    covariance fallback. It excludes empirical inflation and adaptive rank.
    """
    data = _frame(values)
    x = data.to_numpy()
    space = continuous_feature_space_from_columns(data.columns)
    tree = build_tree(data, pdist(x), linkage_method="ward").tree
    tree.populate_node_distributions(data, feature_space=space)
    p = x.shape[1]
    scores = []
    for cut in cuts:
        if not cut.left or not cut.right or sorted((*cut.left, *cut.right)) != list(range(len(x))):
            raise ValueError("Each cut must partition the local row indices exactly once")
        covariance = pooled_within_child_covariance_by_block(
            tree, tree.root(), space, x, [list(cut.left), list(cut.right)]
        )
        statistic, _, _, _ = sibling_divergence_test(
            x[list(cut.left)].mean(0),
            x[list(cut.right)].mean(0),
            len(cut.left),
            len(cut.right),
            branch_length_sum=0.0,
            mean_branch_length=None,
            projection_dimension_from_edge_comparisons=p,
            parent_principal_component_projection=np.eye(p),
            parent_principal_component_eigenvalues=np.ones(p),
            feature_space=space,
            continuous_covariance_by_block=covariance,
        )
        if not np.isfinite(statistic):
            raise ValueError("Nonfinite sibling statistic")
        scores.append(statistic)
    return np.asarray(scores, dtype=float)


def _candidate_scores(points):
    spaces = spectral_subspaces(points)
    if not spaces:
        return (), np.array([], dtype=float)
    # Every smaller-span affine cut lies in the full cumulative span.
    coordinates = max(spaces, key=lambda space: space.rank).coordinates
    coordinates = coordinates / np.max(np.abs(coordinates), axis=0)
    cuts = tuple(
        cut
        for cut in unordered_cuts(len(points))
        if is_separable(coordinates, np.isin(np.arange(len(points)), cut.left))
    )
    return cuts, sibling_cut_scores(points, cuts)


def subtree_alternative_cut_test(
    values, *, subtree_indices, selection_event=None, steps=999, seed=0, alpha=0.05
):
    """Replay all full-span affine candidates and their maximum sibling score.

    Returned cut ranks share the maximum reference and give simultaneous control
    only under the complete common-mean subtree null. Membership is prespecified
    unless the supplied event reproduces its selection. No observed candidate
    count is substituted into a fixed-label Beta tail.
    """
    x, rows = _inputs(values, subtree_indices, alpha)
    observed = {}

    def score(matrix):
        cuts, scores = _candidate_scores(matrix[rows])
        if not observed:
            observed.update(cuts=cuts, scores=scores)
        return float(np.max(scores, initial=0.0))

    result = subtree_gaussian_chain_test(
        x,
        subtree_indices=rows,
        score=score,
        selection_event=selection_event,
        steps=steps,
        seed=seed,
        alpha=alpha,
    )
    result["candidate_family"] = "full_span_strict_affine_cuts"
    result["score"] = "full_space_raw_sibling_zero_branches"
    result["cut_results"] = [
        dict(
            left=rows[list(cut.left)].tolist(),
            right=rows[list(cut.right)].tolist(),
            statistic=float(statistic),
            global_null_adjusted_p_value=(
                sum(value >= statistic for value in result["chain_scores"]) / (steps + 1)
                if result["status"] == "complete"
                else None
            ),
        )
        for cut, statistic in zip(observed.get("cuts", ()), observed.get("scores", ()), strict=True)
    ]
    return result


def ward_subtree_selection(values, *, path=(), keep_split=False):
    """Capture a prespecified Ward node path, ordered by the smallest leaf ID.

    Return (row_indices, event). The event preserves every ancestral child pair;
    keep_split also preserves the tested node's own pair. The path must be chosen
    before observing the data. It does not encode an additional search for a node
    or a TBS gate/traversal decision. All outside rows stay fixed in the sampler.
    """
    path = _validated_ward_path(path)
    if not isinstance(keep_split, (bool, np.bool_)):
        raise ValueError("keep_split must be boolean")
    observed = _ward_description(values, path, keep_split=keep_split)
    if observed is None:
        raise ValueError("The specified path must reach an internal node")
    rows, _, _ = observed
    return np.asarray(
        rows
    ), lambda matrix: _ward_description(matrix, path, keep_split=keep_split) == observed


def _validated_ward_path(path):
    path = tuple(path)
    if any(
        isinstance(side, bool) or not isinstance(side, Integral) or side not in (0, 1)
        for side in path
    ):
        raise ValueError("path must contain only child positions 0 or 1")
    return path


def _ward_description(values, path, *, keep_split=False, keep_subtree=False):
    x = np.asarray(values, dtype=float)
    if x.ndim != 2 or len(x) < 2 or x.shape[1] < 1 or not np.isfinite(x).all():
        raise ValueError("Expected at least two finite rows and one feature")
    # Native build_tree calls this linkage on the same canonical positional IDs.
    rows = [(i,) for i in range(len(x))]
    children_by_node = [()] * len(x)
    for i, j in linkage(pdist(x), method="ward")[:, :2].astype(int):
        pair = tuple(sorted((i, j), key=rows.__getitem__))
        rows.append(tuple(sorted(rows[i] + rows[j])))
        children_by_node.append(pair)
    node, ancestors = len(rows) - 1, []
    for depth in range(len(path) + 1):
        pair = children_by_node[node]
        if len(pair) != 2:
            return None
        if depth < len(path) or keep_split:
            ancestors.append(tuple(rows[child] for child in pair))
        if depth < len(path):
            node = pair[path[depth]]

    family = []
    if keep_subtree:
        pending = [node]
        while pending:
            parent = pending.pop()
            pair = children_by_node[parent]
            if pair:
                if len(pair) != 2:
                    raise ValueError("Expected a binary Ward subtree")
                family.append((rows[parent], *(rows[child] for child in pair)))
                pending.extend(pair)
    return rows[node], tuple(ancestors), tuple(family)


def native_sibling_statistic(values, subtree_indices):
    """Refit Ward, NNLS, edge spectral context and native adaptive raw sibling score.

    This retains the native raw statistic, including branch time and .90
    contrast-energy rank. It does not reuse the native chi-square/inflated p-value,
    impose its eligibility gate, or claim control of the final partition.
    """
    _, _, _, _, _, record = _native_fit_record(values, subtree_indices)
    return float(record.stat)


def _native_fit_record(values, subtree_indices):
    data = _frame(values)
    _, rows = _inputs(values, subtree_indices, 0.05)
    space = continuous_feature_space_from_columns(data.columns)
    result = run_tbs_on_distance(
        data,
        pdist(data.to_numpy()),
        0.05,
        tree_linkage_method="ward",
        feature_space=space,
        adaptive_projection_dimension_energy_fraction=0.90,
        branch_length_optimization_method="fixed_topology_nnls",
        branch_length_data_df=data,
    )
    if result.status != "ok":
        raise ValueError(f"Native TBS replay failed: {result.status}")
    tree = result.extra["tree"]
    edge = result.extra["gate_bundle"].edge_gate_result
    spectral = edge.spectral_context
    dimensions = derive_sibling_projection_dimensions_from_child_edge_comparisons(
        tree, spectral_context=spectral
    )
    records, _ = collect_sibling_pair_records(
        tree,
        edge.annotated_df,
        sibling_projection_dimensions_from_edge_comparisons=dimensions,
        parent_principal_component_projections=spectral.principal_component_projections_by_node,
        parent_principal_component_eigenvalues=spectral.principal_component_eigenvalues_by_node,
        feature_space=space,
        adaptive_projection_dimension_energy_fraction=0.90,
    )
    target = frozenset(data.index[rows])
    clades = tree.compute_descendant_sets()
    record = next((record for record in records if clades[record.parent] == target), None)
    if record is None or not np.isfinite(record.stat):
        raise ValueError("Requested subtree has no finite native raw sibling score")
    return data, space, result, dimensions, spectral, record


def topology_sibling_score(values, topology: Tree):
    """Prespecified topology score: sum of raw full-space sibling scores at its nodes."""
    x = _frame(values).to_numpy()

    def leaves(node):
        if isinstance(node, int):
            return [node]
        if not isinstance(node, tuple) or len(node) != 2:
            raise ValueError("Expected a binary nested-tuple topology")
        return leaves(node[0]) + leaves(node[1])

    if sorted(leaves(topology)) != list(range(len(x))):
        raise ValueError("Topology must contain each local row index exactly once")

    def score(node):
        if isinstance(node, int):
            return 0.0
        left, right = leaves(node[0]), leaves(node[1])
        cut = Cut(tuple(range(len(left))), tuple(range(len(left), len(left) + len(right))))
        local = sibling_cut_scores(x[left + right], (cut,))[0]
        return float(local + score(node[0]) + score(node[1]))

    return score(topology)


def topology_pair_reference(
    values,
    *,
    subtree_indices,
    first: Tree,
    second: Tree,
    selection_event=None,
    steps=999,
    seed=0,
    alpha=0.05,
):
    """Paired absolute score difference for two prespecified labeled topologies.

    A data-selected pair needs its pair-selection event as well. The null is
    homogeneous Gaussian rows, not equal expected accuracy of the two trees.
    """
    x, rows = _inputs(values, subtree_indices, alpha)

    def score(matrix):
        return abs(
            topology_sibling_score(matrix[rows], first)
            - topology_sibling_score(matrix[rows], second)
        )

    result = subtree_gaussian_chain_test(
        x,
        subtree_indices=rows,
        score=score,
        selection_event=selection_event,
        steps=steps,
        seed=seed,
        alpha=alpha,
    )
    return result | {"score": "absolute_difference_of_summed_full_space_sibling_scores"}


def subtree_topology_search_test(
    values, *, subtree_indices, selection_event=None, steps=999, seed=0, alpha=0.05
):
    """Replay the union of tied-Ward families over whole spectral spans (at most six rows).

    The geometry rule fixes the relative tie tolerance at 1e-10. Counts are
    descriptive; the maximum is recalculated over the new family on every draw.
    Scores sum raw full-space sibling scores, using the same feature space for
    every candidate. They are a declared diagnostic criterion, not a likelihood.
    """
    x, rows = _inputs(values, subtree_indices, alpha)
    if len(rows) > 6:
        raise ValueError("Tied-Ward topology search requires at most six subtree leaves")

    def candidates(matrix):
        points = matrix[rows]
        trees = set()
        for space in spectral_subspaces(points):
            trees.update(ward_topology_family(space.coordinates, rtol=1e-10))
        ordered = sorted(trees, key=repr)
        return ordered, [topology_sibling_score(points, tree) for tree in ordered]

    observed = {}

    def score(matrix):
        trees, scores = candidates(matrix)
        if not observed:
            observed.update(trees=trees, scores=scores)
        return max(scores, default=0.0)

    result = subtree_gaussian_chain_test(
        x,
        subtree_indices=rows,
        score=score,
        selection_event=selection_event,
        steps=steps,
        seed=seed,
        alpha=alpha,
    )
    result.update(candidate_family="union_of_spectral_tied_ward_topologies", tie_rtol=1e-10)
    result["topology_results"] = [
        dict(
            topology=tree,
            statistic=statistic,
            global_null_adjusted_p_value=(
                sum(value >= statistic for value in result["chain_scores"]) / (steps + 1)
                if result["status"] == "complete"
                else None
            ),
        )
        for tree, statistic in zip(
            observed.get("trees", ()), observed.get("scores", ()), strict=True
        )
    ]
    return result
