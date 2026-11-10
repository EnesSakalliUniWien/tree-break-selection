"""Bounded subtree candidate selection and all-pair objective checks."""

import numpy as np
import pandas as pd
import pytest
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist

from benchmarks.diagnostics.calibration.selected.hierarchy.subtree_topology_candidate import (
    MAX_CANDIDATE_LEAVES,
    select_subtree_topology,
    topology_from_linkage,
)
from benchmarks.diagnostics.calibration.selected.hierarchy.subtree_topology_geometry import (
    fit_additive_distances,
    topology_linkage,
    tree_clades,
)
from tree_break_selection.tree.optimized_branch_lengths import (
    BRANCH_LENGTH_TARGET_MEAN_SQUARED_STANDARDIZED_EUCLIDEAN,
    _mean_squared_standardized_targets,
)


@pytest.mark.parametrize(
    "metric", ["euclidean", BRANCH_LENGTH_TARGET_MEAN_SQUARED_STANDARDIZED_EUCLIDEAN]
)
def test_compressed_all_pair_objective_matches_dense_nnls(metric):
    x = np.random.default_rng(138).normal(size=(6, 3))
    x[:, 2] = 0.0
    current = ((0, 1), ((2, 3), (4, 5)))
    result = select_subtree_topology(
        x, current_topology=current, target_metric=metric, exhaustive=True
    )
    assert result.coverage == "exhaustive_rooted" and result.candidate_count == 945
    assert result.pair_count == 15
    if metric == "euclidean":
        target = pdist(x)
    else:
        target = _mean_squared_standardized_targets(pd.DataFrame(x), *np.triu_indices(len(x), k=1))
    for candidate in result.candidates:
        expected = fit_additive_distances(candidate.topology, target)
        np.testing.assert_allclose(
            candidate.fit.branch_lengths, expected.branch_lengths, atol=1e-12
        )
        np.testing.assert_allclose(
            candidate.fit.fitted_distances, expected.fitted_distances, atol=1e-12
        )
        assert candidate.fit.rmse == pytest.approx(expected.rmse, abs=1e-12)


def test_simplex_equivalence_retains_current_linkage_and_exposes_zero_edges():
    x = np.eye(4)
    current = ((0, 1), (2, 3))
    z = topology_linkage(current) * np.array([1, 1, 2.5, 1])
    result = select_subtree_topology(
        x, current_topology=current, current_linkage=z, exhaustive=True
    )
    assert result.selected_topology == current
    np.testing.assert_array_equal(result.selected_linkage, z)
    assert result.spectral_span_count == 1
    assert all(candidate.zero_internal_edge_count == 1 for candidate in result.candidates)
    assert all(candidate.equivalent_to_selected_metric for candidate in result.candidates)


def test_rectangle_rotation_and_row_relabel_preserve_labeled_geometry():
    x = np.array([[-3, -1], [-3, 1], [3, -1], [3, 1]], dtype=float)
    square = x / [3, 1]
    current = ((0, 2), (1, 3))
    wanted = ((0, 1), (2, 3))
    a = select_subtree_topology(x, current_topology=current)
    assert a.selected_topology == wanted
    assert a.spectral_span_count == 3
    assert topology_from_linkage(a.selected_linkage) == wanted
    assert a.selected_linkage[-1, 2] == pytest.approx(np.sqrt(2) * 6)
    rotation = np.array([[0.6, -0.8], [0.8, 0.6]])
    first = select_subtree_topology(square, current_topology=current)
    second = select_subtree_topology(square @ rotation, current_topology=current)
    assert first.spectral_span_count == second.spectral_span_count == 1
    assert first.selected_topology == second.selected_topology
    assert [c.topology for c in first.candidates] == [c.topology for c in second.candidates]

    order = np.array([2, 3, 1, 0])
    b = select_subtree_topology(x[order], current_topology=((0, 3), (1, 2)), leaf_labels=order)
    relabeled = {frozenset(order[list(c)]) for c in tree_clades(b.selected_topology)}
    assert relabeled == tree_clades(a.selected_topology)
    np.testing.assert_allclose(b.selected_linkage[:, 2], a.selected_linkage[:, 2])


def test_rerootings_share_the_same_unrooted_fit():
    x = np.array([[0.0], [1.0], [8.0], [10.0]])
    current = ((0, 1), (2, 3))
    rerooted = (0, (1, (2, 3)))
    result = select_subtree_topology(x, current_topology=current, exhaustive=True)
    fits = {candidate.topology: candidate.fit for candidate in result.candidates}
    assert fits[current] is fits[rerooted]


def test_bounded_candidate_support_and_explicit_exhaustive_limit():
    x = np.random.default_rng(15).normal(size=(32, 5))
    z = linkage(x, method="ward")
    current = topology_from_linkage(z)
    result = select_subtree_topology(x, current_topology=current, current_linkage=z)
    assert result.coverage == "bounded_spectral_ward"
    assert result.pair_count == 32 * 31 // 2
    assert result.spectral_span_count == 9
    assert 1 <= result.candidate_count <= result.spectral_span_count + 1
    with pytest.raises(ValueError, match="Exhaustive"):
        select_subtree_topology(x, current_topology=current, exhaustive=True)
    with pytest.raises(ValueError, match=str(MAX_CANDIDATE_LEAVES)):
        select_subtree_topology(np.ones((MAX_CANDIDATE_LEAVES + 1, 1)), current_topology=current)


@pytest.mark.parametrize(
    ("row", "column", "value", "message"),
    [
        (0, 2, np.nan, "finite"),
        (0, 0, 0.25, "integer"),
        (0, 3, 2.5, "integer"),
        (0, 3, 3, "count"),
        (0, 2, 5, "height"),
    ],
)
def test_invalid_current_linkage_is_rejected(row, column, value, message):
    current = ((0, 1), (2, 3))
    z = topology_linkage(current)
    z[row, column] = value
    with pytest.raises(ValueError, match=message):
        select_subtree_topology(
            np.arange(4.0)[:, None], current_topology=current, current_linkage=z
        )


@pytest.mark.parametrize(
    "current", [((0, 1), 2), ((0, 1), (2, 2)), ((0, 1), (2, 4)), ((False, 1), (2, 3))]
)
def test_current_topology_requires_every_local_leaf_exactly_once(current):
    with pytest.raises(ValueError, match="local leaf"):
        select_subtree_topology(np.arange(4.0)[:, None], current_topology=current)
