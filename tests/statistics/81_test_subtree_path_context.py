"""Paths must distinguish structural absence from surviving nodes off the path."""

from math import comb

import numpy as np
import pytest

from benchmarks.diagnostics.calibration.selected.hierarchy.subtree_path_context import (
    compare_surviving_paths,
)


def path_row(reference, alternative, start, end):
    return next(
        row
        for row in compare_surviving_paths(reference, alternative)
        if row.start == start and row.end == end
    )


def test_identical_tree_includes_every_pair_and_both_relations():
    tree = ((0, 1), (2, 3))
    rows = list(compare_surviving_paths(tree, tree))
    assert len(rows) == comb(7, 2)
    assert {r.relation for r in rows} == {"ancestor_descendant", "cross_branch"}
    for row in rows:
        assert row.reference_path == row.alternative_path
        assert not row.missing_intermediates
        assert not row.surviving_intermediates_off_path
        assert not row.inserted_intermediates
        assert not row.missing_neighborhood
        assert not row.surviving_neighborhood_displaced


def test_missing_chain_and_neighbor_clades_between_surviving_endpoints():
    reference = ((((0, 1), (2, 3)), (4, 5)), (6, 7))
    alternative = ((0, 1), ((2, 4), (3, (5, (6, 7)))))
    row = path_row(reference, alternative, tuple(range(8)), (0, 1))
    assert row.relation == "ancestor_descendant"
    assert row.alternative_path == (tuple(range(8)), (0, 1))
    assert row.missing_intermediates == (tuple(range(6)), tuple(range(4)))
    assert row.missing_neighborhood == ((2, 3), (4, 5))
    assert (6, 7) in row.surviving_neighborhood_displaced
    assert not row.surviving_intermediates_off_path


def test_changed_lca_can_survive_outside_the_new_path():
    reference = (((0, 1), (2, 3)), (4, 5))
    alternative = (((0, 2), (1, 3)), (4, 5))
    row = path_row(reference, alternative, (0,), (2,))
    assert row.relation == "cross_branch"
    assert row.reference_lca == (0, 1, 2, 3)
    assert row.alternative_lca == (0, 2)
    assert row.missing_intermediates == ((0, 1), (2, 3))
    assert row.surviving_intermediates_off_path == ((0, 1, 2, 3),)
    assert row.inserted_intermediates == ((0, 2),)


def test_child_order_does_not_change_path_comparisons():
    tree = (((0, 1), 2), 3)
    reversed_tree = (3, (2, (1, 0)))
    assert list(compare_surviving_paths(tree, tree)) == list(
        compare_surviving_paths(tree, reversed_tree)
    )


@pytest.mark.parametrize("alternative", [((0, 1), 3), ((0, 1), 1), ((0, 1), True)])
def test_invalid_leaf_universe_is_rejected(alternative):
    with pytest.raises(ValueError):
        list(compare_surviving_paths(((0, 1), 2), alternative))


def test_full_mean_path_telescopes_but_standardized_squares_do_not_add():
    path = [(0,), (0, 1), (0, 1, 2, 3), (2, 3), (2,)]
    weights = np.zeros((len(path), 4))
    for row, clade in enumerate(path):
        weights[row, list(clade)] = 1 / len(clade)
    increments = np.diff(weights, axis=0)
    gram = increments @ increments.T
    values = np.array([0, 2, 3, 7])
    np.testing.assert_allclose(increments.sum(axis=0), weights[-1] - weights[0])
    assert gram[1, 2] == pytest.approx(0.25)
    assert np.sum((increments @ values) ** 2 / np.diag(gram)) == pytest.approx(42)
    endpoint = weights[-1] - weights[0]
    assert (endpoint @ values) ** 2 / (endpoint @ endpoint) == pytest.approx(4.5)


def test_saved_path_profiles_reconstruct_actual_masks_and_lengths():
    from reports.subtree_context_statistics_20261109 import run as native
    from reports.subtree_path_statistics_20261109.run import SavedFit, row_for

    values, truth, _ = native.dataset("crossed", 0)
    _, z, coordinates, basis, _ = next(native.proposals(values))
    result = native.fit_tree(values, truth[2], z, coordinates, basis)
    saved = {"linkage": z, **result}
    fit = SavedFit(values, saved)
    lengths = {(p, c): length for p, c, length in result["trace"]["branches"]}
    for child, parent in fit.parents.items():
        if parent is None:
            continue
        expected = values[list(child)[::2]].mean(axis=0) - fit.parent_means[parent]
        np.testing.assert_allclose(
            fit.potential[child] - fit.potential[parent], expected, atol=1e-12
        )
        assert fit.path_metrics((parent, child), parent)["half_continuity_error"] < 1e-12
    checked = set()
    for context in compare_surviving_paths(fit.tree, fit.tree):
        if context.relation in checked or len(context.reference_path) < 3:
            continue
        path = context.reference_path
        expected_length = sum(
            lengths.get((fit.ids[u], fit.ids[v]), lengths.get((fit.ids[v], fit.ids[u])))
            for u, v in zip(path, path[1:])
        )
        metric = fit.path_metrics(path, context.reference_lca)
        assert metric["branch_sum"] == pytest.approx(expected_length)
        increments = []
        for u, v in zip(path, path[1:]):
            increments.append(fit.input_mean(v, u) - fit.input_mean(u, v))
        np.testing.assert_allclose(
            np.sum(increments, axis=0), fit.potential[path[-1]] - fit.potential[path[0]], atol=1e-12
        )
        row = row_for(context, fit, fit)
        assert row["same_path"] and row["same_neighborhood"]
        assert row["half_path_increment_change"] == 0
        checked.add(context.relation)
        if len(checked) == 2:
            break
    assert checked == {"ancestor_descendant", "cross_branch"}
