"""Deterministic spectral alignment and reversible diffusion checks."""

import numpy as np
import pytest
from scipy import sparse

from benchmarks.diagnostics.calibration.selected.hierarchy.subtree_spectral_alignment import (
    compare_parent_child,
    compare_reflection,
    compare_subspaces,
    covariance_decomposition,
    covariance_spectrum,
    orthonormal_basis,
    restrict_parent_basis,
    reversible_diffusion_spectrum,
)
from tree_break_selection.space_separation.diffusion import compute_diffusion_coordinates


def test_sign_basis_rotation_and_unequal_rank():
    a = np.eye(4)[:, :2]
    b = a @ np.array([[0.6, 0.8], [0.8, -0.6]])
    same = compare_subspaces(a, b)
    assert same.same_span
    np.testing.assert_allclose(same.singular_overlaps, [1, 1], atol=1e-14)
    np.testing.assert_allclose(same.principal_angles, [0, 0], atol=1e-14)
    assert same.projector_frobenius < 1e-14
    unequal = compare_subspaces(a[:, :1], b)
    assert (unequal.rank_a, unequal.rank_b) == (1, 2)
    assert unequal.a_contained_in_b and not unequal.b_contained_in_a and not unequal.same_span
    assert unequal.projector_frobenius == pytest.approx(1)
    np.testing.assert_allclose(unequal.principal_angles, [0], atol=1e-14)


def test_zero_and_dependent_spans_are_explicit():
    x = np.array([[1, 2], [2, 4], [0, 0]], dtype=float)
    q = orthonormal_basis(x)
    assert q.shape == (3, 1)
    assert compare_subspaces(x, 1e-250 * x).same_span
    empty = np.empty((3, 0))
    result = compare_subspaces(empty, x)
    assert (result.rank_a, result.rank_b) == (0, 1)
    assert result.principal_angles.size == 0
    assert result.a_contained_in_b and not result.b_contained_in_a
    assert result.projector_frobenius == pytest.approx(1)
    assert compare_subspaces(empty, np.zeros((3, 2))).same_span


def test_equal_covariance_eigenvalues_can_have_different_eigenspaces():
    a = np.array([[2, 0, 0], [-2, 0, 0], [0, 2, 0], [0, -2, 0], [0, 0, 1], [0, 0, -1]])
    b = a[:, [2, 1, 0]]
    sa, sb = covariance_spectrum(a), covariance_spectrum(b)
    np.testing.assert_allclose(sa.eigenvalues, sb.eigenvalues)
    assert sa.blocks == sb.blocks == ((0, 2), (2, 3))
    result = compare_subspaces(sa.eigenvectors[:, :2], sb.eigenvectors[:, :2])
    np.testing.assert_allclose(result.singular_overlaps, [1, 0], atol=1e-14)
    assert result.projector_frobenius == pytest.approx(np.sqrt(2))
    assert not result.same_span


def test_prescribed_reflection_preserves_or_changes_separate_objects():
    reflection = np.diag([-1.0, 1])
    line = np.array([[1.0], [1]])
    covariance = np.array([[3.0, 1], [1, 2]])
    broken = compare_reflection(
        line, line, reflection, covariance_a=covariance, covariance_b=covariance
    )
    assert not broken.subspaces.same_span and not broken.covariance_matches
    reflected = compare_reflection(
        line,
        reflection @ line,
        reflection,
        covariance_a=covariance,
        covariance_b=reflection @ covariance @ reflection.T,
    )
    assert reflected.subspaces.same_span and reflected.covariance_matches
    axis = np.array([[1.0], [0]])
    assert compare_reflection(axis, axis, reflection).subspaces.same_span
    assert compare_reflection(axis, axis, reflection).covariance_matches is None


@pytest.mark.parametrize("ddof", [0, 1])
def test_covariance_decomposition_weights_child_scatter_and_mean_shifts(ddof):
    x = np.array([[0, 0], [2, 1], [4, 3], [10, 0], [14, -2], [21, 2]], dtype=float)
    groups = [0, 0, 0, 1, 1, 2]
    parts = covariance_decomposition(x, groups, ddof=ddof)
    mean = x.mean(axis=0)
    expected_within = np.zeros((2, 2))
    expected_between = np.zeros((2, 2))
    for label in np.unique(groups):
        child = x[np.asarray(groups) == label]
        centered = child - child.mean(axis=0)
        expected_within += centered.T @ centered
        diff = child.mean(axis=0) - mean
        expected_between += len(child) * np.outer(diff, diff)
    np.testing.assert_allclose(parts.within, expected_within / (len(x) - ddof))
    np.testing.assert_allclose(parts.between, expected_between / (len(x) - ddof))
    np.testing.assert_allclose(parts.total, np.cov(x, rowvar=False, ddof=ddof))
    np.testing.assert_allclose(parts.total, parts.within + parts.between)


def test_irregular_diffusion_bases_need_degree_transform():
    weights = np.array([[9, 1, 3, 0], [0, 8, 2, 1], [0, 2, 7, 1], [0, 1, 0, 6]], dtype=float)
    result = reversible_diffusion_spectrum(sparse.csr_matrix(weights))
    expected = np.maximum(weights, weights.T)
    np.fill_diagonal(expected, 0)
    np.testing.assert_array_equal(result.weights, expected)
    np.testing.assert_allclose(result.markov_operator.sum(axis=1), 1)
    assert result.eigenvalues[-1] < 0
    for error in result.residuals.values():
        assert error < 1e-12
    leading = result.block_alignment(result.blocks[0])
    assert not leading["symmetric_vs_right"].same_span
    assert not leading["right_vs_left"].same_span
    assert leading["symmetric_vs_degree_right"].same_span
    assert leading["symmetric_vs_inverse_degree_left"].same_span
    for operator in ("symmetric", "markov"):
        fitted = compute_diffusion_coordinates(
            weights, diffusion_time=2, n_components=3, coordinate_operator=operator
        )
        source = result.symmetric_vectors if operator == "symmetric" else result.right_vectors
        retained = source[:, 1:] * np.maximum(result.eigenvalues[1:], 0) ** 2
        assert compare_subspaces(fitted, retained).same_span


def test_row_restriction_detects_rank_loss_and_keeps_child_order():
    parent = np.eye(5)[:, :3]
    child_rows = [4, 1, 0]
    restricted = restrict_parent_basis(parent, child_rows)
    assert (restricted.parent_rank, restricted.restricted_rank, restricted.rank_loss) == (3, 2, 1)
    child = np.eye(3)[:, 1:]
    assert compare_parent_child(parent, child, child_rows).same_span
    assert restrict_parent_basis(parent, [3, 4]).restricted_rank == 0


@pytest.mark.parametrize(
    "call",
    [
        lambda: compare_subspaces(np.ones((2, 1)), np.ones((3, 1))),
        lambda: orthonormal_basis([[np.nan]]),
        lambda: orthonormal_basis([[1j]]),
        lambda: orthonormal_basis([[1]], rtol=-1),
        lambda: covariance_decomposition([[1], [2]], [0]),
        lambda: covariance_spectrum([[1]], ddof=1),
        lambda: reversible_diffusion_spectrum([[0, -1], [1, 0]]),
        lambda: reversible_diffusion_spectrum([[0, 0], [0, 0]]),
        lambda: restrict_parent_basis(np.eye(3), [0, 0]),
        lambda: restrict_parent_basis(np.eye(3), [1.0]),
        lambda: compare_reflection(np.eye(2), np.eye(2), [[0, -1], [1, 0]]),
        lambda: compare_reflection(np.eye(2), np.eye(2), np.eye(2), covariance_a=np.eye(2)),
    ],
)
def test_invalid_inputs_fail_explicitly(call):
    with pytest.raises(ValueError):
        call()


def test_diffusion_keeps_whole_zero_block_and_negative_mode():
    w = np.zeros((4, 4))
    w[0, 1:] = [1, 2, 3]
    w[1:, 0] = w[0, 1:]
    result = reversible_diffusion_spectrum(w)
    assert result.blocks == ((0, 1), (1, 3), (3, 4))
    np.testing.assert_allclose(result.eigenvalues, [1, 0, 0, -1], atol=1e-14)
    for block in result.blocks:
        assert result.block_alignment(block)["symmetric_vs_degree_right"].same_span
    with pytest.raises(ValueError, match="whole"):
        result.block_alignment((1, 2))


def test_rebuilt_child_modes_can_disagree_with_the_restricted_parent_modes():
    w = np.array([[0, 1, 3], [1, 0, 1], [3, 1, 0]], dtype=float)
    parent = reversible_diffusion_spectrum(w)
    child_rows = [0, 1]
    child = reversible_diffusion_spectrum(w[np.ix_(child_rows, child_rows)])
    result = compare_parent_child(
        parent.symmetric_vectors[:, :1], child.symmetric_vectors[:, :1], child_rows
    )
    assert not result.same_span
    aligned = compare_parent_child(
        parent.right_vectors[:, :1], child.right_vectors[:, :1], child_rows
    )
    assert aligned.same_span
