"""Common-mean conditional Gaussian reference checks."""

import numpy as np
import pytest
from scipy.linalg import helmert
from scipy.stats import beta

from benchmarks.diagnostics.calibration.selected.hierarchy.subtree_gaussian_reference import (
    fixed_subtree_contrast_test,
    subtree_gaussian_chain_test,
)


def test_chain_preserves_subtree_mean_gram_and_outside_rows_at_full_row_rank():
    x = np.random.default_rng(4).normal(size=(10, 12))
    rows = [8, 1, 3, 4, 7]
    outside = sorted(set(range(len(x))) - set(rows))
    centered = x[rows] - x[rows].mean(0)

    def score(values):
        np.testing.assert_array_equal(values[outside], x[outside])
        np.testing.assert_allclose(values[rows].mean(0), x[rows].mean(0), atol=1e-14)
        y = values[rows] - values[rows].mean(0)
        np.testing.assert_allclose(y.T @ y, centered.T @ centered, atol=1e-13)
        return values[rows[0], 0]

    r = subtree_gaussian_chain_test(x, subtree_indices=rows, score=score, steps=39, seed=17)
    assert r["status"] == "complete"
    assert r["accepted_moves"] == r["completed_steps"] == 39
    assert r["chain_scores"][r["anchor_index"]] == score(x)
    assert r["p_value"] == np.mean(np.asarray(r["chain_scores"]) >= score(x))
    assert r["anchor_index"] == np.random.default_rng(17).integers(40)
    assert len(r["chain_moves"]) == 40
    assert {move[3] for move in r["chain_moves"] if move is not None} == {-1, 1}
    sorted_rows = r["subtree_indices"]
    sorted_centered = x[sorted_rows] - x[sorted_rows].mean(0)
    basis = helmert(len(rows), full=False)
    anchor = r["anchor_index"]
    for positions in (range(anchor - 1, -1, -1), range(anchor + 1, 40)):
        state = basis @ sorted_centered
        for position in positions:
            i, j, angle, sign = r["chain_moves"][position]
            a, b = state[[i, j]]
            c, sine = np.cos(angle), np.sin(angle)
            state[i] = c * a - sine * b
            state[j] = sign * (sine * a + c * b)
            reconstructed = x[sorted_rows].mean(0) + basis.T @ state
            assert r["chain_scores"][position] == pytest.approx(
                reconstructed[sorted_rows.index(rows[0]), 0]
            )
    assert r["selection_event"] == "fixed_membership"
    assert r == subtree_gaussian_chain_test(x, subtree_indices=rows, score=score, steps=39, seed=17)


def test_selection_rejections_preserve_the_last_state_and_do_not_evaluate_its_score():
    x = np.arange(12, dtype=float).reshape(6, 2)
    rows = [0, 1, 2, 3]

    def event(values):
        return values[0, 0] <= 0.2

    def score(values):
        assert event(values)
        return values[0, 0]

    r = subtree_gaussian_chain_test(
        x, subtree_indices=rows, score=score, selection_event=event, steps=59, seed=18
    )
    assert r["status"] == "complete"
    assert r["selection_event"] == "caller_supplied"
    assert 0 < r["accepted_moves"] < 59
    assert max(r["chain_scores"]) <= 0.2
    for i, (move, value) in enumerate(zip(r["chain_moves"], r["chain_scores"], strict=True)):
        if i != r["anchor_index"] and move is None:
            assert value == r["chain_scores"][i + (1 if i < r["anchor_index"] else -1)]


def test_two_row_orbit_moves_by_reflection_and_upper_tail_ties_are_retained():
    x = np.array([[1.0], [3.0], [99.0]])
    r = subtree_gaussian_chain_test(
        x, subtree_indices=[0, 1], score=lambda _: 1.0, steps=59, seed=11
    )
    assert r["decision"] == "do_not_reject" and r["p_value"] == 1.0
    assert 0 < r["accepted_moves"] < 59
    assert r["chain_scores"] == [1.0] * 60


@pytest.mark.parametrize(
    ("event", "score", "phase"),
    [
        (lambda _: False, lambda _: 1.0, "observation"),
        (lambda _: 1, lambda _: 1.0, "observation"),
        (lambda _: True, lambda _: np.nan, "observation"),
        (lambda x: True if x[0, 0] == 0 else 1, lambda _: 1.0, "sampling"),
        (lambda _: True, lambda x: 0.0 if x[0, 0] == 0 else np.inf, "sampling"),
    ],
)
def test_callback_failures_are_explicit(event, score, phase):
    r = subtree_gaussian_chain_test(
        np.arange(12).reshape(6, 2),
        subtree_indices=[0, 1, 2, 3],
        score=score,
        selection_event=event,
        steps=59,
        seed=3,
    )
    assert r["status"] == r["decision"] == "unresolved"
    assert r["p_value"] is None and r["error_phase"] == phase


def test_stall_and_insufficient_rank_resolution_are_distinct():
    x = np.arange(12).reshape(6, 2)
    stalled = subtree_gaussian_chain_test(
        x,
        subtree_indices=[0, 1, 2],
        score=lambda _: 0.0,
        selection_event=lambda z: np.array_equal(x, z),
        steps=39,
    )
    assert stalled["decision"] == "unresolved"
    assert stalled["sampling_status"] == "stalled" and stalled["p_value"] == 1.0
    moved = subtree_gaussian_chain_test(x, subtree_indices=[0, 1, 2], score=lambda _: 0.0, steps=1)
    assert moved["decision"] == "unresolved"
    assert moved["sampling_status"] == "moved" and not moved["resolution_sufficient"]


@pytest.mark.parametrize("indices", [[], [0], [0, 0], [-1, 0], [0, 6], [True, 2], [0.5, 1]])
def test_invalid_subtree_indices_raise(indices):
    with pytest.raises(ValueError, match="subtree_indices"):
        subtree_gaussian_chain_test(np.ones((6, 2)), subtree_indices=indices, score=lambda _: 0.0)


def test_beta_statistic_uses_unit_row_contrast_and_all_labels_count():
    x = np.random.default_rng(12).normal(size=(10, 3))
    rows, left = [1, 2, 3, 5, 8, 9], [2, 3]
    a = fixed_subtree_contrast_test(x, subtree_indices=rows, left_indices=left)
    b = fixed_subtree_contrast_test(
        x, subtree_indices=rows, left_indices=[1, 5, 8, 9], all_labels=True
    )
    u, _, _ = np.linalg.svd(x[rows] - x[rows].mean(0), full_matrices=False)
    e = np.array([-0.25, 0.5, 0.5, -0.25, -0.25, -0.25])
    q = np.linalg.norm(u.T @ e) ** 2 / (e @ e)
    assert a["statistic"] == pytest.approx(q)
    assert a["p_value"] == pytest.approx(beta.sf(q, 1.5, 1.0))
    assert b["number_of_labels"] == 31
    assert b["p_value"] == pytest.approx(min(1.0, 31 * a["p_value"]))


def test_beta_gram_rank_and_degenerate_support():
    rng = np.random.default_rng(14)
    x = rng.normal(size=(6, 10))
    r = fixed_subtree_contrast_test(x, subtree_indices=range(6), left_indices=[0, 1])
    assert r["decision"] == "unresolved" and r["p_value"] is None
    r = fixed_subtree_contrast_test(x, subtree_indices=range(6), left_indices=[0, 1], rank=2)
    assert r["status"] == "complete" and r["rank"] == 2
    r = fixed_subtree_contrast_test(np.ones((6, 1)), subtree_indices=range(6), left_indices=[0])
    assert r["decision"] == "unresolved"
    with pytest.raises(ValueError, match="left_indices"):
        fixed_subtree_contrast_test(x, subtree_indices=[0, 1, 2], left_indices=[5])


@pytest.mark.parametrize(
    "options",
    [
        {"steps": 0},
        {"steps": True},
        {"seed": -1},
        {"alpha": 1.0},
        {"alpha": np.nan},
        {"selection_event": 2},
    ],
)
def test_invalid_chain_options_raise(options):
    with pytest.raises(ValueError):
        subtree_gaussian_chain_test(
            np.ones((6, 2)), subtree_indices=[0, 1, 2], score=lambda _: 0.0, **options
        )


def test_beta_leading_rank_cannot_cut_a_tied_gram_eigenspace():
    x = helmert(4, full=False).T[:, :2]
    r = fixed_subtree_contrast_test(x, subtree_indices=range(4), left_indices=[0], rank=1)
    assert r["decision"] == "unresolved" and "tied Gram eigenspace" in r["error"]


def test_subtree_membership_order_cannot_change_the_chain_kernel():
    x = np.random.default_rng(5).normal(size=(8, 4))

    def score(values):
        return values[0, 0]

    forward = subtree_gaussian_chain_test(x, subtree_indices=[0, 2, 3, 5], score=score, steps=19)
    reverse = subtree_gaussian_chain_test(x, subtree_indices=[5, 3, 2, 0], score=score, steps=19)
    assert forward == reverse
