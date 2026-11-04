"""Root-only chain ranks, preserved selection, and explicit unresolved decisions."""

import numpy as np
import pandas as pd
import pytest

from benchmarks.diagnostics.calibration.selected.hierarchy import (
    diffusion_root_conditional_mc as mc,
)


def test_random_anchor_runs_both_directions_and_counts_upper_ties(monkeypatch):
    class FixedDraws:
        def __init__(self):
            self.draws = iter([-0.2, 0.5, -0.1, 0.1])

        def integers(self, high):
            assert high == 5
            return 2

        def normal(self):
            return next(self.draws)

        def random(self):
            return 1e-12

    monkeypatch.setattr(mc.np.random, "default_rng", lambda seed: FixedDraws())
    result = mc._sample_chi_chain(2.0, 1, lambda r: True, steps=4, seed=1)
    np.testing.assert_allclose(result["chain_radii"], [2.3, 1.8, 2.0, 1.9, 2.0])
    assert result["anchor_index"] == 2
    assert result["p_value"] == 3 / 5
    assert result["accepted_moves"] == 4


def test_stalled_chain_has_rank_one_and_no_movement():
    result = mc._sample_chi_chain(2.0, 3, lambda r: r == 2.0, steps=39, seed=5)
    assert result["status"] == "complete"
    assert result["p_value"] == 1
    assert result["accepted_moves"] == 0
    assert result["unique_states"] == 1
    assert result["sampling_status"] == "stalled"


def test_event_failure_never_becomes_a_rejected_proposal():
    def broken_event(radius):
        raise RuntimeError("geometry failed")

    result = mc._sample_chi_chain(2.0, 3, broken_event, steps=99, seed=5)
    assert result["status"] == "unresolved"
    assert result["p_value"] is None
    assert "geometry failed" in result["error"]
    assert result["failed_radius"] > 0
    assert result["completed_steps"] < 99


def test_root_partition_ignores_commuting_merges_but_keeps_membership():
    first = ((0, 1), (2, 3), (4, 5))
    reordered = ((2, 3), (0, 1), (4, 5))
    changed = ((0, 2), (1, 3), (4, 5))
    assert mc._root_partition(first, 4) == ((0, 1), (2, 3))
    assert mc._root_partition(reordered, 4) == ((0, 1), (2, 3))
    assert mc._root_partition(changed, 4) != ((0, 1), (2, 3))


@pytest.mark.parametrize(
    "kwargs",
    [
        {"steps": 0},
        {"steps": -1},
        {"steps": True},
        {"steps": 1.5},
        {"seed": -1},
        {"seed": True},
        {"seed": 1.5},
        {"alpha": 0},
        {"alpha": 1},
        {"alpha": np.nan},
        {"alpha": True},
        {"sigma": 0},
        {"sigma": -1},
        {"sigma": np.nan},
        {"sigma": True},
    ],
)
def test_invalid_options_fail_before_geometry(kwargs):
    options = dict(sigma=1, steps=3, seed=0, alpha=0.05)
    options.update(kwargs)
    with pytest.raises(ValueError):
        mc.selected_diffusion_root_chain_test([[0], [1], [2]], **options)


@pytest.mark.parametrize("values", [[], [[np.nan], [1]], [[np.inf], [1]], [1, 2]])
def test_invalid_observations_raise(values):
    with pytest.raises(ValueError):
        mc.selected_diffusion_root_chain_test(values, sigma=1, steps=3)


def test_known_sigma_is_required():
    with pytest.raises(TypeError):
        mc.selected_diffusion_root_chain_test([[0], [1], [2]])


def test_native_chain_preserves_labels_data_and_observed_root():
    values = pd.DataFrame(
        np.random.default_rng(20261105).normal(size=(9, 1)), index=[f"S{i}" for i in range(9)]
    )
    original = values.copy(deep=True)
    first = mc.selected_diffusion_root_chain_test(values, sigma=1, steps=59, seed=17)
    second = mc.selected_diffusion_root_chain_test(
        values.sample(frac=1, random_state=8), sigma=1, steps=59, seed=17
    )
    assert first == second
    pd.testing.assert_frame_equal(values, original)
    assert first["status"] == "complete" and first["selection_event"] == "root_split"
    assert len(first["chain_radii"]) == 60
    assert first["source_hashes"] and first["minimum_p_value"] == 1 / 60
    assert first["numerical_event_certified"] is False
    assert set(sum(first["root_children"], ())) == set(values.index)
    assert first["decision"] in {"reject", "do_not_reject", "unresolved"}


def test_native_strong_case_matches_saved_candidate_and_decision_threshold():
    values = np.random.default_rng(2026110401).normal(size=(9, 3))
    values[:4, 0] -= 2
    values[4:, 0] += 2
    result = mc.selected_diffusion_root_chain_test(
        values, sigma=1, steps=999, seed=2026140401, alpha=0.01
    )
    assert result["status"] == "complete" and result["decision"] == "reject"
    assert result["p_value"] == 0.001
    assert result["accepted_moves"] == 144


def test_unsupported_observation_returns_no_decision():
    result = mc.selected_diffusion_root_chain_test(np.zeros((9, 1)), sigma=1, steps=3)
    assert result["status"] == "unresolved" and result["decision"] == "unresolved"
    assert result["p_value"] is None and result["error_phase"] == "observation"


def test_native_proposal_failure_returns_no_partial_decision(monkeypatch):
    values = np.random.default_rng(20261105).normal(size=(9, 1))

    def broken_geometry(*args):
        raise RuntimeError("native proposal failure")

    monkeypatch.setattr(mc, "_geometry", broken_geometry, raising=False)
    result = mc.selected_diffusion_root_chain_test(values, sigma=1, steps=99, seed=17)
    assert result["status"] == "unresolved" and result["decision"] == "unresolved"
    assert result["p_value"] is None and result["error_phase"] == "sampling"
    assert "native proposal failure" in result["error"]
