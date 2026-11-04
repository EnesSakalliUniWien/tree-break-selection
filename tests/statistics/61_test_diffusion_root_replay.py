"""Observable actual-map replay, without interval certification."""

import numpy as np
import pandas as pd
import pytest

from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_replay import (
    make_diffusion_root_path,
)


def test_actual_diffusion_reentry_and_point_margins():
    x = np.random.default_rng(20261105).normal(size=(9, 1))
    path = make_diffusion_root_path(x)
    np.testing.assert_allclose(path.offset + path.observed_radius * path.slope, x, atol=1e-14)
    np.testing.assert_allclose(path.e @ path.offset, 0.0, atol=1e-14)
    for fraction, status in [(1.0, "selected"), (1.1, "nonselected"), (1.2, "selected")]:
        point = path.evaluate(path.observed_radius * fraction)
        assert point.status == status
        assert len(point.margins) == sum(k * (k - 1) // 2 - 1 for k in range(2, 10))
        assert (min(point.margins) > 0) == (status == "selected")
        assert point.state["kde_dimension"] == 1
        assert point.state["retained_dimension"] == 1
        assert point.state["neighbors"].shape == (9, 8)
        assert point.state["eigenvalues"].shape == (9,)
    assert path.observed_signature == (
        (2, 6),
        (3, 5),
        (7, 8),
        (4, 10),
        (1, 9),
        (0, 13),
        (11, 14),
        (12, 15),
    )
    # Reentry is full-history membership; these commuting merges preserve clades.
    assert (
        path.evaluate(path.observed_radius).topology_signature
        == path.evaluate(1.1 * path.observed_radius).topology_signature
    )
    assert path.source_hashes


def test_labels_canonical_and_caller_data_unchanged():
    x = pd.DataFrame(
        np.random.default_rng(20261105).normal(size=(9, 1)), index=[f"S{i}" for i in range(9)]
    )
    original = x.copy(deep=True)
    path = make_diffusion_root_path(x.sample(frac=1.0, random_state=2))
    expected = make_diffusion_root_path(x)
    assert path.labels == expected.labels
    np.testing.assert_allclose(path.offset, expected.offset)
    pd.testing.assert_frame_equal(x, original)


def test_unsupported_point_and_invalid_radius_are_explicit():
    path = make_diffusion_root_path(np.random.default_rng(20261105).normal(size=(9, 1)))
    for radius in [-1.0, np.nan, np.inf]:
        with pytest.raises(ValueError):
            path.evaluate(radius)
    with pytest.raises(ValueError, match="observed"):
        make_diffusion_root_path(np.zeros((9, 1)))


def test_observed_ties_are_not_treated_as_selected():
    with pytest.raises(ValueError, match="tie"):
        make_diffusion_root_path(np.arange(9.0)[:, None])


def test_failed_radius_replay_is_not_treated_as_nonselected(monkeypatch):
    from benchmarks.diagnostics.calibration.selected.hierarchy import diffusion_root_replay

    path = make_diffusion_root_path(np.random.default_rng(20261105).normal(size=(9, 1)))

    def fail(*args, **kwargs):
        raise ValueError("geometry unavailable at this radius")

    monkeypatch.setattr(diffusion_root_replay, "_geometry", fail)
    point = path.evaluate(path.observed_radius)
    assert point.status == "unsupported"
    assert "geometry unavailable" in point.error
    assert point.margins.size == 0
