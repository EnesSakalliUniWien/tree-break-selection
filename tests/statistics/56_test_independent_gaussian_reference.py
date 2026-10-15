"""Independent mathematical and data-provenance tests for Gaussian discovery gates."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy.integrate import quad
from scipy.stats import chi2

from tree_break_selection.hierarchy_analysis.statistics.independent_gaussian import (
    fit_independent_gaussian_tree,
    test_independent_gaussian_tree,
)
from tree_break_selection.hierarchy_analysis.statistics.projection.projected_wald.energy_selected_reference import (
    energy_selected_gaussian_sf,
)


@pytest.mark.parametrize("q", [0.1, 0.5, 0.9, 0.9999])
@pytest.mark.parametrize("t", [0.0, 1.0, 7.0, 40.0])
def test_two_dimensional_tail_against_independent_angular_integral(q, t):
    bound = np.arccos(np.sqrt(q))
    expected = (
        quad(lambda theta: np.exp(-t / (2 * np.cos(theta) ** 2)), 0, bound, epsabs=1e-13)[0] / bound
    )
    actual = energy_selected_gaussian_sf(
        t,
        candidate_dimension=2,
        projection_dimension=1,
        previous_group_dimension=0,
        energy_fraction=q,
    )
    assert actual == pytest.approx(expected, rel=2e-9, abs=2e-12)


def test_wolfram_checkpoint():
    assert energy_selected_gaussian_sf(
        7,
        candidate_dimension=2,
        projection_dimension=1,
        previous_group_dimension=0,
        energy_fraction=0.9,
    ) == pytest.approx(
        0.02679366463202124077335995555965,
        rel=1e-10,
    )


@pytest.mark.parametrize("dimensions", [(8, 3, 2), (8, 4, 2), (6, 2, 0)])
def test_general_group_law_against_independent_gaussian_sampling(dimensions):
    m, k, previous = dimensions
    rng = np.random.default_rng(8231)
    z2 = rng.normal(size=(300000, m)) ** 2
    shares = z2.cumsum(1) / z2.sum(1)[:, None]
    q = 0.5
    event = shares[:, k - 1] >= q
    if previous:
        event &= shares[:, previous - 1] < q
    retained = z2[event, :k].sum(1)
    assert len(retained) > 10000
    for t in (4, 8, 12):
        exact = energy_selected_gaussian_sf(
            t,
            candidate_dimension=m,
            projection_dimension=k,
            previous_group_dimension=previous,
            energy_fraction=q,
        )
        empirical = np.mean(retained >= t)
        se = np.sqrt(exact * (1 - exact) / len(retained))
        assert abs(empirical - exact) < 6 * se + 1 / len(retained)


def test_full_group_and_fraction_one_keep_radial_law():
    for previous in (0, 3):
        assert energy_selected_gaussian_sf(
            9,
            candidate_dimension=4,
            projection_dimension=4,
            previous_group_dimension=previous,
            energy_fraction=1,
        ) == pytest.approx(chi2.sf(9, 4))
    with pytest.raises(ValueError, match="probability zero"):
        energy_selected_gaussian_sf(
            9,
            candidate_dimension=4,
            projection_dimension=2,
            previous_group_dimension=1,
            energy_fraction=1,
        )


def matrices(dimension=3):
    rng = np.random.default_rng(792)
    labels = [f"S{i}" for i in range(12)]
    discovery = pd.DataFrame(rng.normal(size=(12, dimension)), index=labels)
    inference = pd.DataFrame(rng.normal(size=(12, dimension)), index=labels)
    return discovery, inference


def fitted(discovery):
    return fit_independent_gaussian_tree(
        discovery,
        known_covariance=np.eye(discovery.shape[1]),
        observation_ids=[f"d{i}" for i in range(len(discovery))],
    )


def test_independent_measurements_feed_actual_gate_and_traversal():
    discovery, inference = matrices()
    model = fitted(discovery)
    original_projections = [x.projection.copy() for x in model.nodes]
    result = test_independent_gaussian_tree(
        model, inference, observation_ids=[f"t{i}" for i in range(12)]
    )
    assert len(result.tests) == 11
    assert result.tests.p_value.between(0, 1).all()
    assert result.tests.reference_law.eq("energy_selected_gaussian_conditional").all()
    np.testing.assert_allclose(
        result.annotations.loc[result.tests.index, "Sibling_Divergence_P_Value"],
        result.tests.p_value,
    )
    assert "traversal_trace" in result.decomposition
    for node, original in zip(model.nodes, original_projections, strict=True):
        np.testing.assert_array_equal(node.projection, original)
        assert not node.projection.flags.writeable
    # A second inference matrix never refits the discovery tree or directions.
    test_independent_gaussian_tree(
        model, inference * 2, observation_ids=[f"u{i}" for i in range(12)]
    )
    assert frozenset(model.tree.edges) == model.tree_edges


def test_rejects_reuse_alignment_and_invalid_covariance():
    discovery, inference = matrices()
    model = fitted(discovery)
    with pytest.raises(ValueError, match="disjoint"):
        test_independent_gaussian_tree(
            model, inference, observation_ids=[f"d{i}" for i in range(12)]
        )
    with pytest.raises(ValueError, match="reuses"):
        test_independent_gaussian_tree(
            model, discovery.copy(), observation_ids=[f"t{i}" for i in range(12)]
        )
    with pytest.raises(ValueError, match="match discovery"):
        test_independent_gaussian_tree(model, inference.iloc[::-1], observation_ids=range(12))
    with pytest.raises(ValueError, match="positive definite"):
        fit_independent_gaussian_tree(
            discovery, known_covariance=np.zeros((3, 3)), observation_ids=range(12)
        )


def test_dual_pca_independent_fixed_rank_keeps_exact_chi_square():
    discovery, inference = matrices(dimension=20)
    model = fit_independent_gaussian_tree(
        discovery, known_covariance=np.eye(20), observation_ids=range(12), energy_fraction=None
    )
    result = test_independent_gaussian_tree(model, inference, observation_ids=range(12, 24))
    np.testing.assert_allclose(
        result.tests.p_value, chi2.sf(result.tests.statistic, result.tests.projection_dimension)
    )
