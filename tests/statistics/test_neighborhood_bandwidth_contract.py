from __future__ import annotations

import networkx as nx
import numpy as np
import pytest
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.neighborhood_bandwidth import (
    CoherentSupportDecision,
    SupportRole,
    TauRegionKey,
    build_branch_length_distance_cache,
    classify_support_role,
    estimate_tau_region,
    role_allows_empirical_null_calibration,
    selected_neighborhood_kernel_weights,
)


def test_branch_length_distance_cache_uses_weighted_paths() -> None:
    tree = nx.DiGraph()
    tree.add_edge("root", "left", branch_length=2.0)
    tree.add_edge("root", "right", branch_length=3.0)

    cache = build_branch_length_distance_cache(tree)

    assert cache.status == "branch_length_observed"
    assert cache.distance("left", "right") == 5.0
    assert cache.distance("left", "left") == 0.0


def test_branch_length_distance_cache_rejects_missing_lengths() -> None:
    tree = nx.DiGraph()
    tree.add_edge("root", "left", branch_length=2.0)
    tree.add_edge("root", "right")

    with pytest.raises(ValueError, match="requires explicit finite non-negative branch_length"):
        build_branch_length_distance_cache(tree)


def test_selected_signal_like_rows_do_not_calibrate_empirical_null() -> None:
    excluded = classify_support_role(
        {"topology_support_role": "selected_nonnull", "data_role": "signal"}
    )

    assert excluded == SupportRole.ALGORITHM_SELECTED_SIGNAL_LIKE_EXCLUDED
    assert not role_allows_empirical_null_calibration(excluded)
    assert role_allows_empirical_null_calibration(SupportRole.NULL_ANCHOR)
    assert role_allows_empirical_null_calibration(SupportRole.STOPPED_EDGE_NULL_ANCHOR)


def test_stopped_edge_role_takes_precedence_over_null_like_flag() -> None:
    role = classify_support_role(
        {
            "is_null_like": True,
            "is_edge_blocked": True,
        }
    )

    assert role == SupportRole.STOPPED_EDGE_NULL_ANCHOR


def test_invalid_explicit_support_role_fails_closed() -> None:
    with pytest.raises(ValueError, match="Invalid explicit support role"):
        classify_support_role(
            {
                "support_role": "typo_signal_anchor",
                "data_role": "null",
            }
        )

    with pytest.raises(ValueError, match="Invalid support role"):
        role_allows_empirical_null_calibration("typo_signal_anchor")


def test_tau_region_shrinks_sparse_regions_to_parent() -> None:
    parent_region = TauRegionKey("*", "*", "*", "*", "*", "*")
    child_region = TauRegionKey(
        "nonroot",
        "non_direct",
        "binary",
        "k_low",
        "frontier",
        "short",
    )
    parent = estimate_tau_region(
        region=parent_region,
        stopping_edge_distances=[4.0],
        stable_neighbor_distances=[6.0],
        signal_neighbor_distances=[8.0],
        stable_log_ks=[1.0, 2.0],
        n_support=5,
        n_signal=2,
        effective_support=5.0,
    )
    child = estimate_tau_region(
        region=child_region,
        stopping_edge_distances=[1.0],
        stable_neighbor_distances=[1.0],
        signal_neighbor_distances=[1.0],
        stable_log_ks=[1.0],
        n_support=1,
        n_signal=1,
        effective_support=1.0,
        parent_estimate=parent,
        shrinkage_strength=3.0,
    )

    assert parent.tau_status == "local_stable"
    assert child.tau_status == "borrowed_or_unstable"
    assert child.borrowed_from_region == parent_region
    assert child.tau_t == 0.25 * 1.0 + 0.75 * 6.0
    assert not child.stable_for_promotion


def test_kernel_weights_combine_tree_distance_and_log_dimension() -> None:
    weights = selected_neighborhood_kernel_weights(
        distances=np.array([0.0, 2.0]),
        source_log_k=np.log(np.array([2.0, 4.0])),
        target_log_k=float(np.log(2.0)),
        tau=2.0,
        h_k=1.0,
    )

    assert weights[0] == 1.0
    assert 0.0 < weights[1] < weights[0]


def test_coherent_support_decision_fails_closed_on_root_invalid() -> None:
    decision = CoherentSupportDecision(
        root_usable_or_nonroot=False,
        topology_coherent=True,
        regional_tau_stable=True,
        spectral_flow_supported=True,
        empirical_null_admissible_or_not_required=True,
    )

    assert not decision.promotion_eligible
    assert decision.dominant_blocker == "root_invalid_or_unusable"
    assert decision.method_action == "fail_closed_root_invalid"
