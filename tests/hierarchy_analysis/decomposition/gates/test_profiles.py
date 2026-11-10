"""Contracts for resolving gate profiles."""

from __future__ import annotations

import pytest

from tree_break_selection.hierarchy_analysis.decomposition.gates.profiles import (
    SIBLING_GATE_PROFILES,
    resolve_sibling_gate_profile_config,
)


@pytest.mark.parametrize(
    (
        "profile_id",
        "root_selective_replicates",
        "root_selective_alpha",
        "root_selective_scope",
    ),
    [
        pytest.param("fixed_coordinate_guarded_v1", 0, None, "root", id="guarded"),
        pytest.param(
            "fixed_coordinate_global_passthrough_refined_v1",
            99,
            0.01,
            "global_sibling_min_passthrough_descendant_refined",
            id="refined-global-passthrough",
        ),
    ],
)
def test_fixed_profile_resolves_candidate_constants(
    profile_id: str,
    root_selective_replicates: int,
    root_selective_alpha: float | None,
    root_selective_scope: str,
) -> None:
    resolved_id, settings = resolve_sibling_gate_profile_config(sibling_gate_profile=profile_id)

    assert resolved_id == profile_id
    assert settings == {
        "sibling_gate_method": "fixed_coordinate_bh",
        "sibling_gate_alpha_penalty": 50.0,
        "root_stability_guard_threshold": 0.24,
        "root_stability_subsample_replicates": 12,
        "root_stability_feature_fraction": 0.8,
        "root_stability_seed": 0,
        "root_selective_permutation_guard_replicates": root_selective_replicates,
        "root_selective_permutation_guard_seed": 0,
        "root_selective_permutation_guard_alpha": root_selective_alpha,
        "root_selective_permutation_guard_scope": root_selective_scope,
    }
    assert SIBLING_GATE_PROFILES[profile_id].status == "diagnostic_candidate"


def test_fixed_profile_allows_explicit_root_selective_guard_when_disabled() -> None:
    _profile_id, settings = resolve_sibling_gate_profile_config(
        sibling_gate_profile="fixed_coordinate_guarded_v1",
        root_selective_permutation_guard_replicates=5,
        root_selective_permutation_guard_seed=123,
        root_selective_permutation_guard_alpha=0.01,
    )

    assert settings["root_selective_permutation_guard_replicates"] == 5
    assert settings["root_selective_permutation_guard_seed"] == 123
    assert settings["root_selective_permutation_guard_alpha"] == 0.01
    assert settings["root_selective_permutation_guard_scope"] == "root"


def test_profile_rejects_conflicting_guard_replicates() -> None:
    with pytest.raises(
        ValueError,
        match="conflicts with explicit root_selective_permutation_guard_replicates",
    ):
        resolve_sibling_gate_profile_config(
            sibling_gate_profile="fixed_coordinate_global_passthrough_refined_v1",
            root_selective_permutation_guard_replicates=5,
        )


def test_fixed_profile_rejects_conflicting_explicit_method() -> None:
    with pytest.raises(ValueError, match="conflicts with explicit sibling_gate_method"):
        resolve_sibling_gate_profile_config(
            sibling_gate_profile="fixed_coordinate_guarded_v1",
            sibling_gate_method="fixed_block_bh",
        )


def test_fixed_profile_rejects_explicit_default_method() -> None:
    # An explicit value is a conflict even when it equals the unprofiled default.
    with pytest.raises(ValueError, match="conflicts with explicit sibling_gate_method"):
        resolve_sibling_gate_profile_config(
            sibling_gate_profile="fixed_coordinate_guarded_v1",
            sibling_gate_method="projected_wald_inflation",
        )
