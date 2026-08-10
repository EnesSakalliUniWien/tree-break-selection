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
    resolved = resolve_sibling_gate_profile_config(sibling_gate_profile=profile_id)

    assert resolved[:11] == (
        profile_id,
        "fixed_coordinate_bh",
        50.0,
        0.24,
        12,
        0.8,
        0,
        root_selective_replicates,
        0,
        root_selective_alpha,
        root_selective_scope,
    )
    assert SIBLING_GATE_PROFILES[profile_id].status == "diagnostic_candidate"


def test_fixed_profile_allows_explicit_root_selective_guard_when_disabled() -> None:
    (
        _profile_id,
        _method,
        _penalty,
        _threshold,
        _stability_replicates,
        _fraction,
        _stability_seed,
        root_selective_replicates,
        root_selective_seed,
        root_selective_alpha,
        root_selective_scope,
    ) = resolve_sibling_gate_profile_config(
        sibling_gate_profile="fixed_coordinate_guarded_v1",
        root_selective_permutation_guard_replicates=5,
        root_selective_permutation_guard_seed=123,
        root_selective_permutation_guard_alpha=0.01,
    )

    assert root_selective_replicates == 5
    assert root_selective_seed == 123
    assert root_selective_alpha == 0.01
    assert root_selective_scope == "root"


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
