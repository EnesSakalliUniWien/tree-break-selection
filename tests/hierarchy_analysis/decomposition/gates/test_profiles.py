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
        pytest.param("fixed_coordinate_selective_root_v1", 99, 0.01, "root", id="selective-root"),
        pytest.param(
            "fixed_coordinate_selective_passthrough_v1",
            99,
            0.01,
            "passthrough_descendant",
            id="selective-passthrough",
        ),
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


def test_fixed_profile_resolves_spectral_transport_passthrough_constants() -> None:
    (
        profile_id,
        method,
        penalty,
        threshold,
        stability_replicates,
        fraction,
        stability_seed,
        root_selective_replicates,
        root_selective_seed,
        root_selective_alpha,
        root_selective_scope,
        spectral_guard,
        spectral_max_cost,
        spectral_require_mp_blocks,
        spectral_block_log_tolerance,
        spectral_unmatched_mode_penalty,
    ) = resolve_sibling_gate_profile_config(
        sibling_gate_profile="fixed_coordinate_spectral_transport_passthrough_v1",
    )

    assert profile_id == "fixed_coordinate_spectral_transport_passthrough_v1"
    assert method == "fixed_coordinate_bh"
    assert penalty == 50.0
    assert threshold == 0.24
    assert stability_replicates == 12
    assert fraction == 0.8
    assert stability_seed == 0
    assert root_selective_replicates == 99
    assert root_selective_seed == 0
    assert root_selective_alpha == 0.01
    assert root_selective_scope == ("global_sibling_min_passthrough_descendant_refined")
    assert spectral_guard is True
    assert spectral_max_cost == 1.2
    assert spectral_require_mp_blocks is True
    assert spectral_block_log_tolerance == 0.05
    assert spectral_unmatched_mode_penalty == 1.0
    assert SIBLING_GATE_PROFILES[profile_id].status == "opt_in_candidate_not_default"


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
        *_spectral_transport_config,
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


def test_selective_root_profile_rejects_conflicting_guard_replicates() -> None:
    with pytest.raises(
        ValueError,
        match="conflicts with explicit root_selective_permutation_guard_replicates",
    ):
        resolve_sibling_gate_profile_config(
            sibling_gate_profile="fixed_coordinate_selective_root_v1",
            root_selective_permutation_guard_replicates=5,
        )


def test_fixed_profile_rejects_conflicting_explicit_method() -> None:
    with pytest.raises(ValueError, match="conflicts with explicit sibling_gate_method"):
        resolve_sibling_gate_profile_config(
            sibling_gate_profile="fixed_coordinate_guarded_v1",
            sibling_gate_method="fixed_block_bh",
        )
