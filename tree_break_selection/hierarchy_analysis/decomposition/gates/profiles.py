"""Named sibling-gate profile contracts.

Profiles are reusable strategy presets for sibling-gate annotation. Keeping the
registry here separates auditable method choices from the annotation pipeline
that executes them.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class SiblingGateProfile:
    """Named sibling-gate configuration for auditable method candidates."""

    profile_id: str
    sibling_gate_method: str
    sibling_gate_alpha_penalty: float
    root_stability_guard_threshold: float | None
    root_stability_subsample_replicates: int
    root_stability_feature_fraction: float
    root_stability_seed: int
    status: str
    description: str
    root_selective_permutation_guard_replicates: int = 0
    root_selective_permutation_guard_seed: int = 0
    root_selective_permutation_guard_alpha: float | None = None
    root_selective_permutation_guard_scope: str = "root"


GLOBAL_PASSTHROUGH_REFINED_REPLICATES = 999


SIBLING_GATE_PROFILES: dict[str, SiblingGateProfile] = {
    "fixed_coordinate_guarded_v1": SiblingGateProfile(
        profile_id="fixed_coordinate_guarded_v1",
        sibling_gate_method="fixed_coordinate_bh",
        sibling_gate_alpha_penalty=50.0,
        root_stability_guard_threshold=0.24,
        root_stability_subsample_replicates=12,
        root_stability_feature_fraction=0.8,
        root_stability_seed=0,
        status="diagnostic_candidate",
        description=(
            "Same-data fixed coordinate-BH sibling gate with selected-topology "
            "alpha penalty and selected-root feature-subsample stability guard."
        ),
    ),
    "fixed_coordinate_global_passthrough_refined_v1": SiblingGateProfile(
        profile_id="fixed_coordinate_global_passthrough_refined_v1",
        sibling_gate_method="fixed_coordinate_bh",
        sibling_gate_alpha_penalty=50.0,
        root_stability_guard_threshold=0.24,
        root_stability_subsample_replicates=12,
        root_stability_feature_fraction=0.8,
        root_stability_seed=0,
        status="diagnostic_candidate",
        description=(
            "Same-data fixed coordinate-BH sibling gate with selected-topology "
            "alpha penalty, selected-root feature-subsample stability guard, and "
            "a global selected-family sibling-min permutation guard that reruns "
            "Monte Carlo floor pass-through families at higher resolution."
        ),
        root_selective_permutation_guard_replicates=99,
        root_selective_permutation_guard_seed=0,
        root_selective_permutation_guard_alpha=0.01,
        root_selective_permutation_guard_scope=(
            "global_sibling_min_passthrough_descendant_refined"
        ),
    ),
}


def resolve_sibling_gate_profile(
    sibling_gate_profile: str | SiblingGateProfile | None,
) -> SiblingGateProfile | None:
    """Resolve a named or explicit sibling-gate profile."""
    if sibling_gate_profile is None:
        return None
    if isinstance(sibling_gate_profile, SiblingGateProfile):
        return sibling_gate_profile
    profile_id = str(sibling_gate_profile)
    try:
        return SIBLING_GATE_PROFILES[profile_id]
    except KeyError as exc:
        raise ValueError(
            f"Unknown sibling_gate_profile {profile_id!r}; "
            f"allowed={tuple(sorted(SIBLING_GATE_PROFILES))!r}."
        ) from exc


def _profile_value(
    *,
    field_name: str,
    current: object,
    default: object,
    profile_value: object,
) -> object:
    if current == default or current == profile_value:
        return profile_value
    raise ValueError(
        f"sibling_gate_profile conflicts with explicit {field_name}: "
        f"profile has {profile_value!r}, explicit value is {current!r}."
    )


def resolve_sibling_gate_profile_config(
    *,
    sibling_gate_profile: str | SiblingGateProfile | None = None,
    sibling_gate_method: str = "projected_wald_inflation",
    sibling_gate_alpha_penalty: float = 1.0,
    root_stability_guard_threshold: float | None = None,
    root_stability_subsample_replicates: int = 0,
    root_stability_feature_fraction: float = 0.8,
    root_stability_seed: int = 0,
    root_selective_permutation_guard_replicates: int = 0,
    root_selective_permutation_guard_seed: int = 0,
    root_selective_permutation_guard_alpha: float | None = None,
    root_selective_permutation_guard_scope: str = "root",
) -> tuple[
    str | None,
    str,
    float,
    float | None,
    int,
    float,
    int,
    int,
    int,
    float | None,
    str,
]:
    """Resolve profile and explicit sibling-gate settings to concrete values."""
    profile = resolve_sibling_gate_profile(sibling_gate_profile)
    if profile is None:
        return (
            None,
            str(sibling_gate_method),
            float(sibling_gate_alpha_penalty),
            (
                None
                if root_stability_guard_threshold is None
                else float(root_stability_guard_threshold)
            ),
            int(root_stability_subsample_replicates),
            float(root_stability_feature_fraction),
            int(root_stability_seed),
            int(root_selective_permutation_guard_replicates),
            int(root_selective_permutation_guard_seed),
            (
                None
                if root_selective_permutation_guard_alpha is None
                else float(root_selective_permutation_guard_alpha)
            ),
            str(root_selective_permutation_guard_scope),
        )

    method = _profile_value(
        field_name="sibling_gate_method",
        current=str(sibling_gate_method),
        default="projected_wald_inflation",
        profile_value=profile.sibling_gate_method,
    )
    penalty = _profile_value(
        field_name="sibling_gate_alpha_penalty",
        current=float(sibling_gate_alpha_penalty),
        default=1.0,
        profile_value=float(profile.sibling_gate_alpha_penalty),
    )
    threshold = _profile_value(
        field_name="root_stability_guard_threshold",
        current=(
            None
            if root_stability_guard_threshold is None
            else float(root_stability_guard_threshold)
        ),
        default=None,
        profile_value=profile.root_stability_guard_threshold,
    )
    replicates = _profile_value(
        field_name="root_stability_subsample_replicates",
        current=int(root_stability_subsample_replicates),
        default=0,
        profile_value=int(profile.root_stability_subsample_replicates),
    )
    fraction = _profile_value(
        field_name="root_stability_feature_fraction",
        current=float(root_stability_feature_fraction),
        default=0.8,
        profile_value=float(profile.root_stability_feature_fraction),
    )
    seed = _profile_value(
        field_name="root_stability_seed",
        current=int(root_stability_seed),
        default=0,
        profile_value=int(profile.root_stability_seed),
    )
    if int(profile.root_selective_permutation_guard_replicates) > 0:
        root_selective_replicates = _profile_value(
            field_name="root_selective_permutation_guard_replicates",
            current=int(root_selective_permutation_guard_replicates),
            default=0,
            profile_value=int(profile.root_selective_permutation_guard_replicates),
        )
        root_selective_seed = _profile_value(
            field_name="root_selective_permutation_guard_seed",
            current=int(root_selective_permutation_guard_seed),
            default=0,
            profile_value=int(profile.root_selective_permutation_guard_seed),
        )
        root_selective_alpha = _profile_value(
            field_name="root_selective_permutation_guard_alpha",
            current=(
                None
                if root_selective_permutation_guard_alpha is None
                else float(root_selective_permutation_guard_alpha)
            ),
            default=None,
            profile_value=profile.root_selective_permutation_guard_alpha,
        )
        root_selective_scope = _profile_value(
            field_name="root_selective_permutation_guard_scope",
            current=str(root_selective_permutation_guard_scope),
            default="root",
            profile_value=str(profile.root_selective_permutation_guard_scope),
        )
    else:
        root_selective_replicates = int(root_selective_permutation_guard_replicates)
        root_selective_seed = int(root_selective_permutation_guard_seed)
        root_selective_alpha = (
            None
            if root_selective_permutation_guard_alpha is None
            else float(root_selective_permutation_guard_alpha)
        )
        root_selective_scope = str(root_selective_permutation_guard_scope)
    return (
        profile.profile_id,
        str(method),
        float(penalty),
        None if threshold is None else float(threshold),
        int(replicates),
        float(fraction),
        int(seed),
        int(root_selective_replicates),
        int(root_selective_seed),
        None if root_selective_alpha is None else float(root_selective_alpha),
        str(root_selective_scope),
    )


__all__ = [
    "GLOBAL_PASSTHROUGH_REFINED_REPLICATES",
    "SIBLING_GATE_PROFILES",
    "SiblingGateProfile",
    "resolve_sibling_gate_profile",
    "resolve_sibling_gate_profile_config",
]
