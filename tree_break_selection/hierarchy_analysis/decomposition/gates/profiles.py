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


# Settings a profile controls, with the value and type used when neither the
# caller nor the profile sets them.
_GATE_SETTING_DEFAULTS: dict[str, tuple[object, type]] = {
    "sibling_gate_method": ("projected_wald_inflation", str),
    "sibling_gate_alpha_penalty": (1.0, float),
    "root_stability_guard_threshold": (None, float),
    "root_stability_subsample_replicates": (0, int),
    "root_stability_feature_fraction": (0.8, float),
    "root_stability_seed": (0, int),
    "root_selective_permutation_guard_replicates": (0, int),
    "root_selective_permutation_guard_seed": (0, int),
    "root_selective_permutation_guard_alpha": (None, float),
    "root_selective_permutation_guard_scope": ("root", str),
}


def resolve_sibling_gate_profile_config(
    *,
    sibling_gate_profile: str | SiblingGateProfile | None = None,
    sibling_gate_method: str | None = None,
    sibling_gate_alpha_penalty: float | None = None,
    root_stability_guard_threshold: float | None = None,
    root_stability_subsample_replicates: int | None = None,
    root_stability_feature_fraction: float | None = None,
    root_stability_seed: int | None = None,
    root_selective_permutation_guard_replicates: int | None = None,
    root_selective_permutation_guard_seed: int | None = None,
    root_selective_permutation_guard_alpha: float | None = None,
    root_selective_permutation_guard_scope: str | None = None,
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
    """Resolve profile and explicit sibling-gate settings to concrete values.

    ``None`` means "not given": the profile value applies, else the default in
    ``_GATE_SETTING_DEFAULTS``. An explicit value the profile also sets must
    match it, even when it equals that default. A profile without a selective
    permutation guard leaves the ``root_selective_*`` settings to the caller.
    """
    explicit = {
        "sibling_gate_method": sibling_gate_method,
        "sibling_gate_alpha_penalty": sibling_gate_alpha_penalty,
        "root_stability_guard_threshold": root_stability_guard_threshold,
        "root_stability_subsample_replicates": root_stability_subsample_replicates,
        "root_stability_feature_fraction": root_stability_feature_fraction,
        "root_stability_seed": root_stability_seed,
        "root_selective_permutation_guard_replicates": (
            root_selective_permutation_guard_replicates
        ),
        "root_selective_permutation_guard_seed": root_selective_permutation_guard_seed,
        "root_selective_permutation_guard_alpha": root_selective_permutation_guard_alpha,
        "root_selective_permutation_guard_scope": root_selective_permutation_guard_scope,
    }
    profile = resolve_sibling_gate_profile(sibling_gate_profile)
    resolved: dict[str, object] = {}
    for name, (default, cast) in _GATE_SETTING_DEFAULTS.items():
        value = None if explicit[name] is None else cast(explicit[name])
        profile_sets = profile is not None and (
            not name.startswith("root_selective_")
            or profile.root_selective_permutation_guard_replicates > 0
        )
        if not profile_sets:
            resolved[name] = default if value is None else value
            continue
        profile_value = getattr(profile, name)
        if value is not None and value != profile_value:
            raise ValueError(
                f"sibling_gate_profile conflicts with explicit {name}: "
                f"profile has {profile_value!r}, explicit value is {value!r}."
            )
        resolved[name] = None if profile_value is None else cast(profile_value)
    return (None if profile is None else profile.profile_id, *resolved.values())  # type: ignore[return-value]


__all__ = [
    "GLOBAL_PASSTHROUGH_REFINED_REPLICATES",
    "SIBLING_GATE_PROFILES",
    "SiblingGateProfile",
    "resolve_sibling_gate_profile",
    "resolve_sibling_gate_profile_config",
]
