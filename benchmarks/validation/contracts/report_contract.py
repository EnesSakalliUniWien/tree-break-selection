"""Shared validation-report contract helpers."""

from __future__ import annotations

import math
import subprocess
from datetime import datetime, timezone


def wilson_interval(
    successes: int,
    total: int,
    *,
    z_value: float = 1.959963984540054,
) -> tuple[float, float]:
    """Return the bounded Wilson score interval for a binomial proportion."""

    if total <= 0:
        raise ValueError("total must be positive.")
    proportion = successes / total
    denominator = 1.0 + z_value**2 / total
    center = (proportion + z_value**2 / (2.0 * total)) / denominator
    half_width = (
        z_value
        * math.sqrt((proportion * (1.0 - proportion) + z_value**2 / (4.0 * total)) / total)
        / denominator
    )
    return max(0.0, center - half_width), min(1.0, center + half_width)


def read_git_commit() -> str:
    """Return the checked-out Git commit."""

    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    )
    commit = result.stdout.strip()
    if not commit:
        raise RuntimeError("git rev-parse HEAD returned an empty commit.")
    return commit


def read_git_worktree_status() -> list[str]:
    """Return non-empty porcelain status lines for the current worktree."""

    result = subprocess.run(
        ["git", "status", "--short"],
        check=True,
        capture_output=True,
        text=True,
    )
    return [line for line in result.stdout.splitlines() if line]


def utc_now() -> str:
    """Return the current UTC timestamp in second-resolution ISO format."""

    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


__all__ = [
    "read_git_commit",
    "read_git_worktree_status",
    "utc_now",
    "wilson_interval",
]
