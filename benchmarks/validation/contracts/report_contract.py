"""Shared validation-report contract helpers."""

from __future__ import annotations

import subprocess
from datetime import datetime, timezone

from scipy.stats import binomtest


def wilson_interval(
    successes: int,
    total: int,
    *,
    confidence_level: float = 0.95,
) -> tuple[float, float]:
    """Return the two-sided Wilson score interval for a binomial proportion."""

    if total <= 0:
        raise ValueError("total must be positive.")
    interval = binomtest(int(successes), int(total)).proportion_ci(
        confidence_level=confidence_level, method="wilson"
    )
    return float(interval.low), float(interval.high)


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
