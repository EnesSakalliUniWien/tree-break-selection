"""Git provenance recorded in benchmark manifests."""

from __future__ import annotations

import os
import subprocess


def current_git_state() -> dict[str, object]:
    """Return git provenance without hiding command failures."""
    state: dict[str, object] = {
        "build_commit": os.environ.get("TBS_GIT_COMMIT", "unknown"),
        "build_branch": os.environ.get("TBS_GIT_BRANCH", "unknown"),
        "build_dirty": os.environ.get("TBS_GIT_DIRTY", "unknown"),
    }
    commands = {
        "commit": ("git", "rev-parse", "HEAD"),
        "branch": ("git", "branch", "--show-current"),
        "status_short": ("git", "status", "--short"),
    }
    for key, command in commands.items():
        try:
            completed = subprocess.run(command, capture_output=True, text=True)
        except OSError as exc:
            state[key] = f"unavailable:{exc}"
            continue
        if completed.returncode != 0:
            state[key] = {
                "returncode": completed.returncode,
                "stderr": completed.stderr.strip(),
            }
        else:
            state[key] = completed.stdout.strip()
    return state


__all__ = ["current_git_state"]
