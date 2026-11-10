"""Canonical benchmark execution statuses."""

from __future__ import annotations

from enum import Enum


class BenchmarkRunStatus(str, Enum):
    """Mutually exclusive scientific outcomes for one benchmark method run."""

    OK = "ok"
    SKIP = "skip"
    UNSUPPORTED = "unsupported"
    # The method raised; skip_reason holds the exception so other methods' rows survive.
    ERROR = "error"

    def __str__(self) -> str:
        """Return the stable wire value used in CSV and DataFrame outputs."""
        return self.value


__all__ = ["BenchmarkRunStatus"]
