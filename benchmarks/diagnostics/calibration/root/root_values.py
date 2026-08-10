"""Shared value parsing and column validation for selected-root diagnostics."""

from __future__ import annotations

import pandas as pd

from benchmarks.diagnostics.calibration.values import finite_float


def require_columns(frame: pd.DataFrame, columns: set[str], label: str) -> None:
    """Validate that a diagnostic input table has the required columns."""
    missing = columns - set(frame.columns)
    if missing:
        raise ValueError(f"{label} missing required columns: {sorted(missing)!r}.")


__all__ = [
    "finite_float",
    "require_columns",
]
