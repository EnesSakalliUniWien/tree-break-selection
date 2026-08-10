"""Shared value parsing and transforms for selected-root diagnostics."""

from __future__ import annotations

import math

import pandas as pd

from benchmarks.diagnostics.calibration.values import finite_float


def require_columns(frame: pd.DataFrame, columns: set[str], label: str) -> None:
    """Validate that a diagnostic input table has the required columns."""
    missing = columns - set(frame.columns)
    if missing:
        raise ValueError(f"{label} missing required columns: {sorted(missing)!r}.")


def safe_log1p(value: object) -> float:
    """Return log1p of the nonnegative part of a numeric diagnostic value."""
    numeric = finite_float(value)
    if not math.isfinite(numeric):
        return math.nan
    return float(math.log1p(max(numeric, 0.0)))


def spectral_excess_log(value: object) -> float:
    """Return the nonnegative log spectral excess used by root diagnostics."""
    numeric = finite_float(value)
    if not math.isfinite(numeric):
        return math.nan
    return float(max(math.log(max(numeric, 1e-12)), 0.0))


__all__ = [
    "finite_float",
    "require_columns",
    "safe_log1p",
    "spectral_excess_log",
]
