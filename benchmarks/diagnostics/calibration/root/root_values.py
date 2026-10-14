"""Shared value parsing and column validation for selected-root diagnostics."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd

from benchmarks.diagnostics.calibration.values import finite_float


def require_columns(frame: pd.DataFrame, columns: set[str], label: str) -> None:
    """Validate that a diagnostic input table has the required columns."""
    missing = columns - set(frame.columns)
    if missing:
        raise ValueError(f"{label} missing required columns: {sorted(missing)!r}.")


def finite_median(values: pd.Series) -> float:
    """Median of the finite numeric values, or nan when there are none."""
    numeric = pd.to_numeric(values, errors="coerce")
    numeric = numeric[np.isfinite(numeric)]
    return float(numeric.median()) if not numeric.empty else math.nan


__all__ = [
    "finite_float",
    "finite_median",
    "require_columns",
]
