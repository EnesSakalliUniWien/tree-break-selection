"""Scalar coercion shared by every calibration diagnostic category.

These helpers read one field out of a diagnostic row and return a total value,
so panels can build result tables without repeating the same guards. They are
deliberately permissive: a missing column, an unparseable value, or a pandas
missing marker yields the neutral result rather than raising, because a
diagnostic row is evidence and a single bad field must not abort a panel.

`root/root_values.py` builds on `finite_float` for selected-root diagnostics
that also need root-specific transforms.
"""

from __future__ import annotations

import math

import pandas as pd


def finite_float(value: object) -> float:
    """Return a finite float or ``nan`` for invalid diagnostic values."""
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return math.nan
    return numeric if math.isfinite(numeric) else math.nan


def string_value(
    row: pd.Series | dict[str, object],
    column: str,
    default: str = "",
) -> str:
    """Return a non-missing string field from a diagnostic row."""
    if column not in row:
        return default
    value = row[column]
    if pd.isna(value):
        return default
    return str(value)


__all__ = [
    "finite_float",
    "string_value",
]
