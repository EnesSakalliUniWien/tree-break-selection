"""Scalar coercion shared by every calibration diagnostic category.

These helpers read one field out of a diagnostic row and return a total value,
so panels can build result tables without repeating the same guards. They are
deliberately permissive: a missing column, an unparseable value, or a pandas
missing marker yields the neutral result rather than raising, because a
diagnostic row is evidence and a single bad field must not abort a panel.

`root/root_values.py` pairs `finite_float` with selected-root input-column
validation.
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


def bool_value(value: object) -> bool:
    """Read one table cell as a boolean: missing is false; text accepts 1/true/t/yes/y."""
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "t", "yes", "y"}
    if pd.isna(value):
        return False
    return bool(value)


def bool_series(series: pd.Series) -> pd.Series:
    """Read a table column cell by cell with the bool_value rules."""
    return series.map(bool_value).astype(bool)


def parse_csv_list(raw: str) -> tuple[str, ...]:
    """Split a comma-separated CLI value into stripped, non-empty items."""
    return tuple(part.strip() for part in str(raw).split(",") if part.strip())


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
    "bool_series",
    "bool_value",
    "finite_float",
    "parse_csv_list",
    "string_value",
]
