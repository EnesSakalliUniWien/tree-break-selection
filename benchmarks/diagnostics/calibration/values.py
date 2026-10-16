"""Scalar coercion shared by every calibration diagnostic category.

These helpers read one field out of a diagnostic row and return a total value,
so panels can build result tables without repeating the same guards. They are
deliberately permissive: a missing column, an unparseable value, or a pandas
missing marker yields the neutral result rather than raising, because a
diagnostic row is evidence and a single bad field must not abort a panel.
"""

from __future__ import annotations

import math

import numpy as np
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


def finite_int(value: object, default: int = 0) -> int:
    """Return an int, or ``default`` when the value is not a finite number."""
    number = finite_float(value)
    return int(number) if math.isfinite(number) else int(default)


def is_finite(value: object) -> bool:
    """Return whether the value parses as a finite number."""
    return math.isfinite(finite_float(value))


def finite_min(values: object) -> float:
    """Minimum of the finite numeric values (``inf`` excluded), or ``nan`` when there are none."""
    numeric = pd.to_numeric(pd.Series(values), errors="coerce")
    numeric = numeric[np.isfinite(numeric)]
    return float(numeric.min()) if not numeric.empty else math.nan


def numeric_column(rows: pd.DataFrame, column: str) -> pd.Series:
    """Read a required table column as numbers; unparseable cells become ``nan``."""
    return pd.to_numeric(rows[column], errors="coerce")


def split_path(value: object) -> list[str]:
    """Split a ``;``-joined node path cell; a missing cell is an empty path."""
    if pd.isna(value):
        return []
    return [part for part in str(value).split(";") if part]


def output_data_role(data_role: object) -> str:
    """Return a CSV-stable role label.

    Pandas reads the literal string ``null`` as missing, so outputs use
    ``selected_null`` while ``null`` stays the accepted CLI/input role.
    """
    return "selected_null" if str(data_role) == "null" else str(data_role)


def input_data_role(data_role: object) -> str:
    """Map an output role label back to the CLI/input role."""
    return "null" if str(data_role) == "selected_null" else str(data_role)


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
    "finite_int",
    "finite_min",
    "input_data_role",
    "is_finite",
    "numeric_column",
    "output_data_role",
    "parse_csv_list",
    "split_path",
    "string_value",
]
