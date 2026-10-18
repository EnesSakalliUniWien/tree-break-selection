"""DataFrame adapters for typed benchmark result rows."""

from __future__ import annotations

from dataclasses import fields
from typing import Iterable

import pandas as pd

from benchmarks.shared.result_records.models import BenchmarkResultRow

_ROW_FIELDS = tuple(field for field in fields(BenchmarkResultRow) if field.name != "params_raw")


def _column_name(field_name: str) -> str:
    return "params" if field_name == "params_display" else field_name


RESULT_COLUMNS = [_column_name(field.name) for field in _ROW_FIELDS]
NUMERIC_RESULT_COLUMNS = tuple(
    _column_name(field.name)
    for field in _ROW_FIELDS
    if field.type not in {"str", "BenchmarkRunStatus"}
)


def benchmark_rows_to_dataframe(rows: Iterable[BenchmarkResultRow]) -> pd.DataFrame:
    """Convert typed benchmark rows to the stable output DataFrame schema."""
    records = [
        {
            _column_name(field.name): (
                row.status.value if field.name == "status" else getattr(row, field.name)
            )
            for field in _ROW_FIELDS
        }
        for row in rows
    ]
    return pd.DataFrame(records, columns=RESULT_COLUMNS)


__all__ = [
    "NUMERIC_RESULT_COLUMNS",
    "RESULT_COLUMNS",
    "benchmark_rows_to_dataframe",
]
