from __future__ import annotations

import io
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from benchmarks.cloud.sharding import iter_expected_shard_dirs, validate_replicate_coverage
from benchmarks.shared.util.values import bool_series, bool_sum, bool_value, parse_csv_list


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (math.nan, False),
        (None, False),
        ("false", False),
        ("no", False),
        ("yes", True),
        (" TRUE ", True),
        (1.0, True),
        (0, False),
        (np.bool_(True), True),
    ],
)
def test_bool_value_treats_missing_and_false_text_as_false(value, expected) -> None:
    assert bool_value(value) is expected


def test_bool_series_reads_a_01_column_with_a_blank_cell() -> None:
    column = pd.read_csv(io.StringIO("id,x\na,1\nb,0\nc,\nd,1\n"))["x"]

    assert column.dtype == float
    assert bool_series(column).tolist() == [True, False, False, True]


def test_bool_sum_counts_true_cells_and_a_missing_column_as_zero() -> None:
    frame = pd.DataFrame({"x": [1, 0, None, 1]})

    assert bool_sum(frame, "x") == 2
    assert bool_sum(frame, "absent") == 0


def test_parse_csv_list_strips_and_drops_empty_items() -> None:
    assert parse_csv_list(" a, b,,c ,") == ("a", "b", "c")
    assert parse_csv_list(None) == ()


def test_validate_replicate_coverage_rejects_duplicate_and_missing_indices() -> None:
    validate_replicate_coverage(
        [{"replicate_indices": [0, 2]}, {"replicate_indices": [1]}],
        expected_replicates=3,
    )
    with pytest.raises(ValueError, match="Duplicate"):
        validate_replicate_coverage(
            [{"replicate_indices": [0, 1]}, {"replicate_indices": [1, 2]}],
            expected_replicates=3,
        )
    with pytest.raises(ValueError, match="Missing"):
        validate_replicate_coverage([{"replicate_indices": [0, 2]}], expected_replicates=3)


def test_iter_expected_shard_dirs_uses_zero_padded_shard_folders() -> None:
    assert iter_expected_shard_dirs(Path("out"), 2) == (
        Path("out/shards/shard_0000"),
        Path("out/shards/shard_0001"),
    )
    with pytest.raises(ValueError, match="positive"):
        iter_expected_shard_dirs(Path("out"), 0)
