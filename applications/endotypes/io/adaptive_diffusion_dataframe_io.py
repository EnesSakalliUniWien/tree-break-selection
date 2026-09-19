"""DataFrame file input/output for the adaptive diffusion GO-IC workflow.

Analysis functions accept and return DataFrames. This module owns CSV/TSV
serialization; callers decide when to read or persist those frames.
"""

from pathlib import Path

import pandas as pd

from applications.endotypes._shared import load_binary_feature_matrix


def load_binary_matrix(path: Path) -> pd.DataFrame:
    return load_binary_feature_matrix(
        path,
        drop_zero_columns=True,
        require_nonzero_rows=True,
        non_binary_message="contains values outside {0,1}.",
    )


def read_dataframe(path: Path | str) -> pd.DataFrame:
    """Read an exported analysis table using the established CSV schema."""
    return pd.read_csv(path)


def write_dataframe(frame: pd.DataFrame, path: Path | str) -> None:
    """Write an analysis table without a synthetic index column."""
    frame.to_csv(path, index=False)
