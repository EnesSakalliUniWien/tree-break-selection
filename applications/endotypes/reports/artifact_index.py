"""Shared artifact-index loading for endotype report generators."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def load_artifact_index(experiment_dir: Path) -> pd.DataFrame:
    """Load and normalize an endotype experiment artifact index."""
    path = experiment_dir / "artifact_index.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    frame = pd.read_csv(path)
    frame["_rank"] = pd.to_numeric(frame["specificity_aware_rank"], errors="coerce")
    return frame.sort_values(["_rank", "weighting", "block_name"], na_position="last").reset_index(
        drop=True
    )
