"""Resumed sweeps reuse saved results only for the configuration that wrote them."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pandas as pd
import pytest

from benchmarks.validation.sweeps import alpha_grid_search


def test_alpha_grid_resume_reuses_only_the_same_case_set(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[list[str]] = []

    def fake_benchmark(*, test_cases: list[dict], **_kwargs: object):
        calls.append([str(case["name"]) for case in test_cases])
        rows = len(test_cases)
        result = pd.DataFrame(
            {
                "method": ["tbs"] * rows,
                "status": ["ok"] * rows,
                "found_clusters": [2] * rows,
                "true_clusters": [2] * rows,
                "ari": [1.0] * rows,
                "skip_reason": [None] * rows,
            }
        )
        return result

    monkeypatch.setattr(alpha_grid_search, "benchmark_cluster_algorithm", fake_benchmark)
    config = alpha_grid_search.AlphaGridConfig(
        suite="binary",
        case_names=("binary_low_noise_2c",),
        edge_alphas=(0.001,),
        sibling_alphas=(0.01,),
        output_dir=tmp_path,
        resume=True,
    )

    alpha_grid_search.run_alpha_grid_search(config)
    alpha_grid_search.run_alpha_grid_search(config)
    assert calls == [["binary_low_noise_2c"]]

    with pytest.raises(ValueError, match="different case set"):
        alpha_grid_search.run_alpha_grid_search(replace(config, case_names=("binary_perfect_2c",)))
