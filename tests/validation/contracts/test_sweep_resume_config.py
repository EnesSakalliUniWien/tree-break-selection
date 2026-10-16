"""Resumed sweeps reuse saved results only for the configuration that wrote them."""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import pandas as pd
import pytest

from benchmarks.diagnostics.calibration.statistics import mixed_internal_calibration_sweep
from benchmarks.validation.sweeps import alpha_grid_search, alpha_structure_sweep


def test_mixed_internal_resume_rejects_another_configuration(tmp_path: Path) -> None:
    (tmp_path / "manifest.json").write_text(
        json.dumps(
            {
                "schema_version": mixed_internal_calibration_sweep.SCHEMA_VERSION,
                "suite": "binary",
                "case_names": ["binary_low_noise_2c"],
                "n_replicates": 1,
                "edge_alpha": mixed_internal_calibration_sweep.DEFAULT_EDGE_ALPHA,
                "sibling_alpha": mixed_internal_calibration_sweep.DEFAULT_SIBLING_ALPHA,
            }
        )
    )

    with pytest.raises(ValueError, match="holds results for"):
        mixed_internal_calibration_sweep.run_mixed_internal_calibration_sweep(
            suite="binary",
            case_names=("binary_low_noise_2c",),
            n_replicates=2,
            output_dir=tmp_path,
            resume=True,
        )


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
        return result, None

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


def test_alpha_structure_resume_rejects_other_cases(tmp_path: Path) -> None:
    pair_dir = tmp_path / alpha_grid_search.alpha_pair_id(0.001, 0.01)
    pair_dir.mkdir()
    pd.DataFrame({"case_id": ["binary_perfect_2c"], "method": ["tbs"]}).to_csv(
        pair_dir / "benchmark_metrics.csv", index=False
    )
    pd.DataFrame({"alpha_pair_id": ["stale"]}).to_csv(
        pair_dir / "structural_summary.csv", index=False
    )
    config = alpha_structure_sweep.AlphaStructureSweepConfig(
        suite="binary",
        case_names=("binary_low_noise_2c",),
        methods=("tbs",),
        edge_alphas=(0.001,),
        sibling_alphas=(0.01,),
        output_dir=tmp_path,
        resume=True,
        baseline_edge_alpha=0.001,
        baseline_sibling_alpha=0.01,
    )

    with pytest.raises(ValueError, match="other cases or methods"):
        alpha_structure_sweep.run_alpha_structure_sweep(config)
