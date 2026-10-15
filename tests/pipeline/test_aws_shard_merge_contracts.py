"""AWS shard merges refuse shards whose work did not run."""

import json

import pandas as pd
import pytest

from benchmarks.cloud.aws_selected_tail_equation_study import (
    SHARD_CASE_SUMMARY_NAME,
    SHARD_MANIFEST_NAME,
    SHARD_RECORDS_NAME,
    SIMULATION_ID_COLUMN,
    load_shard_records,
)


def _write_shard(shard_dir, statuses):
    shard_dir.mkdir()
    pd.DataFrame({SIMULATION_ID_COLUMN: ["sim0"]}).to_csv(
        shard_dir / SHARD_RECORDS_NAME, index=False
    )
    (shard_dir / SHARD_MANIFEST_NAME).write_text(json.dumps({"shard_index": 0, "shard_seed": 7}))
    pd.DataFrame(
        {
            "case_id": [f"case_{i}" for i in range(len(statuses))],
            "status": statuses,
            "skip_reason": ["" if status == "ok" else "boom" for status in statuses],
        }
    ).to_csv(shard_dir / SHARD_CASE_SUMMARY_NAME, index=False)


def test_tail_study_shard_with_skipped_case_is_rejected(tmp_path):
    _write_shard(tmp_path / "shard_0000", ["ok", "skip"])

    with pytest.raises(ValueError, match="did not run.*case_1.*boom"):
        load_shard_records(tmp_path / "shard_0000")


def test_tail_study_shard_with_all_cases_ok_loads(tmp_path):
    _write_shard(tmp_path / "shard_0000", ["ok"])

    records = load_shard_records(tmp_path / "shard_0000")

    assert records[SIMULATION_ID_COLUMN].tolist() == ["shard_0000:sim0"]
