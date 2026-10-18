"""AWS shard merges refuse shards whose work did not run or used an old schema."""

import json

import pandas as pd
import pytest

from benchmarks.cloud import aws_traversal_sibling_fdr_null as fdr_null
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


def test_traversal_fdr_merge_rejects_shard_from_old_schema(tmp_path):
    configured = fdr_null.AwsTraversalSiblingFdrConfig(
        layer=fdr_null.FdrLayer.SYNTHETIC_VALID_P,
        case_names=("synthetic",),
        output_dir=tmp_path,
        replicates=2,
        alpha=0.05,
        base_seed=0,
        shard_count=2,
    )
    schemas = (fdr_null.SCHEMA_VERSION, "traversal_sibling_fdr_null/v1")
    for index, schema in enumerate(schemas):
        shard_dir = tmp_path / "shards" / f"shard_{index:04d}"
        shard_dir.mkdir(parents=True)
        (shard_dir / fdr_null.SHARD_MANIFEST_NAME).write_text(
            json.dumps(
                {"replicate_indices": [index], "source_manifest": {"schema_version": schema}}
            )
        )

    with pytest.raises(ValueError, match="shard_0001 must use source schema"):
        fdr_null.merge_shards(configured)

