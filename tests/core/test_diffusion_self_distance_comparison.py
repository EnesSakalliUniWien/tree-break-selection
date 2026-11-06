"""Topology comparisons separate merge ordering from rooted clades."""

import numpy as np
import pytest

from benchmarks.diagnostics.calibration.diffusion_self_distance_comparison import (
    compare_records,
    main,
    observe_dataset,
)


def test_merge_order_change_is_not_mislabeled_as_clade_or_root_change():
    common = dict(
        case="control",
        replicate=0,
        arm="original",
        epsilon=1.0,
        known_identity_statistic=1.0,
        unselected_chi_p=0.5,
        positive_self_count=0,
        clades=[[0, 1], [2, 3], [0, 1, 2, 3]],
        root_groups=[["a", "b"], ["c", "d"]],
    )
    before = common | dict(signature=[[0, 1], [2, 3], [4, 5]])
    after = common | dict(signature=[[2, 3], [0, 1], [4, 5]])
    row = compare_records([before], [after]).iloc[0]
    assert row.full_history_changed
    assert not row.clades_changed
    assert not row.root_partition_changed
    assert row.clade_symmetric_difference == 0


def test_native_comparison_observes_both_arms_without_mutating_data():
    values = np.random.default_rng(20261105).normal(size=(9, 1))
    before = values.copy()
    rows = observe_dataset(values, [f"S{i}" for i in range(9)], "small", 0, 20261105)
    assert [row["arm"] for row in rows] == ["original", "reconstructed"]
    assert all(row["n"] == 9 and row["p"] == 1 for row in rows)
    assert all(len(row["signature"]) == 8 for row in rows)
    assert all(row["n_left"] + row["n_right"] == 9 for row in rows)
    assert all(0 <= row["unselected_chi_p"] <= 1 for row in rows)
    assert all(row["positive_self_count"] == 0 for row in rows)
    assert all(row["native_vs_zero_self_reference_operator_max_abs"] < 1e-14 for row in rows)
    np.testing.assert_array_equal(values, before)


@pytest.mark.parametrize(
    "stage,filename", [("baseline", "baseline_runs.json"), ("compare", "post_fix_runs.json")]
)
def test_cli_preserves_existing_raw_evidence(tmp_path, monkeypatch, stage, filename):
    target = tmp_path / filename
    target.write_text("preserve this captured evidence")
    monkeypatch.setattr("sys.argv", ["comparison", "--stage", stage, "--output", str(tmp_path)])
    with pytest.raises(SystemExit) as error:
        main()
    assert error.value.code == 2
    assert target.read_text() == "preserve this captured evidence"
