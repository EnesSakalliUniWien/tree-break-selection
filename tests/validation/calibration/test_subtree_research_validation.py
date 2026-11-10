import json

import numpy as np
import pandas as pd
import pytest

from benchmarks.validation.sweeps.subtree_research_validation import (
    calibration_summary,
    paired_recovery_summary,
)


def _calibration():
    return pd.DataFrame(
        [
            (0, 0, "ok", 0.0, "reject"),
            (0, 1, "ok", 0.02, "do_not_reject"),
            (0, 2, "unresolved", None, "unresolved"),
            (0, 3, "error", None, None),
        ],
        columns=["cell", "replicate", "status", "p_value", "decision"],
    )


def _recovery(cells=2, replicates=4):
    rows = []
    for cell in range(cells):
        for rep in range(replicates):
            seed = replicates * cell + rep
            for method, ami in (("base", 0.30), ("test", 0.45 + 0.1 * (cell % 2))):
                rows.append(
                    {
                        "cell": cell,
                        "replicate": rep,
                        "method": method,
                        "status": "ok",
                        "data_hash": f"{1 + seed:064x}",
                        "seed": seed,
                        "ami": ami,
                        "matched_macro_f1": ami,
                    }
                )
    return pd.DataFrame(rows)


def test_calibration_denominators_holm_and_unresolved_rows():
    outcome = calibration_summary(_calibration(), expected_cells=2, expected_replicates=4)
    assert not outcome["complete"] and not outcome["passed"]
    assert outcome["counts"]["unresolved"] == outcome["counts"]["error"] == 1
    first, second, empty, _ = outcome["comparisons"]
    assert (first["attempts"], first["rejections"], first["rate"]) == (4, 1, 0.25)
    assert (second["attempts"], second["rejections"], second["rate"]) == (4, 2, 0.5)
    assert first["excess_p_holm"] >= first["excess_p"]
    assert empty["attempts"] == 0 and empty["interval"] is None
    json.dumps(outcome, allow_nan=False)
    unresolved = _calibration().assign(status="unresolved", p_value=None, decision="unresolved")
    result = calibration_summary(unresolved, expected_cells=1, expected_replicates=4)
    assert result["complete"] and not result["all_successful"] and not result["passed"]


@pytest.mark.parametrize(
    "transform, reason",
    [
        (lambda f: pd.concat([f, f.iloc[:1]]), "duplicate"),
        (lambda f: f.assign(cell=9), "cell"),
        (lambda f: f.assign(p_value="broken"), "p_value"),
        (lambda f: f.assign(p_value=1.1), "p_value"),
        (lambda f: f.assign(status="unknown"), "status"),
        (lambda f: f.assign(decision="typo"), "decision"),
    ],
)
def test_calibration_rejects_invalid_rows(transform, reason):
    with pytest.raises(ValueError, match=reason):
        calibration_summary(transform(_calibration()), expected_cells=1, expected_replicates=4)


def test_recovery_equal_cell_blocks_and_support_denominator():
    frame = _recovery()
    result = paired_recovery_summary(frame, "base", "test", expected_cells=2, expected_replicates=4)
    assert result["complete"] and result["all_successful"] and not result["passed"]
    assert result["aggregate"]["blocks"] == 4
    assert result["aggregate"]["mean_ami_difference"] == pytest.approx(0.2)
    assert result["aggregate"]["interval"] == pytest.approx([0.2, 0.2])
    assert result["cells"][0]["auxiliary"]["matched_macro_f1"]["mean_difference"] == pytest.approx(
        0.15
    )
    frame.loc[0, ["status", "ami"]] = ["error", np.nan]
    result = paired_recovery_summary(
        frame.iloc[:-1], "base", "test", expected_cells=2, expected_replicates=4
    )
    assert not result["complete"] and not result["all_successful"] and not result["passed"]
    assert result["coverage"]["paired_successful"] == 6
    assert result["coverage"]["paired_share"] == 0.75
    assert result["aggregate"]["blocks"] == 2
    assert result["cells"][0]["baseline_status_counts"]["error"] == 1
    assert result["cells"][1]["candidate_status_counts"]["missing"] == 1
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize(
    "column,value,reason",
    [
        ("seed", 900, "seed"),
        ("data_hash", "a" * 64, "data_hash"),
        ("method", "other", "method"),
        ("ami", 2.0, "ami"),
    ],
)
def test_paired_recovery_rejects_mismatched_or_malformed_rows(column, value, reason):
    frame = _recovery()
    frame.loc[1, column] = value
    with pytest.raises(ValueError, match=reason):
        paired_recovery_summary(frame, "base", "test", expected_cells=2, expected_replicates=4)


def test_constant_paired_losses_are_flagged_even_if_overall_ami_improves():
    frame = _recovery()
    frame.loc[(frame.cell == 0) & (frame.method == "test"), "ami"] = 0.25
    result = paired_recovery_summary(frame, "base", "test", expected_cells=2, expected_replicates=4)
    assert result["aggregate"]["interval"][0] > 0
    assert result["cells"][0]["loss_p_holm"] == 0
    assert result["cells"][1]["loss_p_holm"] == 1


def test_all_one_calibration_is_degenerate():
    frame = _calibration().assign(status="ok", p_value=1.0, decision="do_not_reject")
    result = calibration_summary(frame, expected_cells=1, expected_replicates=4)
    assert result["degenerate_cells"] == [0]
    assert not result["excess_detected"] and not result["passed"]
    assert result["comparisons"][0]["one_p_values"] == 4


def test_unresolved_decision_cannot_count_as_completed_even_with_ok_status():
    frame = _calibration().assign(status="ok", p_value=0, decision="unresolved")
    result = calibration_summary(frame, expected_cells=1, expected_replicates=4)
    assert result["complete"] and not result["all_successful"]
    assert result["counts"]["ok"] == 0 and result["counts"]["unresolved"] == 4
    assert all(row["rejections"] == row["zero_p_values"] == 0 for row in result["comparisons"])


def test_recovery_rejects_reuse_of_dataset_seed_or_hash():
    frame = _recovery()
    frame.loc[2:3, "seed"] = frame.loc[0, "seed"]
    with pytest.raises(ValueError, match="seed.*reused"):
        paired_recovery_summary(frame, "base", "test", expected_cells=2, expected_replicates=4)


def test_frozen_calibration_gate_requires_complete_nondegenerate_attempts():
    p_values = (np.arange(5_000) + 0.5) / 5_000
    frame = pd.DataFrame(
        {
            "cell": np.repeat(np.arange(24), 5_000),
            "replicate": np.tile(np.arange(5_000), 24),
            "status": "ok",
            "p_value": np.tile(p_values, 24),
            "decision": "do_not_reject",
            "phase": 2,
        }
    )
    frame.loc[frame.p_value <= 0.05, "decision"] = "reject"
    summary = calibration_summary(frame)
    assert summary["passed"] and len(summary["comparisons"]) == 48
    assert all(row["excess_p_holm"] is not None for row in summary["comparisons"])
    assert not calibration_summary(frame.assign(phase=1))["passed"]
    unverified = calibration_summary(frame.drop(columns="phase"))
    assert not unverified["passed"] and unverified["design"] == "development_or_custom"
    frame.loc[frame.cell == 0, ["p_value", "decision"]] = [1.0, "do_not_reject"]
    summary = calibration_summary(frame)
    assert (
        summary["all_successful"] and summary["degenerate_cells"] == [0] and not summary["passed"]
    )


def test_frozen_recovery_gate_requires_explicit_phase_two():
    frame = _recovery(cells=16, replicates=200)
    assert not paired_recovery_summary(frame, "base", "test")["passed"]
    assert not paired_recovery_summary(frame.assign(phase=1), "base", "test")["passed"]
    summary = paired_recovery_summary(frame.assign(phase=2), "base", "test")
    assert summary["passed"] and summary["aggregate"]["blocks"] == 200
