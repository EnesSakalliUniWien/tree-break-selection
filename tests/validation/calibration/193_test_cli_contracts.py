from __future__ import annotations

import sys
from pathlib import Path

import pytest
from benchmarks.diagnostics.calibration.cli import parse_null_calibration_panel_args
from benchmarks.diagnostics.calibration.selected.hierarchy.cli import (
    parse_selected_hierarchy_args,
)
from benchmarks.diagnostics.calibration.traversal.cli import (
    parse_traversal_audit_args,
)


def test_traversal_audit_cli_contract(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output_dir = tmp_path / "out"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "traversal-audit",
            "--output-dir",
            str(output_dir),
            "--suite",
            "binary",
            "--case-names",
            "case_a,case_b",
            "--methods",
            "tbs,tbs_nnls",
            "--significance-level",
            "0.02",
            "--edge-alpha",
            "0.003",
        ],
    )

    args = parse_traversal_audit_args(
        description="test audit",
        default_case_names=("default_case",),
        default_methods=("tbs",),
        default_significance_level=0.01,
        default_edge_alpha=0.001,
    )

    assert args.output_dir == output_dir
    assert args.suite == "binary"
    assert args.case_names == ("case_a", "case_b")
    assert args.methods == ("tbs", "tbs_nnls")
    assert args.significance_level == pytest.approx(0.02)
    assert args.edge_alpha == pytest.approx(0.003)


def test_null_calibration_panel_cli_contract(tmp_path: Path) -> None:
    panel_path = tmp_path / "panel.csv"
    output_dir = tmp_path / "out"

    args = parse_null_calibration_panel_args(
        [
            "--panel",
            str(panel_path),
            "--output-dir",
            str(output_dir),
            "--alpha",
            "0.02",
            "--tolerance",
            "0.03",
            "--min-rows",
            "40",
        ],
        description="test null panel",
    )

    assert args.panel == panel_path
    assert args.output_dir == output_dir
    assert args.alpha == pytest.approx(0.02)
    assert args.tolerance == pytest.approx(0.03)
    assert args.min_rows == 40


def test_selected_hierarchy_cli_contract(tmp_path: Path) -> None:
    output_dir = tmp_path / "out"

    args = parse_selected_hierarchy_args(
        [
            "--case-names",
            "case_a,case_b",
            "--n-replicates",
            "12",
            "--seed",
            "37",
            "--target-mode",
            "non_root_strongest",
            "--context-match",
            "projection_parent_size_depth",
            "--output-dir",
            str(output_dir),
        ],
        description="test hierarchy panel",
        default_case_names=("default_case",),
        default_seed=1,
    )

    assert args.case_names == "case_a,case_b"
    assert args.n_replicates == 12
    assert args.seed == 37
    assert args.target_mode == "non_root_strongest"
    assert args.context_match == "projection_parent_size_depth"
    assert args.output_dir == output_dir
