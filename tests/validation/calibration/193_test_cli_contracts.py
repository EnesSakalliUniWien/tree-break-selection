from __future__ import annotations

from pathlib import Path

from benchmarks.diagnostics.calibration.selected.hierarchy.cli import (
    parse_selected_hierarchy_args,
)


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
