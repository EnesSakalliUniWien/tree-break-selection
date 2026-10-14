import json

import numpy as np
from benchmarks.diagnostics.calibration.sibling.known_branch_skew_diagnostic import (
    KnownBranchConfig,
    _fixed_tree,
    _traversal_row,
    expected_node_contrasts,
    population_means,
    run_diagnostic,
    summarize_tables,
    write_diagnostic,
)


def test_named_patterns_preserve_population_contrasts_and_null_membership() -> None:
    expected = {
        "global_null": {"R": True, "A": True, "B": True},
        "root_only": {"R": False, "A": True, "B": True},
        "A_deep_only": {"R": True, "A": False, "B": True},
        "separate_deep_branches": {"R": True, "A": False, "B": False},
        "strong_root_weak_A": {"R": False, "A": False, "B": True},
    }
    for pattern, node_nulls in expected.items():
        means = population_means(pattern, 4)
        contrasts = expected_node_contrasts(means)
        assert {node: is_null for node, (_contrast, is_null) in contrasts.items()} == node_nulls

    weak_a = expected_node_contrasts(population_means("strong_root_weak_A", 4))["A"][0]
    assert np.allclose(weak_a, [0.0, 0.5, 0.0, 0.0])


def test_fixed_tree_traversal_counts_only_reached_splits() -> None:
    truths = expected_node_contrasts(population_means("global_null", 2))
    tests = {
        "R": {"statistic": 0.2, "degrees_of_freedom": 2, "p_value": 0.8},
        "A": {"statistic": 15.0, "degrees_of_freedom": 2, "p_value": 0.001},
        "B": {"statistic": 0.3, "degrees_of_freedom": 2, "p_value": 0.7},
    }
    no_passthrough = _traversal_row(
        _fixed_tree(), tests, truths, alpha=0.01, passthrough=False
    )
    with_passthrough = _traversal_row(
        _fixed_tree(), tests, truths, alpha=0.01, passthrough=True
    )

    assert no_passthrough["n_reached_splits"] == 0
    assert no_passthrough["n_reached_false_splits"] == 0
    assert with_passthrough["n_reached_splits"] == 1
    assert with_passthrough["false_split_node_ids"] == "A"
    injected = json.loads(with_passthrough["injected_local_tests"])
    assert injected["A"] == {"degrees_of_freedom": 2.0, "p_value": 0.001, "statistic": 15.0}
    decisions = json.loads(with_passthrough["traversal_decisions"])
    assert decisions["R"] == "pass_through"
    assert decisions["A"] == "split"


def test_small_replicable_run_has_paired_stages_and_summary_intervals(tmp_path) -> None:
    config = KnownBranchConfig(
        output_dir=tmp_path / "diagnostic",
        replicates=3,
        base_seed=17,
        n_per_group=12,
        dimension=4,
        patterns=("global_null", "strong_root_weak_A"),
        covariance_profiles=("identity",),
    )
    tables = run_diagnostic(config)
    summary = summarize_tables(tables, alpha=config.alpha)

    assert len(tables["local"]) == 2 * 3 * 3 * 4
    assert len(tables["selected_root"]) == 2 * 3
    assert len(tables["traversal"]) == 2 * 3 * 2
    assert tables["local"]["dataset_sha256"].nunique() == 6
    assert tables["local"].groupby(["pattern", "replicate"])["dataset_sha256"].nunique().eq(1).all()
    assert tables["local"]["population_contrast"].str.len().gt(2).all()
    local_summaries = summary[summary["table"] == "local"]
    assert local_summaries["null_n"].ge(0).all()
    assert local_summaries["null_rejection_ci95_low"].notna().any()
    assert tables["traversal"]["injected_local_tests"].str.len().gt(2).all()

    output_paths = write_diagnostic(config, tables, run_command="test diagnostic")
    manifest = json.loads(output_paths["manifest"].read_text(encoding="utf-8"))
    assert manifest["config"]["base_seed"] == 17
    assert manifest["runtime"]["python"]
    assert manifest["source_sha256"]
    assert output_paths["report"].is_file()


def test_writer_refuses_to_overwrite_nonempty_output_directory(tmp_path) -> None:
    config = KnownBranchConfig(output_dir=tmp_path, replicates=1)
    (tmp_path / "keep.txt").write_text("existing", encoding="utf-8")
    try:
        write_diagnostic(config, {"local": None, "selected_root": None, "traversal": None}, run_command="test")
    except FileExistsError:
        pass
    else:
        raise AssertionError("expected populated-output guard")
