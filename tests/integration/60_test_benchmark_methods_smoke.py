"""Smoke test for benchmark method runners."""

import pytest
from benchmarks.shared.cases import SMALL_TEST_CASES, get_default_test_cases
from benchmarks.shared.pipeline import benchmark_cluster_algorithm


def test_benchmark_graph_and_density_methods_smoke():
    """Run one graph method and one density method on a single test case."""
    case = SMALL_TEST_CASES[0].copy()
    df_results, _ = benchmark_cluster_algorithm(
        test_cases=[case],
        verbose=False,
        plot_umap=False,
        methods=["leiden", "dbscan"],
    )

    assert len(df_results) == 2
    assert set(df_results["method"]) == {"leiden", "dbscan"}

    dbscan_row = df_results[df_results["method"] == "dbscan"].iloc[0]
    assert dbscan_row["status"] == "ok"
    assert dbscan_row["labels_length"] == dbscan_row["samples"]

    ok_rows = df_results[df_results["status"] == "ok"]
    assert (ok_rows["labels_length"] == ok_rows["samples"]).all()


def test_benchmark_louvain_and_adaptive_diffusion_methods_smoke():
    """Separate a working comparator from fail-closed selected-hierarchy TBS."""
    case = SMALL_TEST_CASES[0].copy()
    df_results, _ = benchmark_cluster_algorithm(
        test_cases=[case],
        verbose=False,
        plot_umap=False,
        methods=["louvain", "tbs_diffusion_adaptive"],
    )

    assert len(df_results) == 2
    assert set(df_results["method"]) == {"louvain", "tbs_diffusion_adaptive"}
    louvain_row = df_results[df_results["method"] == "louvain"].iloc[0]
    assert louvain_row["status"] == "ok"
    assert louvain_row["labels_length"] == louvain_row["samples"]
    tbs_row = df_results[df_results["method"] == "tbs_diffusion_adaptive"].iloc[0]
    assert tbs_row["status"] == "unsupported"
    assert tbs_row["labels_length"] == 0
    assert tbs_row["unsupported_reason_code"] == (
        "empirical_null_unvalidated_reference_law"
    )
    assert tbs_row["unsupported_admissible_support_count"] > 0


@pytest.mark.optional
@pytest.mark.graphtools
def test_benchmark_graphtools_diffusion_method_smoke(require_optional_dependencies):
    """Run the optional graphtools diffusion backend when it is installed."""
    require_optional_dependencies("graphtools")
    case = SMALL_TEST_CASES[0].copy()
    with pytest.warns(RuntimeWarning, match="Detected zero distance between samples"):
        df_results, _ = benchmark_cluster_algorithm(
            test_cases=[case],
            verbose=False,
            plot_umap=False,
            methods=[
                "tbs_diffusion_graphtools",
                "tbs_diffusion_graphtools_adaptive_nnls",
            ],
        )

    assert set(df_results["method"]) == {
        "tbs_diffusion_graphtools",
        "tbs_diffusion_graphtools_adaptive_nnls",
    }
    assert len(df_results) == 9
    assert df_results["run_id"].nunique() == 9
    assert (
        df_results[df_results["method"] == "tbs_diffusion_graphtools_adaptive_nnls"][
            "benchmark_grid"
        ]
        == "graphtools_adaptive_k_tree_strategy"
    ).all()
    assert set(df_results["status"]) == {"unsupported"}
    assert (df_results["labels_length"] == 0).all()
    assert df_results["unsupported_reason_code"].eq(
        "empirical_null_unvalidated_reference_law"
    ).all()
    assert (df_results["unsupported_admissible_support_count"] > 0).all()


def test_hamming_diffusion_rejects_continuous_benchmark_input():
    """The Hamming diffusion method must explicitly skip continuous cases."""
    case = next(
        case.copy()
        for case in get_default_test_cases()
        if case["name"] == "gauss_clear_medium_continuous"
    )
    df_results, _ = benchmark_cluster_algorithm(
        test_cases=[case],
        verbose=False,
        plot_umap=False,
        methods=["tbs_diffusion"],
    )

    assert len(df_results) == 1
    row = df_results.iloc[0]
    assert row["status"] == "skip"
    assert row["labels_length"] == 0
    assert "requires binary or one-hot" in row["skip_reason"]
