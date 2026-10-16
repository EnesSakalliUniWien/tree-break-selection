"""Smoke test for benchmark method runners."""

import pytest

from benchmarks.shared.cases import get_default_test_cases
from benchmarks.shared.pipeline import benchmark_cluster_algorithm

SMOKE_CASE = {
    "name": "clear",
    "generator": "blobs",
    "n_samples": 24,
    "n_features": 12,
    "n_clusters": 3,
    "cluster_std": 0.4,
    "seed": 0,
    "category": "small_smoke",
}


def test_benchmark_graph_and_density_methods_smoke():
    """Run one graph method and one density method on a single test case."""
    case = SMOKE_CASE.copy()
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
    """Run the previously nonfunctional methods through one benchmark case."""
    case = SMOKE_CASE.copy()
    df_results, _ = benchmark_cluster_algorithm(
        test_cases=[case],
        verbose=False,
        plot_umap=False,
        methods=["louvain", "tbs_diffusion_adaptive"],
    )

    assert len(df_results) == 2
    assert set(df_results["method"]) == {"louvain", "tbs_diffusion_adaptive"}
    assert set(df_results["status"]) == {"ok"}
    assert (df_results["labels_length"] == df_results["samples"]).all()


@pytest.mark.optional
@pytest.mark.graphtools
def test_benchmark_graphtools_diffusion_method_smoke(require_optional_dependencies):
    """Run the optional graphtools diffusion backend when it is installed."""
    require_optional_dependencies("graphtools")
    case = SMOKE_CASE.copy()
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
    # 1 graphtools row plus the 5-linkage adaptive-K tree-strategy grid.
    assert len(df_results) == 6
    assert df_results["run_id"].nunique() == 6
    assert (
        df_results[df_results["method"] == "tbs_diffusion_graphtools_adaptive_nnls"][
            "benchmark_grid"
        ]
        == "graphtools_adaptive_k_tree_strategy"
    ).all()
    assert set(df_results["status"]) == {"ok"}
    assert (df_results["labels_length"] == df_results["samples"]).all()


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
