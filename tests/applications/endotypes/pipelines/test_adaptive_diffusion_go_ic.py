"""Contracts across the modular adaptive diffusion GO-IC workflow."""

from argparse import Namespace
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from applications.endotypes.analysis.go_ic.adaptive_diffusion_go_ic_ranking import (
    add_specificity_aware_rank,
    rank_go_ic_results,
)
from applications.endotypes.analysis.go_ic.adaptive_diffusion_go_ic_terms import (
    component_feature_loadings,
)
from scipy.spatial.distance import pdist
from tree_break_selection.space_separation import SpectralBlock

from applications.endotypes.pipelines import adaptive_diffusion_go_ic_workflow as workflow


def test_feature_loadings_reconstruct_normalized_features():
    values = np.array([[1.0, 0.0, 1.0], [0.0, 1.0, 1.0], [1.0, 1.0, 0.0]])
    normalized = values / np.linalg.norm(values, axis=1, keepdims=True)
    eigenvalues, eigenvectors = np.linalg.eigh(normalized @ normalized.T)
    loadings = component_feature_loadings(values, eigenvalues, eigenvectors)
    np.testing.assert_allclose(
        (eigenvectors * np.sqrt(eigenvalues)) @ loadings.T, normalized, atol=1e-14
    )


def test_specificity_order_preserves_raw_go_ic_rank_and_failed_rows():
    rows = [
        dict(
            run_id=name,
            status="ok",
            n_clusters=5,
            largest_cluster_fraction=0.3,
            singleton_gene_fraction=0.0,
            coherent_cluster_fraction=1.0,
            weighted_mean_within_tfidf_cosine=0.8,
            go_bic_active=bic,
            go_bic_active_per_gene=bic / 10,
        )
        for name, bic in [("lower_bic", 10.0), ("more_specific", 20.0)]
    ]
    rows.append(dict(run_id="failed", status="failed_gate"))
    coherence = pd.DataFrame(
        [
            dict(
                run_id=name, top_term_prevalence_delta=delta, cluster_size=10, coherent_by_rule=True
            )
            for name, delta in [("lower_bic", 0.1), ("more_specific", 0.8)]
        ]
    )
    result = add_specificity_aware_rank(rank_go_ic_results(rows), coherence).set_index("run_id")
    assert result.loc["more_specific", "display_rank"] == 1
    assert result.loc["more_specific", "raw_go_ic_rank"] == 2
    assert result.loc["more_specific", "old_display_rank"] == 2
    assert pd.isna(result.loc["failed", "display_rank"])
    assert result.loc["failed", "status"] == "failed_gate"


@pytest.mark.parametrize("failure", ["failed_diffusion", "failed_gate"])
def test_workflow_keeps_success_and_failure_artifacts(tmp_path, monkeypatch, failure):
    data = pd.DataFrame(
        np.random.default_rng(42).integers(0, 2, (16, 8)), index=[f"gene_{i}" for i in range(16)]
    )
    data[0] = 1
    input_path = tmp_path / "feature_matrix_test.tsv"
    data.to_csv(input_path, sep="\t")
    args = Namespace(
        input=input_path,
        output_dir=tmp_path / "out",
        dataset_label="test",
        edge_alpha=0.001,
        sibling_alpha=0.01,
        max_rank=6,
        min_segment_length=2,
        max_segments=2,
        diffusion_k_neighbors=5,
        diffusion_time=3,
        diffusion_components=4,
        adaptive_metric="euclidean",
        adaptive_bandwidth_type="-1/(d+2)",
        adaptive_epsilon="median",
        weightings=["binary", "tfidf"],
        block_names=["adaptive_common_mode_01"],
        top_terms_per_axis=2,
    )
    calls = {"geometry": 0, "gate": 0}
    monkeypatch.setattr(
        workflow,
        "adaptive_spectral_blocks",
        lambda *args, **kwargs: ([SpectralBlock(0, "adaptive_common_mode_01", 1, 1, "common")], {}),
    )

    def geometry(coords, **kwargs):
        calls["geometry"] += 1
        if failure == "failed_diffusion" and calls["geometry"] == 2:
            raise RuntimeError("controlled diffusion failure")
        return SimpleNamespace(distance_condensed=pdist(data), metadata={"test": True})

    def gate(features, distances, **kwargs):
        calls["gate"] += 1
        assert features.index.tolist() == data.index.tolist()
        if failure == "failed_gate" and calls["gate"] == 2:
            raise RuntimeError("controlled gate failure")
        return np.repeat([0, 1], 8)

    monkeypatch.setattr(workflow, "block_adaptive_diffusion_geometry", geometry)
    monkeypatch.setattr(workflow, "run_current_kl", gate)
    for name in (
        "plot_axis_terms",
        "plot_combined_axis_terms",
        "plot_dendrogram",
        "plot_diffusion_embedding",
        "plot_subspace_embedding",
        "plot_tree_subtree_clusters",
        "write_ranking_plots",
    ):
        monkeypatch.setattr(workflow, name, lambda *args, **kwargs: None)
    monkeypatch.setattr(
        workflow,
        "write_workflow_pdfs",
        lambda output, *args: {
            "method_tree_pages_pdf": str(output / "method.pdf"),
            "all_tree_pages_pdf": str(output / "all.pdf"),
        },
    )
    workflow.run_analysis(args)
    ranking = pd.read_csv(args.output_dir / "rankings" / f"{workflow.RESULT_PREFIX}_ranking.csv")
    assert set(ranking.status) == {"ok", failure}
    good = ranking[ranking.status.eq("ok")].iloc[0]
    assignments = pd.read_csv(good.cluster_assignments)
    assert assignments.gene.tolist() == data.index.tolist()
    assert assignments.cluster_id.tolist() == [0] * 8 + [1] * 8
    assert np.isfinite(good.go_bic_active_per_gene)
    bad = ranking[ranking.status.eq(failure)].iloc[0]
    assert pd.isna(bad.cluster_assignments)
    assert pd.read_csv(bad.failure_status).status.iloc[0] == failure
    assert (args.output_dir / "connected_results_manifest.json").exists()
