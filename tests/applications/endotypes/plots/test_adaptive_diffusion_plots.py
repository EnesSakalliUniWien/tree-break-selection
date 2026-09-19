"""In-memory plotting and DataFrame I/O contracts for GO-IC reports."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from applications.endotypes.analysis.go_ic.adaptive_diffusion_go_ic_terms import (
    build_axis_loadings,
    build_term_metadata,
)
from applications.endotypes.io.adaptive_diffusion_dataframe_io import (
    read_dataframe,
    write_dataframe,
)
from applications.endotypes.plots.embeddings.diffusion_embedding_plot import (
    create_diffusion_embedding,
)
from applications.endotypes.plots.embeddings.embedding_comparison_plot import (
    create_embedding_comparison,
)
from applications.endotypes.plots.embeddings.subspace_embedding_plot import (
    create_subspace_embedding,
)
from applications.endotypes.plots.ranking.subspace_ranking_plot import create_subspace_ranking
from applications.endotypes.plots.terms.axis_term_bar_plot import create_axis_term_bars
from applications.endotypes.plots.terms.axis_term_heatmap import create_axis_term_heatmap
from applications.endotypes.plots.trees.dendrogram_plot import create_dendrogram
from applications.endotypes.plots.trees.tree_cluster_plot import create_tree_clusters
from applications.endotypes.reports.go_ic.adaptive_diffusion_plot_exports import write_figure
from matplotlib.figure import Figure
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist


@pytest.fixture
def loadings():
    data = pd.DataFrame([[1, 0, 1], [1, 1, 0], [0, 1, 1]], columns=["A", "B", "C"])
    return build_axis_loadings(
        term_metadata=build_term_metadata(data),
        loadings=np.array([[0.7, -0.4], [-0.3, 0.8], [0.1, 0.2]]),
        block_start=1,
        block_end=2,
        top_n=2,
    )


def no_io(*args, **kwargs):
    raise AssertionError("Computation or plotting attempted file I/O")


def test_loading_tables_are_pure_and_do_not_mutate_metadata(monkeypatch):
    monkeypatch.setattr(pd.DataFrame, "to_csv", no_io)
    monkeypatch.setattr(pd, "read_csv", no_io)
    metadata = build_term_metadata(pd.DataFrame([[1, 0], [0, 1]], columns=["A", "B"]))
    original = metadata.copy(deep=True)
    all_terms, top_terms = build_axis_loadings(
        term_metadata=metadata,
        loadings=np.array([[0.2, -0.8], [-0.5, 0.1]]),
        block_start=1,
        block_end=2,
        top_n=1,
    )
    pd.testing.assert_frame_equal(metadata, original)
    assert len(all_terms) == 4
    assert top_terms[top_terms.selection.eq("top_absolute")].column.tolist() == ["B", "A"]
    assert list(all_terms.axis.unique()) == [1, 2]


@pytest.mark.parametrize(
    "kind", ["bars", "heatmap", "subspace", "diffusion", "dendrogram", "tree_clusters", "ranking"]
)
def test_each_plot_type_returns_a_figure_without_file_io(kind, loadings, monkeypatch):
    monkeypatch.setattr(pd.DataFrame, "to_csv", no_io)
    monkeypatch.setattr(pd, "read_csv", no_io)
    monkeypatch.setattr(Figure, "savefig", no_io)
    all_terms, _ = loadings
    coords = np.array([[0.0], [1.0], [3.0], [5.0]])
    labels = np.array([0, 0, 1, 1])
    distances = pdist(coords)
    z = linkage(distances, method="average")
    assignments = pd.DataFrame({"gene": ["A", "B", "C", "D"], "cluster_id": labels})
    ranking = pd.DataFrame(
        [
            dict(
                status="ok",
                display_rank=1,
                weighting="binary",
                block_name="test",
                go_bic_active_per_gene=2.0,
                quality_tier=0,
                quality_tier_label="quality_plausible",
                specificity_score=0.7,
            )
        ]
    )
    creators = {
        "bars": lambda: create_axis_term_bars(all_terms[all_terms.axis.eq(1)], "terms", top_n=2),
        "heatmap": lambda: create_axis_term_heatmap(all_terms),
        "subspace": lambda: create_subspace_embedding(coords, labels, "subspace"),
        "diffusion": lambda: create_diffusion_embedding(distances, labels, "diffusion"),
        "dendrogram": lambda: create_dendrogram(z, "tree"),
        "tree_clusters": lambda: create_tree_clusters(z, assignments, "clusters"),
        "ranking": lambda: create_subspace_ranking(ranking),
    }
    original = all_terms.copy(deep=True)
    fig = creators[kind]()
    try:
        assert isinstance(fig, Figure)
        fig.canvas.draw()
        pd.testing.assert_frame_equal(all_terms, original)
        if kind == "tree_clusters":
            expected = labels[dendrogram(z, no_plot=True)["leaves"]]
            np.testing.assert_array_equal(fig.axes[1].images[0].get_array()[0], expected)
    finally:
        plt.close(fig)


def test_heatmap_png_and_pdf_share_data_and_color_limits(loadings):
    all_terms, _ = loadings
    png = create_axis_term_heatmap(all_terms, terms_per_axis=2, max_terms=3)
    pdf = create_axis_term_heatmap(all_terms, terms_per_axis=2, max_terms=3, report=True)
    try:
        np.testing.assert_array_equal(
            png.axes[0].images[0].get_array(), pdf.axes[0].images[0].get_array()
        )
        assert png.axes[0].images[0].get_clim() == pdf.axes[0].images[0].get_clim()
        assert [tick.get_text() for tick in png.axes[0].get_yticklabels()] == [
            tick.get_text() for tick in pdf.axes[0].get_yticklabels()
        ]
    finally:
        plt.close(png)
        plt.close(pdf)


def test_empty_heatmap_preserves_export_and_report_behavior():
    empty = pd.DataFrame(columns=["axis"])
    assert create_axis_term_heatmap(empty) is None
    report = create_axis_term_heatmap(empty, report=True)
    try:
        assert report.axes[0].texts[0].get_text() == "No GO-term loading data available"
    finally:
        plt.close(report)


def test_dataframe_io_roundtrip(tmp_path):
    frame = pd.DataFrame({"gene": ["A", "B"], "cluster_id": [0, 1]})
    path = tmp_path / "assignments.csv"
    write_dataframe(frame, path)
    pd.testing.assert_frame_equal(read_dataframe(path), frame)


def test_figure_export_closes_figure_when_write_fails(monkeypatch, tmp_path):
    fig = plt.figure()

    def fail(*args, **kwargs):
        raise OSError("unwritable")

    monkeypatch.setattr(fig, "savefig", fail)
    with pytest.raises(OSError, match="unwritable"):
        write_figure(fig, tmp_path / "plot.png", dpi=180)
    assert not plt.fignum_exists(fig.number)


@pytest.mark.parametrize("labels", [[0, 0, 1, 1], [2, 2, 7, 7], [0, 0, 0, 0]])
def test_tree_and_embedding_share_discrete_colors_and_exact_cluster_keys(labels):
    labels = np.array(labels)
    coords = np.arange(4.0)[:, None]
    tree = create_tree_clusters(
        linkage(pdist(coords)), pd.DataFrame({"cluster_id": labels}), "tree"
    )
    scatter = create_subspace_embedding(coords, labels, "embedding")
    try:
        tree.canvas.draw()
        scatter.canvas.draw()
        image = tree.axes[1].images[0]
        points = scatter.axes[0].collections[0]
        np.testing.assert_allclose(image.cmap(image.norm(labels)), points.get_facecolors())
        expected = [str(cid) for cid in np.unique(labels)]
        for fig in (tree, scatter):
            legends = [ax.get_legend() for ax in fig.axes if ax.get_legend() is not None]
            assert len(legends) == 1
            assert [text.get_text() for text in legends[0].get_texts()] == expected
    finally:
        plt.close(tree)
        plt.close(scatter)


@pytest.mark.parametrize("mode_count", [6, 39])
def test_comparison_preserves_long_terms_as_readable_native_text(monkeypatch, mode_count):
    monkeypatch.setattr(pd, "read_csv", no_io)
    monkeypatch.setattr(plt, "imread", no_io)
    lines = [
        f"mode {i:02d}: " + "regulation of cellular metabolic process " * (5 if mode_count == 6 else 1)
        for i in range(mode_count)
    ]
    lines[-1] = (
        f"mode {mode_count - 1:02d}: RNA Polymerase II Transcription Regulatory Region Sequence-Specific DNA "
        "Binding (+0.139); Sequence-Specific DNA Binding (+0.138)"
    )
    image = np.ones((130, 160, 3))
    fig = create_embedding_comparison(image, image, "comparison", lines)
    try:
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        texts = [text for ax in fig.axes for text in ax.texts]
        assert all(text.get_fontsize() >= 12 for text in texts)
        for text in texts:
            box = text.get_window_extent(renderer)
            bounds = text.axes.get_window_extent(renderer)
            assert box.x0 >= bounds.x0 - 1 and box.x1 <= bounds.x1 + 1
            assert box.y0 >= bounds.y0 - 1 and box.y1 <= bounds.y1 + 1
        content = " ".join(text.get_text() for text in texts)
        for i in range(mode_count):
            assert f"mode {i:02d}:" in content
    finally:
        plt.close(fig)


def test_ranking_shows_tiers_and_specificity_without_reordering_by_bic():
    ranking = pd.DataFrame([
        dict(status="ok", display_rank=2, weighting="tfidf", block_name="b",
             quality_tier=1, quality_tier_label="broad_coherent",
             specificity_score=0.8, go_bic_active_per_gene=4.0),
        dict(status="ok", display_rank=1, weighting="binary", block_name="a",
             quality_tier=0, quality_tier_label="quality_plausible",
             specificity_score=0.5, go_bic_active_per_gene=8.0),
    ])
    original = ranking.copy(deep=True)
    fig = create_subspace_ranking(ranking)
    try:
        fig.canvas.draw()
        tier, specificity, bic = fig.axes
        assert [t.get_text() for t in tier.texts] == ["0  quality_plausible", "1  broad_coherent"]
        assert [p.get_width() for p in specificity.patches] == [0.5, 0.8]
        assert [p.get_width() for p in bic.patches] == [8.0, 4.0]
        assert [t.get_text() for t in specificity.texts] == ["0.500", "0.800"]
        pd.testing.assert_frame_equal(ranking, original)
    finally:
        plt.close(fig)
