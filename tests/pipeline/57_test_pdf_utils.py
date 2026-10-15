from pathlib import Path

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt
from pypdf import PdfReader

from benchmarks.shared.plots.cover_page import write_case_manifest_pages_to_pdf
from benchmarks.shared.util.pdf.merge import merge_existing_pdfs


def test_merge_existing_pdfs_uses_declared_python_dependency(tmp_path: Path):
    first = tmp_path / "first.pdf"
    second = tmp_path / "second.pdf"
    output = tmp_path / "merged.pdf"

    for path, title in ((first, "First"), (second, "Second")):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title(title)
        fig.savefig(path)
        plt.close(fig)

    assert merge_existing_pdfs([first, second], output, verbose=False)
    assert output.exists()
    assert len(PdfReader(str(output)).pages) == 2


def test_case_manifest_pages_record_generation_contract(tmp_path: Path):
    output = tmp_path / "manifest.pdf"
    cases = [
        {
            "name": "gaussian_case",
            "category": "improved_gaussian",
            "generator": "blobs",
            "n_samples": 12,
            "n_features": 8,
            "n_clusters": 3,
            "cluster_std": 0.5,
        },
        {
            "name": "binary_case",
            "category": "improved_binary_low_noise",
            "generator": "binary",
            "n_samples": 20,
            "n_features": 10,
            "noise_features": 2,
            "n_clusters": 4,
            "entropy_param": 0.1,
        },
    ]

    write_case_manifest_pages_to_pdf(cases, output, rows_per_page=20)

    assert output.exists()
    assert len(PdfReader(str(output)).pages) == 1
