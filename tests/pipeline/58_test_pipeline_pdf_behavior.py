from pathlib import Path

from benchmarks.shared import pipeline


def test_concat_pdf_streams_and_no_pngs(monkeypatch, tmp_path: Path):
    called = {}

    def fake_generate_benchmark_plots(df_results, computed_results, **kwargs):
        called.update(kwargs)

    monkeypatch.setattr(pipeline, "generate_benchmark_plots", fake_generate_benchmark_plots)

    # With concat_plots_pdf=True, pipeline streams to PdfPages and never emits PNGs.
    pipeline.benchmark_cluster_algorithm(
        test_cases=[],
        verbose=True,
        concat_plots_pdf=True,
        methods=["tbs"],
    )

    assert called.get("include_cover_pages") is True
    assert called.get("include_validation_page") is True
    assert called.get("pdf") is not None


def test_pipeline_preserves_explicit_case_identity(monkeypatch):
    seen_case_numbers: list[int] = []

    def fake_run_single_case(**kwargs):
        seen_case_numbers.append(kwargs["tc"]["test_case_num"])
        return [], []

    monkeypatch.setattr(pipeline, "run_single_case", fake_run_single_case)
    monkeypatch.setattr(
        pipeline,
        "generate_benchmark_plots",
        lambda *args, **kwargs: None,
    )

    pipeline.benchmark_cluster_algorithm(
        test_cases=[
            {"name": "case_a", "seed": 1, "test_case_num": 41},
            {"name": "case_b", "seed": 2, "test_case_num": 42},
        ],
        verbose=False,
        methods=["tbs"],
    )

    assert seen_case_numbers == [41, 42]
