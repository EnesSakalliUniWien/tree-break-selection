"""Application scripts default to fresh results folders, never retained raw/ evidence."""

from applications.output_dirs import PROJECT_ROOT, timestamped_results_dir


def test_timestamped_results_dir_is_under_results_analyses():
    path = timestamped_results_dir("demo_report")

    assert path.parent == PROJECT_ROOT / "results" / "analyses"
    assert path.name.startswith("demo_report_")
    assert "raw" not in path.relative_to(PROJECT_ROOT).parts
