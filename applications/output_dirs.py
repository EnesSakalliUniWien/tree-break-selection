"""Default output folders for application scripts.

New results go to ``results/analyses/<slug>_<YYYYMMDD_HHMMSS>/``. Retained evidence
under ``raw/assets/`` is written only when the caller passes it as ``--output-dir``.
"""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def timestamped_results_dir(slug: str) -> Path:
    """Return a new ``results/analyses/<slug>_<timestamp>`` path under the project root."""
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return PROJECT_ROOT / "results" / "analyses" / f"{slug}_{stamp}"


def parse_output_dir(description: str | None) -> Path | None:
    """Parse the single ``--output-dir`` option of a report script."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Default: a new results/analyses/<report>_<timestamp>/ folder.",
    )
    return parser.parse_args().output_dir


def parse_run_dir(description: str | None) -> Path:
    """Parse the required ``--output-dir`` of a script that reads and writes one run in place."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Benchmark run directory to read and write; pass a raw/assets/ path explicitly.",
    )
    return parser.parse_args().output_dir


__all__ = ["PROJECT_ROOT", "parse_output_dir", "parse_run_dir", "timestamped_results_dir"]
