"""Rerun canonical `tbs` on the cases that failed in the 2026-10-16 full run.

A case failed when TBS returned `unsupported` or an ARI below 0.5. The run
directory under `benchmarks/results/` is local (gitignored); its comparison
table supplies the case list.
"""

import sys
import warnings
from pathlib import Path

import pandas as pd

from benchmarks.shared.cases import get_test_cases_by_suite
from benchmarks.shared.pipeline import benchmark_cluster_algorithm

warnings.filterwarnings("ignore")

RUN_DIR = Path("benchmarks/results/run_20261016_163716Z_full")
OUT = Path(__file__).parent / "tbs_rerun_failing_cases.csv"

old = pd.read_csv(RUN_DIR / "full_benchmark_comparison.csv")
old = old[old.method == "tbs"].set_index("case_id")
failing = set(old.index[(old.status != "ok") | (old.ari < 0.5)])
cases = [case for case in get_test_cases_by_suite("full") if case["name"] in failing]
print(f"{len(failing)} failing cases; {len(cases)} found in the suite", flush=True)
results = benchmark_cluster_algorithm(
    test_cases=cases,
    methods=["tbs"],
    verbose=False,
    plot_umap=False,
    include_cover_pages=False,
    include_validation_page=False,
)
results.to_csv(sys.argv[1] if len(sys.argv) > 1 else OUT, index=False)
