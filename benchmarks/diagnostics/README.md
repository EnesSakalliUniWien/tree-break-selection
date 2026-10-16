# Benchmark Diagnostics

This directory contains investigation tools for benchmark behavior. These
scripts may read benchmark outputs or simulate calibration diagnostics, but
they are not part of the production clustering method.

Use `benchmarks/shared/` for reusable benchmark execution contracts. Add files
here only when they are diagnostic entrypoints or diagnostic-only helpers.

## Layout

- `calibration/`: responsibility-grouped calibration diagnostics; see
  `calibration/README.md` for the source and mirrored-test map.
- `spectral/`: adaptive-cosine, tree-strategy stability, and Vandermonde diagnostics.
- `math_trace/`: trace-schema validation and deterministic mathematical
  failure attribution for benchmark run artifacts.
- `failure/`: benchmark failure-report tracing used by the full benchmark.
- `analysis/`: post-run result analysis and durable diagnostic notes.

Maintained entrypoints:

- `analysis/analyze_relationships.py`
