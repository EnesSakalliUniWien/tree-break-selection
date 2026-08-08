# Tests Overview

The suite is organized by supported product and artifact contracts so you can validate fast
structural behavior first, then statistics, traversal, pipelines, integration, and visualization.

## Directory Structure

```text
tests/
├── conftest.py
├── core/
├── tree/
├── statistics/
├── localization/
├── validation/
│   ├── calibration/
│   └── contracts/
├── pipeline/
├── applications/
├── integration/
├── visualization/
└── wiki/
```

## Current Suite Map

The ordered gate currently collects 635 test cases:

| Stage | Responsibility | Tests |
| --- | --- | ---: |
| 1 | Core structure and decomposition (`core/`, `tree/`) | 140 |
| 2 | Statistical engines and calibration (`statistics/`) | 237 |
| 3 | Localization and post-hoc merge behavior (`localization/`) | 33 |
| 4 | Validation artifact contracts (`validation/`) | 11 |
| 5 | Pipeline and application contracts (`pipeline/`, `applications/`) | 171 |
| 6 | Integration smoke and visualization (`integration/`, `visualization/`) | 39 |
| 7 | Wiki memory contracts (`wiki/`) | 4 |

Stage 4 is intentionally narrow: eight cases protect the diagnostic reporting
and CLI contracts, while three protect the method-constant manifest. Research
simulations, calibration panels, sweeps, and one-off diagnostic studies are
evidence-generating tools rather than production contracts and are not part of
the project test gate.

`conftest.py` is the only shared root test module. Tests otherwise belong to
the directory matching the supported production or artifact responsibility
they protect. Direct tests of production statistical behavior live in
`statistics/`; end-to-end method behavior lives in `integration/`.

## Recommended Test Execution Order

Use the lean development environment for targeted suites:

```bash
uv sync --extra dev --extra benchmark --extra viz --locked
```

Use the full test environment before running the full suite. The full suite
includes scRNA and optional GPL paths backed by `scanpy`, `anndata`, and
`graphtools`.

```bash
uv sync --extra all --extra experimental-gpl --locked
```

Use the helper script when you want the staged order the repository expects:

```bash
uv run python scripts/run_tests_ordered.py --list
uv run python scripts/run_tests_ordered.py --stage 1
uv run python scripts/run_tests_ordered.py
```

Or run the suites directly:

```bash
# 1) Core structure + decomposition
uv run pytest tests/core/ tests/tree/

# 2) Statistical engines + calibration
uv run pytest tests/statistics/

# 3) Traversal and gate behavior
uv run pytest tests/localization/

# 4) Validation artifact contracts
uv run pytest tests/validation/

# 5) Pipeline + application contracts and reporting artifacts
uv run pytest tests/pipeline/ tests/applications/

# 6) Integration smoke + visualization
uv run pytest tests/integration/ tests/visualization/

# 7) Wiki memory contracts
uv run pytest tests/wiki/

# Full suite
uv run pytest
```

## Notes

- Grouping is by intent, not by a strict import dependency graph.
- If you touched decomposition logic, start with `tests/core/`, `tests/statistics/`, and `tests/localization/`.
- For benchmark or artifact-generation changes, prioritize `tests/pipeline/` and `tests/integration/`.
- Research diagnostics belong in `benchmarks/diagnostics/`. Promote a test into
  the project gate only when it protects a supported production or artifact
  contract independent of the diagnostic implementation.
