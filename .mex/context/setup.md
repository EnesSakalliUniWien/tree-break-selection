---
name: setup
description: Setup pointers for mex; root project files remain authoritative.
triggers:
  - "setup"
  - "install"
  - "run"
  - "environment"
edges:
  - target: ../pyproject.toml
    condition: dependencies, optional groups, pytest, and ruff config
  - target: ../Makefile
    condition: check, lint, audit, test and wiki-lint targets
  - target: ../README.md
    condition: project overview and implementation map
last_updated: 2026-10-15
---

# Setup

## Prerequisites

- Python 3.11 or newer.
- `uv`; dependencies come from `pyproject.toml` and the locked `uv.lock`.

## Common Commands

- `uv sync --extra dev --extra benchmark --extra viz --locked` - locked development environment.
- `make check` - lint, dependency and dead-code audit, wiki lint, generated-file check and the staged test suite.
- `make test` - staged test suite (`scripts/run_tests_ordered.py`).
- `uv run pytest tests/wiki/test_memory_contract.py -q` - memory and lookup contract tests.
- `make wiki-lint` - structural wiki validation.
- `make lint` - `ruff check .` and `ruff format --check .`.

## Common Issues

- Mex scanner output currently omits Python dependency details; read `pyproject.toml`.
- Wiki lint failures usually mean a source citation no longer points to existing local evidence.
