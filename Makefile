.PHONY: audit check lint test wiki-lint

check: lint audit wiki-lint test

lint:
	uv run --no-sync ruff check .

audit:
	uv run --no-sync deptry applications tree_break_selection benchmarks
	uv run --no-sync vulture applications benchmarks scripts tree_break_selection --min-confidence 90

test:
	uv run --no-sync python scripts/run_tests_ordered.py

wiki-lint:
	python3 scripts/wiki/lint.py
