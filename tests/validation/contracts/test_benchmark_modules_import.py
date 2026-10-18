"""Every benchmark module imports, so removed helpers cannot leave dangling imports."""

import importlib
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[3]
# Module-level scripts that run a benchmark when imported.
SCRIPT_MODULES = {"benchmarks.smoke.run_subset"}
MODULES = sorted(
    ".".join(path.relative_to(REPO_ROOT).with_suffix("").parts).removesuffix(".__init__")
    for path in (REPO_ROOT / "benchmarks").rglob("*.py")
    if "__pycache__" not in path.parts and "results" not in path.parts
)


@pytest.mark.parametrize("module", [m for m in MODULES if m not in SCRIPT_MODULES])
def test_benchmark_module_imports(module):
    importlib.import_module(module)
