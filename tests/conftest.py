import importlib
import importlib.util

import pytest


@pytest.fixture(scope="module")
def require_optional_dependencies():
    """Import optional modules, skipping only when they are not installed.

    An installed module whose import fails is a broken optional environment and
    must remain visible as a test failure instead of being silently skipped.
    """

    def require(*module_names: str) -> tuple[object, ...]:
        missing = [name for name in module_names if importlib.util.find_spec(name) is None]
        if missing:
            pytest.skip(f"Missing optional dependencies: {', '.join(missing)}")
        return tuple(importlib.import_module(name) for name in module_names)

    return require
