from __future__ import annotations

import numpy as np
import pandas as pd
from benchmarks.shared.types import MethodRunResult


def _toy_dataframe() -> pd.DataFrame:
    return pd.DataFrame(
        np.array(
            [
                [0.0, 0.0],
                [0.1, 0.0],
                [1.0, 1.0],
                [0.9, 1.0],
            ],
            dtype=float,
        ),
        columns=["f0", "f1"],
    )


def _successful_method_result() -> MethodRunResult:
    return MethodRunResult(
        labels=np.array([0, 0, 1, 1], dtype=int),
        found_clusters=2,
        report_df=None,
        status="ok",
        skip_reason=None,
        extra={},
    )


def _capturing_runner(captured: dict[str, object], *, include_args: bool = True):
    def runner(*args, **kwargs):
        if include_args:
            captured["args"] = args
        captured["kwargs"] = kwargs
        return _successful_method_result()

    return runner

def _appending_runner(captured_calls: list[dict[str, object]]):
    def runner(*args, **kwargs):
        captured_calls.append({"args": args, "kwargs": kwargs})
        return _successful_method_result()

    return runner
