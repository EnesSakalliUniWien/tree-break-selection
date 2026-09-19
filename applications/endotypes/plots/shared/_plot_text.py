"""Shared text formatting for GO analysis plots."""

from __future__ import annotations

import textwrap

import pandas as pd


def wrap_labels(values: pd.Series, width: int = 38) -> list[str]:
    return ["\n".join(textwrap.wrap(str(value), width=width)) for value in values]


def wrap_annotation_lines(lines: list[str], *, width: int = 58) -> str:
    wrapped: list[str] = []
    for line in lines:
        wrapped.extend(textwrap.wrap(line, width=width, subsequent_indent="  "))
        wrapped.append("")
    return "\n".join(wrapped).rstrip()
