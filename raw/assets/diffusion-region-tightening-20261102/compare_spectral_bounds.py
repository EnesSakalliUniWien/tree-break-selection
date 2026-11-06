"""Captured same-box before/after diagnostic; run from repository root with uv."""

import gzip
import hashlib
import json
import sys
import types
from pathlib import Path

import numpy as np

from benchmarks.diagnostics.calibration.selected.hierarchy import diffusion_root_cells as new
from benchmarks.diagnostics.calibration.selected.hierarchy import diffusion_root_kernel as kernel
from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_mass import (
    chi_interval_mass,
)
from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_replay import (
    make_diffusion_root_path,
)

out = Path(__file__).resolve().parent
old = types.ModuleType("captured_gram_before")
sys.modules[old.__name__] = old
exec(
    compile(
        gzip.decompress((out / "gram_engine_before.py.gz").read_bytes()),
        "<captured prior Gram engine>",
        "exec",
    ),
    old.__dict__,
)
middle = types.ModuleType("captured_global_factor_before")
sys.modules[middle.__name__] = middle
exec(
    compile(
        gzip.decompress((out / "factor_global_before.py.gz").read_bytes()),
        "<captured prior global factor engine>",
        "exec",
    ),
    middle.__dict__,
)
files = [
    Path(new.__file__),
    Path(kernel.__file__),
    out / "gram_engine_before.py.gz",
    out / "factor_global_before.py.gz",
]
hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
path = make_diffusion_root_path(np.random.default_rng(20261105).normal(size=(9, 1)))
rows = []
for width in [1e-12, 1e-10, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3]:
    lo, hi = path.observed_radius - width, path.observed_radius + width
    try:
        box = kernel.enclose_diffusion_operator(path, lo, hi)
        midpoint = (np.array(box.lower, dtype=float) + np.array(box.upper, dtype=float)) / 2
    except ValueError as error:
        rows.append(dict(half_width=width, kernel_status="unresolved", reason=str(error)))
        continue
    row = dict(half_width=width, kernel_status="enclosed", same_operator_box=True)
    for name, engine in [
        ("before_gram", old),
        ("intermediate_global_factor", middle),
        ("after_restricted_factor", new),
    ]:
        try:
            reference = engine._spectral_reference(midpoint, box.retained_count_target, box.time)
            lower, upper, details = engine._distance_enclosures(reference, box.lower, box.upper)
            status, margins = engine._classify_linkage(
                lower, upper, path.observed_signature, path.tie_tolerance
            )
            row[name] = dict(
                status=status,
                certified_selected_mass_lower=(
                    str(chi_interval_mass(lo, hi, path.dimension)[0])
                    if status == "selected"
                    else "0"
                ),
            )
        except ValueError as error:
            row[name] = dict(status="unresolved", reason=str(error))
    rows.append(row)
    print(
        width,
        row["before_gram"]["status"],
        row["intermediate_global_factor"]["status"],
        row["after_restricted_factor"]["status"],
        flush=True,
    )
result = dict(
    fixture_seed=20261105,
    n=9,
    p=1,
    observed_radius=path.observed_radius,
    source_sha256=hashes,
    contract="same exact-real operator box for both spectral methods",
    sources_unchanged_during_run=all(
        hashlib.sha256(Path(p).read_bytes()).hexdigest() == h for p, h in hashes.items()
    ),
    rows=rows,
)
(out / "spectral_postkernel_comparison.json").write_text(json.dumps(result, indent=2) + "\n")
