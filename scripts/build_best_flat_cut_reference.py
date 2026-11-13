"""Build the null table of the best flat cut's ``rho`` for one standard Gaussian population.

Command (repository root, deterministic; the table depends on --datasets only):

    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 uv run --no-sync --offline python \\
        scripts/build_best_flat_cut_reference.py --datasets 19999 --workers 64

``--k-min K`` computes only the cells with k >= K and copies the others from the existing table at ``--out``
(same n grid, B and seed required). Every cell has its own seed, so this gives the same file contents as a full
rebuild; the table still depends on --datasets only.

For every (n, k) it draws ``datasets`` standard Gaussian samples of n cells in k coordinates (seed
``[20261110, 81, n, k, block]``, blocks of ``BLOCK`` datasets), searches each for its best flat cut and keeps
the largest ``ceil(0.1 * datasets)`` values of ``rho``. Cells with k > 10 and n below
``MIN_CELLS_ABOVE_10_DIRECTIONS`` = 50 are not built and stored as NaN (at n = 15 the whitening fails for k = 15
and ``rho`` explodes for k = 11 and 13); the module fails closed there. The search is the one of the
recorded evidence (``raw/inbox/choir-cell-lines/gaussian-best-flat-cut-reference-scripts-20261111.py.txt``, ``lib.search``
with top = 50 and single-point moves) and is checked against brute force in
``tests/statistics/82_test_best_flat_cut_reference.py``.
"""

from __future__ import annotations

import argparse
import math
import time
from functools import lru_cache
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.best_flat_cut_reference import (
    MIN_CELLS_ABOVE_10_DIRECTIONS,
    TABLE_PATH,
    to_rho,
    whiten,
)

SEED = 20261110
N_GRID = (15, 20, 30, 40, 50, 75, 100, 150, 200, 300, 400, 800)
K_GRID = tuple(range(1, 16))
BLOCK = 100
KEEP_FRACTION = 0.1


def grid(k: int, m: int, stream: int = 0) -> np.ndarray:
    """m unit directions; k = 1 exact, k = 2 equal angles, k = 3 Fibonacci sphere, k >= 4 random."""
    if k == 1:
        return np.ones((1, 1))
    if k == 2:
        a = np.linspace(0, np.pi, m, endpoint=False)
        return np.c_[np.cos(a), np.sin(a)]
    if k == 3:
        i = np.arange(m) + 0.5
        z, phi = i / m, np.pi * (1 + 5**0.5) * i
        r = np.sqrt(1 - z * z)
        return np.c_[r * np.cos(phi), r * np.sin(phi), z]
    v = np.random.default_rng([SEED, 60, k, m, stream]).standard_normal((m, k))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


@lru_cache(maxsize=None)
def directions(k: int) -> np.ndarray:
    """Search directions: the recorded fixed grids for k <= 6 and ``grid(k, 3000, stream=7)`` above."""
    if k == 1:
        return np.ones((1, 1))
    if k == 2:
        a = np.linspace(0, np.pi, 720, endpoint=False)
        return np.c_[np.cos(a), np.sin(a)]
    if k == 3:
        i = np.arange(2000) + 0.5
        z, phi = i / 2000, np.pi * (1 + 5**0.5) * i
        r = np.sqrt(1 - z * z)
        return np.c_[r * np.cos(phi), r * np.sin(phi), z]
    if k <= 6:
        v = np.random.default_rng([SEED, 7, k]).standard_normal((3000, k))
        return v / np.linalg.norm(v, axis=1, keepdims=True)
    return grid(k, 3000, stream=7)


def _along_one(y, v):
    n = len(y)
    order = np.argsort(y @ v)
    prefix = np.cumsum(y[order], axis=0)[:-1]
    s = np.arange(1, n)
    r2 = (prefix**2).sum(-1) / (s * (n - s))
    i = int(np.argmax(r2))
    mask = np.zeros(n, bool)
    mask[order[: i + 1]] = True
    return r2[i], mask


def resort(y, mask, r2, steps=20):
    """Fisher re-sorts: sort along the mean difference, take the best threshold."""
    for _ in range(steps):
        d = y[mask].mean(0) - y[~mask].mean(0)
        r2_new, m_new = _along_one(y, d / np.linalg.norm(d))
        if r2_new <= r2 + 1e-12:
            break
        r2, mask = r2_new, m_new
    return r2, mask


def swap(y, mask, r2):
    """Single-point moves while any move raises r2; then re-sort; repeat until stable."""
    n = len(y)
    while True:
        improved = False
        while True:
            S, s = y[mask].sum(0), mask.sum()
            sign = np.where(mask, -1.0, 1.0)
            S_new = S + sign[:, None] * y
            s_new = s + sign
            ok = (s_new >= 1) & (s_new <= n - 1)
            cand = np.where(ok, (S_new**2).sum(1) / np.maximum(s_new * (n - s_new), 1), -1)
            i = int(np.argmax(cand))
            if cand[i] <= r2 + 1e-12:
                break
            mask = mask.copy()
            mask[i] = ~mask[i]
            r2, improved = cand[i], True
        r2b, maskb = resort(y, mask, r2)
        if r2b > r2 + 1e-12:
            r2, mask, improved = r2b, maskb, True
        if not improved:
            return r2, mask


def search(x, dirs, top=1, local_swap=False, chunk=2000):
    """``rho`` of the best flat cut found: thresholds on the grid directions, re-sorts, single-point moves."""
    y = whiten(x)
    n = len(y)
    s = np.arange(1, n)[:, None]
    best_r2 = np.empty(len(dirs))
    for a in range(0, len(dirs), chunk):
        g = dirs[a : a + chunk]
        order = np.argsort(y @ g.T, axis=0)
        prefix = np.cumsum(y[order], axis=0)[:-1]
        best_r2[a : a + chunk] = ((prefix**2).sum(-1) / (s * (n - s))).max(0)
    picks = np.argsort(-best_r2, kind="stable")[:top]
    out = 0.0
    for j in picks:
        r2, mask = _along_one(y, dirs[j])
        r2, mask = resort(y, mask, r2)
        if local_swap:
            r2, mask = swap(y, mask, r2)
        out = max(out, r2)
    return to_rho(out, n)


def best_flat_cut_rho(x: np.ndarray) -> float:
    """``rho`` of the best flat cut of the sample ``x`` (n x k) as stored in the table."""
    return search(x, directions(x.shape[1]), top=50, local_swap=True)


def block_job(job: tuple[int, int, int, int]) -> tuple[int, int, int, list[float]]:
    n, k, block, size = job
    rng = np.random.default_rng([SEED, 81, n, k, block])
    return n, k, block, [best_flat_cut_rho(rng.standard_normal((n, k))) for _ in range(size)]


def built(n: int, k: int) -> bool:
    """Whether the cell (n, k) is computed; the others are NaN."""
    return k <= 10 or n >= MIN_CELLS_ABOVE_10_DIRECTIONS


def build(datasets: int, workers: int, out: Path, k_min: int = 1) -> None:
    keep = math.ceil(KEEP_FRACTION * datasets)
    upper = np.full((len(N_GRID), len(K_GRID), keep), np.nan, dtype=np.float32)
    if k_min > 1:
        with np.load(out) as old:
            same = (
                old["n_grid"].tolist() == list(N_GRID)
                and old["k_grid"][: k_min - 1].tolist() == list(range(1, k_min))
                and int(old["B"]) == datasets
                and int(old["seed"]) == SEED
            )
            if not same:
                raise SystemExit(f"{out}: n grid, k grid, B or seed differ; rebuild with --k-min 1")
            upper[:, : k_min - 1] = old["upper"][:, : k_min - 1]
    blocks = [(b, min(BLOCK, datasets - b * BLOCK)) for b in range(math.ceil(datasets / BLOCK))]
    jobs = [
        (n, k, b, size)
        for n in N_GRID
        for k in K_GRID
        if k >= k_min and built(n, k)
        for b, size in blocks
    ]
    jobs.sort(key=lambda job: -job[0] * job[1])  # slowest first
    samples: dict[tuple[int, int], list[float]] = {}
    started = time.time()
    with Pool(workers) as pool:
        for done, (n, k, _block, rho) in enumerate(pool.imap_unordered(block_job, jobs), 1):
            samples.setdefault((n, k), []).extend(rho)
            if done % 50 == 0:
                print(f"{done}/{len(jobs)} blocks, {time.time() - started:.0f} s", flush=True)
    for i, n in enumerate(N_GRID):
        for j, k in enumerate(K_GRID):
            if (n, k) not in samples:
                continue
            values = np.sort(np.asarray(samples[n, k]))
            assert len(values) == datasets
            upper[i, j] = values[-keep:]
    np.savez_compressed(
        out,
        upper=upper,
        n_grid=np.array(N_GRID),
        k_grid=np.array(K_GRID),
        B=np.array(datasets),
        seed=np.array(SEED),
    )
    print(f"wrote {out} after {time.time() - started:.0f} s", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--datasets", type=int, default=19999, help="datasets per (n, k) cell, B")
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--out", type=Path, default=TABLE_PATH)
    parser.add_argument(
        "--k-min", type=int, default=1, help="copy the cells with smaller k from --out"
    )
    args = parser.parse_args()
    build(args.datasets, args.workers, args.out, args.k_min)
