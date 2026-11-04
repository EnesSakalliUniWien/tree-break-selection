# Maintained conditional sampler for root decisions

2026-11-04. Implements the root-split random-position chain from the
[candidate comparison](../selective_rebuild_comparison_20261104/README.md) as
an opt-in diagnostic API. It does not replace default production TBS inference.

## Entry point and result

```python
from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_conditional_mc import (
    selected_diffusion_root_chain_test,
)

result = selected_diffusion_root_chain_test(
    values, sigma=known_noise_standard_deviation, steps=999, seed=17, alpha=0.05
)
print(result["decision"], result["p_value"])
```

The input must consist of independent Gaussian rows with known common noise
covariance `sigma**2 * I` in the supplied features. The API requires `sigma`
and does not estimate it. The caller must justify this model; the array checks
do not verify Gaussianity, independence or absence of fitted preprocessing.
All input features enter the contrast; there is no MP/adaptive-rank or branch
length correction in this candidate statistic.

One native diffusion/average tree defines the observed root children. Native
geometry is rebuilt along the Gaussian contrast path for simulation proposals.
The event fixes those children while permitting other merges to change. The
API keeps the canonical labels, observed merge signature and source hashes.
It does not run a separate tree on each child or decide descendant splits.

`decision` is `reject`, `do_not_reject`, or `unresolved`. Inspect it together
with `status`, `sampling_status`, `accepted_moves`, `unique_states` and
`resolution_sufficient`. A stalled chain has p=1 and an unresolved decision;
insufficient rank resolution also gives an unresolved decision. Failed
observation or proposal geometry gives no p-value and reports its phase.
`do_not_reject` does not establish that the two group means are equal.

The full chain and random observation position are retained as `chain_radii`
and `anchor_index`. Geometry errors invalidate the calculation, instead of
silently becoming rejected proposals. Native floating-point tie breaking
defines proposal membership; `numerical_event_certified` remains false.
The existing observation constructor also requires a supported strict replay.

## Statistical scope

The unit Gaussian proposal step and chain length are fixed, and the observation
occupies an independent uniform position on the path. Reversibility justifies
the conservative upper rank under the exact selected Gaussian reference, with
all ties included. There is no burn-in or convergence-based stopping rule.
Lack of mixing can reduce detection. Neither a low p-value nor a high move
count establishes exact-real membership or applicability to the cancer data.
Do not choose seeds or chain lengths by repeatedly checking significance.

The rejection-sampling API remains available and its existing tests pass.
Both methods use shared native replay, which temporarily instruments module
functions and therefore must run serially within a process.

## Validation and reproduction

All 86 runs completed without unresolved decisions. The 24 parity runs match
the saved candidate chain arrays and ranks exactly. Fresh-case results:

| Scenario | Datasets | Rejections at 0.05 | Rejections at 0.01 |
| --- | ---: | ---: | ---: |
| No group difference | 40 | 3 | 0 |
| Moderate signal | 10 | 4 | 1 |
| Strong signal | 10 | 7 | 5 |

The registered 60-by-40 clear Gaussian case gives p=0.005 with sigma=0.6;
the 45-by-30 moderate Gaussian case gives p=0.001 with sigma=1.2. These are
root decisions only. All covariance scales were known from the generators.
Three false alarms among 40 nulls neither establishes nor disproves calibration
at 0.05; the empirical precision is limited. Signals share paired noise across
scenarios, and the earlier development cases are not counted as fresh evidence.

Validation recomputed all 86 ranks, checked 258 native chain memberships and
149 source hashes. Fresh-case chains made at least 44 accepted moves; movement
is a diagnostic, not a convergence certificate. Saved results and timings are
in [results.csv](results.csv), with [summary.csv](summary.csv) and
[verification.json](verification.json). No failures were discarded or rerun.

The [protocol](protocol.md) fixes parity, fresh-seed and registered-case checks.
The maintained implementation first failed 28 new tests because the API and
helpers were absent. After implementation, all 44 focused replay/sampler tests
passed. The full ordered suite passed 1,031 tests across all seven stages,
with one unrelated scikit-learn MDS `n_init` FutureWarning.

```bash
uv run --no-sync --offline pytest -q tests/statistics/61_test_diffusion_root_replay.py tests/statistics/67_test_diffusion_root_conditional_mc.py tests/statistics/68_test_diffusion_root_conditional_chain.py
make test
uv run --no-sync --offline python reports/root_conditional_sampler_20261104/validate_sampler.py
```

The validation runner requires a fresh `generated/records.jsonl` destination
and refuses to overwrite its completed records. It saves each input, full
result, source hashes and compact tables. Reproduction starts from the saved
prior comparison arrays for the parity phase; those are catalogued local
outputs, not promised external backups.

## Remaining limits

This implementation tests the selected root's average-mean contrast only.
Unknown covariance, fitted-PC/count preprocessing, correlated rows,
descendant/multiple-testing decisions, full cancer-data calibration and
exact-real event certification remain outside the supported model. No
dependencies were installed, production default changed or commit made.
