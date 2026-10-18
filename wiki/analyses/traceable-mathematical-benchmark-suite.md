---
title: Traceable Mathematical Benchmark Suite
type: analysis
status: draft
updated: 2026-10-17
sources:
  - raw/inbox/traceable-benchmark-suite-request-20260604.txt
  - benchmarks/shared/cases/method_proof.py
  - benchmarks/diagnostics/math_trace/trace_schema.py
tags:
  - benchmarks
  - traceability
  - calibration
  - validation
---

# Traceable Mathematical Benchmark Suite

## Summary

The method-proof benchmark suite turns Tree-Break Selection benchmark cases into explicit
mathematical probes. Each case declares the contract it is meant to validate or
falsify. `benchmarks/diagnostics/math_trace/trace_schema.py` lists the trace
columns each case requires; the failure classifier that labelled trace rows by
first failed layer was removed on 2026-10-17 because nothing produced the
trace it read.

## Details

The initial `method_proof` suite is registered in
`benchmarks/shared/cases/method_proof.py`. It adds ten cases covering
barycentric balance, selected MP, categorical boundary behavior, continuous
low-rank covariance, support failure, phylogenetic Brownian nulls, and
pass-through replacement.

The current implementation is a scaffold and contract validator. It does not
claim production selected-tail calibration, selected MP validity, or traversal
FDR control. Those claims still depend on held-out precision, support, and
family-specific diagnostics recorded elsewhere in the wiki.

## Evidence

- `raw/inbox/traceable-benchmark-suite-request-20260604.txt` specifies the
  requested benchmark philosophy, case tiers, trace artifacts, failure labels,
  and metamorphic checks.
- `benchmarks/shared/cases/method_proof.py` records mathematical targets and
  required trace columns for the first ten method-proof cases.
- `benchmarks/diagnostics/math_trace/failure_classifier.py` implemented the
  deterministic failure hierarchy. Nothing produced a node-decision trace for
  it, so it was removed on 2026-10-17 (last present at `bc38f16c`).
  `trace_schema.py` keeps the required column list.

## Links

- [[benchmark-pipeline-contract]]
- [[selected-hierarchy-null-support-contract]]
- [[method-proof-web]]
- [[open-mathematical-questions]]

## Open Questions

- Which method-proof cases should be promoted from trace scaffolds to full
  acceptance tests after support and precision thresholds are locked?
- What held-out selected-tail contexts, if any, become production-admissible
  after parent-size and balance transfer checks?
