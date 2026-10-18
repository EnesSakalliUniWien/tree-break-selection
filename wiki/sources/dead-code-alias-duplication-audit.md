---
title: Dead Code Alias and Duplication Audit
type: source
status: reviewed
updated: 2026-10-17
sources:
  - reports/dead_code_alias_duplication_20260928/README.md
  - reports/dead_code_alias_duplication_20260928/audit_commands.json
  - reports/dead_code_alias_duplication_20260928/ast_candidates.json
  - Makefile
  - pyproject.toml
  - benchmarks/diagnostics/calibration/sibling/known_branch_skew_diagnostic.py
  - benchmarks/validation/sweeps/alpha_grid_search.py
  - benchmarks/validation/statistics/traversal_sibling_fdr_null.py
  - applications/scrna/plots/selected_nnls_fit_summary.py
tags:
  - code-audit
  - duplication
  - aliases
---

# Dead Code Alias and Duplication Audit

## Summary

The audit identified and verified existing repository tools
before auditing dead code, aliases and duplication. Configured deptry, Vulture
and Ruff checks passed. Supplemental static analysis found one repository-
unused private helper and two actionable groups of live duplicate helpers.
The initial audit made no runtime changes. The user subsequently authorized
removal of its three assignment aliases; their callers now use the original
constants directly, preserving defaults and output columns.

## Key Points

- The configured Vulture threshold is 90; it reported no findings. A separate
  threshold-60 scan produced 62 candidate lines, including dynamic-registry,
  serialized-field and pytest-hook false positives. Candidates are not proof
  of deadness.
- Custom AST analysis covered 623 Python files. Thirteen exact-body groups
  were candidates, not thirteen approved refactors. No clone checker is
  configured in the current manifest/lock/scripts.
- `_node_pvalue` in the known-branch diagnostic has no observed calls,
  exports, decorators or dynamic lookup in the inspected repository paths.
  Its removal is a separate cleanup candidate.
- Two `current_git_state` implementations and four timestamp renderers have
  identical bodies and live callers. Consolidation must preserve their
  manifest/error and rendering contracts.
- The three inspected assignment aliases were active, then explicitly removed
  at the user's request: `TOPOLOGIES` became `TREE_CONSENSUS_STRATEGIES`,
  `BASE_CONTEXT_COLUMNS` became `TAIL_LAW_CONTEXT_COLUMNS`, and
  `CANDIDATE_METHOD_CONTRAST_COLUMNS` became `METHOD_CONTRAST_COLUMNS` at all
  executable use sites. Dynamic method-registry runners and the optional
  tested bootstrap-consensus API remain.
- Existing dirty/untracked work was included. Historical provenance, vendor,
  ignored workspace and mirrored fixture duplication were excluded from
  production cleanup findings; external consumers remain outside the proof.

## Evidence

The report retains exact configured commands, versions, exit codes, outputs,
the supplemental AST algorithm description and candidates. Retained
locations and usages were verified. This is a bounded static audit,
not a whole-program reachability proof or an exhaustive alias census.
The removal follow-up passed three focused tests, changed-file Ruff, import/
default/schema checks and AST comparison against the pre-edit source; no
unrelated executable changes were introduced in the three affected modules.

## Links

- [[redundant-and-legacy-code-map-20260623]]
- [[naming-and-claims-review]]
- [[eigenvector-calibration-construction-review]]
