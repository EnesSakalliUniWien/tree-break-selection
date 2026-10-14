# Dead code, aliases and duplication audit

**Follow-up:** the user subsequently authorized removal of the three assignment
aliases. All callers now use `TREE_CONSENSUS_STRATEGIES`,
`TAIL_LAW_CONTEXT_COLUMNS` and `METHOD_CONTRAST_COLUMNS` directly. The audit
below and its saved outputs describe the pre-removal snapshot. The private
helper and duplicate implementations were outside this removal request.

One Sol agent, using the cavecrew-investigator workflow, first identified the
repository's existing analysis tools and verified their availability. It then
ran the configured checks and a bounded supplemental static analysis. The
coordinating agent checked the retained candidates against source and usages.

Result: one high-confidence repository-unused private helper and two groups
of live duplicate implementations are cleanup candidates. Three inspected
assignment aliases have active uses. No runtime code was changed or removed;
this report and its wiki records preserve the findings for a separate cleanup.

## Tools found before analysis

The Makefile and pyproject configure Vulture, deptry and Ruff; executable
versions were checked against uv.lock. No current clone checker was found in
the manifest, lockfile or scripts. Historical jscpd/Pylint reports do not make
those tools part of the current configured workflow.

| Tool/check | Observed result |
|---|---|
| `uv run --no-sync deptry applications tree_break_selection benchmarks` | Exit 0; 480 files; no dependency issues. |
| `uv run --no-sync vulture applications benchmarks scripts tree_break_selection --min-confidence 90` | Exit 0; no findings at the repository threshold. |
| `uv run --no-sync ruff check --no-cache .` | Exit 0; configured E/F/I checks pass. |
| Supplemental Vulture, adding tests and using `--min-confidence 60` | Exit 3; 62 candidate lines, including false positives. This is detector output, not 62 confirmed dead items. |
| Supplemental Python AST plus `rg` | 623 Python files; 13 exact function-body groups of at least eight source lines; one unreferenced private-function candidate; three top-level assignment aliases. |

All `uv` checks used `UV_OFFLINE=true`; no dependency installation or network
resolution was performed. Versions: Python 3.11.12, uv 0.9.5, deptry 0.25.1,
Vulture 2.16, Ruff 0.15.14, rg 15.2.0. Exact configured commands and exit codes
are in [audit_commands.json](audit_commands.json).

The AST comparison strips leading docstrings and compares function bodies
without renaming identifiers. It does not require matching signatures; matches
therefore need manual contract checks before consolidation. This supplemental
analysis is not an installed clone checker or a whole-program reachability
proof. It is not comparable to the earlier weak-mode clone census.

## Cleanup candidates

### D1 — Repository-unused private helper

[`_node_pvalue`:550](../../benchmarks/diagnostics/calibration/sibling/known_branch_skew_diagnostic.py#L550)
has no observed calls, exports, decorators or dynamic lookup. The scanned code
contains no name, attribute, imported-name or exact-string reference to it.
The module's CLI and
[focused tests](../../tests/validation/calibration/195_test_known_branch_skew_diagnostic.py)
reach the active calculations through `run_diagnostic` and `_local_test`.
The helper also assigns unused `left_name` and `right_name` variables.

Assessment: high-confidence unused code within the audited repository, not a
claim about unknown external consumers. Candidate action: remove this private
helper and run the known-branch diagnostic tests. No removal was performed.

### D2 — Duplicate Git-provenance implementation

`current_git_state` has identical function bodies in
[alpha_grid_search.py:234](../../benchmarks/validation/sweeps/alpha_grid_search.py#L234)
and
[traversal_sibling_fdr_null.py:514](../../benchmarks/validation/statistics/traversal_sibling_fdr_null.py#L514).
Both are live manifest producers, including build-environment fields and
explicit subprocess failure records.

Candidate action: consolidate the implementation while preserving the exact
output and error contract. Removing one implementation without updating its
callers would break a live path.

### D3 — Duplicate timestamp rendering

The same `fig.text` body is used in four live helpers:

- [selected_nnls_report.py:189](../../applications/scrna/plots/selected_nnls_report.py#L189), `add_generated_at`.
- [selected_nnls_fit_summary.py:198](../../applications/scrna/plots/selected_nnls_fit_summary.py#L198), `_add_generated_at`.
- [analyze_goncalves_progenitors.py:279](../../applications/scrna/analysis/analyze_goncalves_progenitors.py#L279), timestamp helper.
- [MNIST plot_report.py:34](../../applications/mnist/plot_report.py#L34), timestamp helper.

Candidate action: share the renderer while preserving placement, style and
timestamp text. These functions are not dead. Similar-looking plots with
different contracts were not automatically grouped as cleanup targets.

## Aliases and dynamic APIs to retain

| Symbol | Evidence and assessment |
|---|---|
| `TOPOLOGIES = TREE_CONSENSUS_STRATEGIES` in [family_metric_nnls_grid.py:62](../../benchmarks/validation/sweeps/family_metric_nnls_grid.py#L62) | Used by defaults, CLI and validation; equality with the shared constant is tested. Active alias. |
| `BASE_CONTEXT_COLUMNS = TAIL_LAW_CONTEXT_COLUMNS` in [selected_tail_topology_refinement.py:36](../../benchmarks/diagnostics/calibration/selected/tail/selected_tail_topology_refinement.py#L36) | Shared grouping/default/manifest contract. Active alias. |
| `CANDIDATE_METHOD_CONTRAST_COLUMNS = METHOD_CONTRAST_COLUMNS` in [selected_neighborhood_distribution_panel.py:306](../../benchmarks/diagnostics/calibration/selected/neighborhood/selected_neighborhood_distribution_panel.py#L306) | Both report builders intentionally share the same output columns. Active alias. |
| `_import_runner` in [method_registry.py:38](../../benchmarks/shared/runners/method_registry.py#L38) | Dynamically resolves six runners flagged by the lower-threshold scan. Preserve registered runners. |
| `bootstrap_consensus.py:81` (removed 2026-10-13) | Documented optional API with focused tests. No application consumer found does not mean dead code. |

No removable alias was established among the three inspected assignment
aliases. This was not an exhaustive census of import aliases, wrappers or
method IDs.

## Scope and evidence

The supplemental scan covers non-ignored Python under applications,
benchmarks, scripts, tree_break_selection and tests, including dirty/untracked
work. Vendor, manuscript, saved reports/raw provenance and ignored analysis
workspace are outside its candidate census. Mirrored test fixtures are not
reported as production duplication. Dynamic runtime graphs and external
consumers were not executed or proven absent. Remaining lower-confidence
candidates are unresolved, not approved for deletion.

[AST candidates](ast_candidates.json), [AST summary](ast-summary.txt),
[deptry output](deptry.txt), [Ruff output](ruff.txt),
[configured Vulture output](vulture-90.txt), and
[supplemental Vulture candidates](vulture-60.txt) preserve the audit evidence.
The custom AST algorithm is described in the command record; its transient
heredoc is not packaged as a reusable repository tool.

Wiki lint passed for 209 pages before the audit record was added. The new
report/wiki links and post-edit wiki lint are checked separately. No unit-test
or benchmark rerun is claimed for this source review; no code, alias or
dependency was removed, and nothing was committed or pushed.

## Follow-up verification: removal of the three aliases

Removed `TOPOLOGIES`, `BASE_CONTEXT_COLUMNS` and
`CANDIDATE_METHOD_CONTRAST_COLUMNS` from their defining benchmark modules and
replaced their executable usages with the original constants. Updated the
existing tree-default regression test to check the CLI's topology defaults.
Parameter names, CLI flags, output column names and constant values remain
unchanged.

The three focused tests pass before and after removal; changed-file Ruff
passes. Import checks confirm the old attributes are absent, and verify CLI
defaults, refinement defaults and candidate-report columns. AST comparison
against the immediate pre-edit snapshot confirms that the three implementation
files differ only by removing the assignments and substituting the original
constant names. Exact-name search finds no remaining executable references in
package, applications, benchmarks, scripts or tests. Historical audit evidence
retains the original names. No commit or push was performed.
