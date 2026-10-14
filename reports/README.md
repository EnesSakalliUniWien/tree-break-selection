# Reports

Track reviewable reports, scripts, protocols, compact results and provenance.
Classify detailed generated outputs in `artifacts.json` before archiving or
ignoring them. See [the generated-file policy](../docs/generated-files.md).

Use this directory only for retained evidence that is intentionally cited by the
manuscript, wiki, or a review note. Routine logs, profiling CSVs, coverage
files, notebook image exports, and duplicate-code scans should stay local or be
regenerated from the source command. New bulky outputs belong under
`reports/<study>/generated/`; existing cited paths remain unchanged.
