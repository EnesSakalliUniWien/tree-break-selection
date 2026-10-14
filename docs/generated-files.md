# Generated files and Git

Keep scripts, protocols, mathematical explanations, configuration, seeds,
source/environment hashes, compact results, verification and selected figures
in Git. Full numerical outputs have an explicit storage decision in
`reports/artifacts.json`. Size alone is not proof that a file is disposable.
The initial reviewed selection is in `docs/generated-output-inventory.md`.

## Layout and catalog

For new runs use `reports/<study>/generated/` for bulky outputs. Keep the
reviewable report and summaries outside that directory. Existing studies keep
their paths so reproduction commands and historical evidence are unchanged.
The catalog lists each file, its study, role and reason:

- `source`: editable code, protocol or report; existence is checked.
- `retained_result`: immutable numerical summary, provenance, verification,
  historical source snapshot or figure; exact size and SHA-256 are checked.
- `generated_output`: detailed output with exact size and SHA-256. Storage is
  `pending`, `local` (no separate backup), or `archive` (recorded archive URI).

Each study points to its existing reproduction documentation. Original run
manifests remain unchanged; the catalog supplements them. Intentional changes
to a frozen result require reviewing its producing run and updating its catalog
entry explicitly. The checker never silently updates checksums.

`reports/.gitignore` has exact paths for finalized generated outputs. It does
not ignore whole historical report folders or file extensions. Existing tracked
evidence stays tracked. No command here stages, commits, untracks or deletes
repository files. Existing local logs remain explicitly local-only.

## Commands

```sh
make generated-files-check
```

`make check` runs it; a new output with `pending` storage fails it. It
detects unclassified report files, changed retained results, missing source files, unresolved local Markdown
references, incorrect ignore rules, force-added generated outputs, and staged
blobs larger than 5 MiB. The staged-size check reads Git's index, not the
potentially different working file. Already tracked legacy evidence is not
automatically migrated; modifying/staging an oversized legacy file still triggers
the size check. The 5 MiB value is a review guard, not a Git service limit.

Archive to a user-chosen directory **outside** the repository:

```sh
.venv/bin/python scripts/generated_files.py archive --archive-root /chosen/archive
.venv/bin/python scripts/generated_files.py check --verify-archive
```

The first command verifies original checksums, copies into
`sha256/<first-two-characters>/<full-hash>`, verifies each copy, then records
`file://` locations and exact ignore rules. Identical bytes share one object.
Existing archive objects are never overwritten. Original outputs remain in
place. A local archive on the same disk is not an off-machine backup.

If the explicit decision is to keep files locally without a separate backup:

```sh
.venv/bin/python scripts/generated_files.py keep-local
```

This finalizes only pending files, verifies them and writes exact ignore rules.
It does not claim they are archived, and it does not remove them.

Restore archived files to their original report paths:

```sh
.venv/bin/python scripts/generated_files.py restore reports/example/generated/rows.csv
```

Omitting file paths selects all generated outputs; local-only files without an
archive cannot be restored by this command. Use `--archive-root` when the same
content-addressed archive is mounted at a different location. Restore refuses
to overwrite a different local file. Use the study's reproduction documentation
when no archive exists.

## Fresh clones and evidence availability

The ordinary check validates locally present bytes. It does not require a
mounted archive or full benchmark outputs on every clone, and it reports how
many catalogued outputs are absent. `--verify-archive` actually reads archive
bytes and fails on missing/corrupt archive objects. Neither mode validates
numerical correctness; that remains the benchmark's own verification job.

Report references and wiki source paths may resolve through a finalized
catalog entry when full output bytes are absent. This preserves the recorded
evidence identity; it does not imply the data were downloaded or inspected.
Fetch or reproduce the outputs before running analyses that consume them.

New report files require explicit catalog entries, including files in ignored
`generated/` directories. Existing `reports/audits/generated/` remains under
its older local-output rule. Caches such as `__pycache__/` are not research
outputs. Historical tracked `raw/assets/` evidence is outside this migration.
