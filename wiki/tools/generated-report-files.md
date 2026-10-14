---
title: Generated Report Files
type: tool
status: reviewed
updated: 2026-10-14
sources:
  - docs/generated-files.md
  - docs/generated-output-inventory.md
  - reports/artifacts.json
  - scripts/generated_files.py
  - tests/validation/contracts/test_generated_files.py
  - Makefile
tags:
  - reports
  - provenance
  - git
---

# Generated Report Files

## Summary

Keep report sources and compact results in Git; record an explicit storage
decision, byte count and SHA-256 for detailed generated outputs. Existing
tracked evidence is not automatically migrated. The initial 36 large outputs
(547 MB) are archived outside the repository with `storage: archive` and
exact ignore rules; the four run logs stay local-only.

## Usage

Run `make generated-files-check` to check classification, frozen-result hashes,
local references, exact ignore rules and the 5 MiB staged-file limit. This
command is part of `make check`; a pending storage decision fails the check.

Use `scripts/generated_files.py archive --archive-root /chosen/archive` with
the project Python after choosing a destination outside the repository. Copies
are verified before catalog entries and ignore rules are updated; originals
remain. `keep-local` explicitly records local retention without a backup.
`restore` recovers archived files without replacing different local bytes.
See `docs/generated-files.md` for complete commands and limitations.

For new studies, put bulky files under `reports/<study>/generated/` and classify
them in `reports/artifacts.json`. Ignoring a directory does not exempt its files
from classification. Preserve producing scripts and reproduction documentation.

Missing finalized generated outputs may resolve wiki citations through their
catalog entries. This identifies the recorded evidence without claiming the
bytes are available. Restore or reproduce them before numerical verification;
use `check --verify-archive` to verify archive bytes.

## Evidence

- `reports/artifacts.json` records classification and output identities.
- `docs/generated-output-inventory.md` lists the 36 large files awaiting a
  storage decision; no deletion is performed by the checker.
- `tests/validation/contracts/test_generated_files.py` checks archive integrity,
  restoration, Git-index checks, local retention and citation behavior.

## Links

- [[maintenance]]
- [[schema]]
- [[branch-distribution-variance-benchmark]]
