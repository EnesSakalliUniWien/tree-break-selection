# Package build verification — 2026-11-03

The user authorized building the source distribution and wheel, inspecting the
self-distance correction and metadata, and testing an installed wheel outside
the checkout in a temporary environment.

## Build and artifacts

Executed from the repository with the existing Python 3.12.11 environment and
`flit_core` 3.12.0:

```bash
uv build --python .venv/bin/python --no-build-isolation --offline --out-dir dist/verification-20261103-gmv1fe71
```

The source distribution and wheel built successfully; uv built the wheel from
the source distribution. Artifacts remain in the ignored local directory
`dist/verification-20261103-gmv1fe71/`:

- `tree_break_selection-0.1.0-py3-none-any.whl`, 188398 bytes,
  SHA256 `7f4ba0d9ef7bc18b32cc6bd7c55d7e6368e9a411a3a58816451dad6b3e9ef54c`.
- `tree_break_selection-0.1.0.tar.gz`, 138355 bytes,
  SHA256 `da63a3172712b6e8f70b35b3c53384f23073b4cc5d36d52d4f8ffecc793a1cda`.

All 109 production Python files occur in the wheel. All 114 package files match
the checkout and both artifacts byte-for-byte. The constructor contains the
identity-diagonal reset before epsilon and the `self_distance_policy` metadata.
Its SHA256 is `81a51a5bec93fda2891b871b7a03b4a0dfc140e054eb3a61e4430ddad6f62565`.
The wheel contains the library, not the repository's benchmark/application modules.

## Installed-wheel verification

A fresh temporary environment used the same Python version. Dependencies were
constrained by a frozen export of `uv.lock`, with pytest 8.4.2 for verification.
The offline installation could not find some dependency artifacts in cache;
registry access then completed installation of 26 packages. `uv pip check`
reported all installed packages compatible.

The unchanged self-distance regression file was copied outside the repository.
Using the temporary interpreter with `-I`, pytest's `--import-mode=importlib`,
and the temporary root/confcutdir, all three direct, feature-family and block
regressions passed in 35.00 seconds. The checkout baseline also passed all three.
Separate assertions verified that diffusion imported from temporary site-packages,
the checkout was absent from `sys.path`, and installed source matched the wheel.

The temporary environment was removed after verification. Artifact inspection,
installed versions and verification JSON remain beside the archives.
`pyproject.toml` and `uv.lock` hashes stayed unchanged. No source edits,
full-suite rerun, publication, commit or push. This verifies packaging and the
correction, not selected-root calibration.
