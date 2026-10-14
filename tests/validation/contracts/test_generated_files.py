"""Artifact preservation and Git checks in isolated repositories."""
import json
import subprocess

import pytest
from scripts.generated_files import (
    CATALOG,
    archive,
    archive_path,
    check,
    declared_outputs,
    keep_local,
    load_catalog,
    restore,
    sha256,
)


@pytest.fixture
def repo(tmp_path):
    root = tmp_path / 'repo'
    study = root / 'reports/example'
    study.mkdir(parents=True)
    subprocess.run(['git', 'init', '-q', str(root)], check=True)
    (study / 'README.md').write_text('[Detailed rows](rows.csv)\n')
    (study / 'summary.csv').write_text('count\n3\n')
    (study / 'rows.csv').write_text('value\n1\n2\n3\n')
    data = {'schema_version': 1, 'max_staged_bytes': 1024, 'studies': {
        'reports/example': {'description': 'Example', 'reproduction_source': 'reports/example/README.md'}}, 'files': []}
    for name, role in [('README.md', 'source'), ('summary.csv', 'retained_result'), ('rows.csv', 'generated_output')]:
        p = study / name
        item = {'study': 'reports/example', 'path': p.relative_to(root).as_posix(),
                'role': role, 'reason': 'Test fixture'}
        if role != 'source':
            item.update(bytes=p.stat().st_size, sha256=sha256(p))
        data['files'].append(item)
    (root / CATALOG).write_text(json.dumps(data))
    return root, load_catalog(root), tmp_path / 'archive'


def test_archive_restore_and_local_only_missing_output(repo):
    root, data, destination = repo
    original = root / 'reports/example/rows.csv'
    before = original.read_bytes()
    assert archive(root, data, destination) == 1
    assert original.read_bytes() == before
    assert not check(root, data, verify_archive=True)
    subprocess.run(['git', '-C', str(root), 'check-ignore', '-q', str(original)], check=True)
    original.unlink()
    assert not check(root, data, verify_archive=True)
    assert restore(root, data, [], None) == 1
    assert original.read_bytes() == before


def test_corrupt_archive_never_overwritten(repo):
    root, data, destination = repo
    archive(root, data, destination)
    target = archive_path(data['files'][-1])
    target.write_text('corruption')
    with pytest.raises(ValueError, match='changed artifact'):
        archive(root, data, destination)
    assert target.read_text() == 'corruption'
    assert any('changed artifact' in e for e in check(root, data, verify_archive=True))


def test_changed_summary_fails_checksum(repo):
    root, data, destination = repo
    archive(root, data, destination)
    (root / 'reports/example/summary.csv').write_text('count\n4\n')
    assert any('summary.csv' in e for e in check(root, data))


def test_unclassified_generated_file_is_not_hidden(repo):
    root, data, destination = repo
    archive(root, data, destination)
    p = root / 'reports/example/generated/new.csv'
    p.parent.mkdir()
    p.write_text('new')
    (root / '.gitignore').write_text('reports/*/generated/\n')
    assert any('Unclassified report file' in e for e in check(root, data))


def test_staged_size_checks_index_not_working_file(repo):
    root, data, destination = repo
    archive(root, data, destination)
    big = root / 'big.bin'
    big.write_bytes(b'x' * 1025)
    subprocess.run(['git', '-C', str(root), 'add', 'big.bin'], check=True)
    big.write_bytes(b'small now')
    assert any('Oversized staged file' in e for e in check(root, data))


def test_staged_unclassified_output_is_detected(repo):
    root, data, destination = repo
    archive(root, data, destination)
    p = root / 'reports/unclassified.csv'
    p.write_text('unknown')
    subprocess.run(['git', '-C', str(root), 'add', 'reports/unclassified.csv'], check=True)
    assert any('Unclassified staged report file' in e for e in check(root, data))


def test_archive_rejects_tracked_file_and_internal_destination(repo):
    root, data, destination = repo
    with pytest.raises(ValueError, match='outside'):
        archive(root, data, root / 'archive')
    subprocess.run(['git', '-C', str(root), 'add', 'reports/example/rows.csv'], check=True)
    with pytest.raises(ValueError, match='already tracked'):
        archive(root, data, destination)
    assert not destination.exists()


def test_catalog_rejects_path_escape(repo):
    root, data, _ = repo
    data['files'][0]['path'] = '../outside'
    (root / CATALOG).write_text(json.dumps(data))
    with pytest.raises(ValueError, match='Unsafe'):
        load_catalog(root)


def test_missing_archive_is_explicit_in_strict_verification(repo):
    root, data, destination = repo
    archive(root, data, destination)
    archive_path(data['files'][-1]).unlink()
    assert not check(root, data)
    assert check(root, data, verify_archive=True)


def test_missing_reference_and_modified_ignore_rules_are_reported(repo):
    root, data, destination = repo
    archive(root, data, destination)
    (root / 'reports/example/README.md').write_text('[Lost](lost.csv)\n')
    (root / 'reports/.gitignore').write_text('')
    errors = check(root, data)
    assert any('Unresolved report reference' in e for e in errors)
    assert any('.gitignore' in e for e in errors)


def test_restore_refuses_to_replace_modified_original(repo):
    root, data, destination = repo
    archive(root, data, destination)
    original = root / 'reports/example/rows.csv'
    original.write_text('new experiment')
    with pytest.raises(ValueError, match='changed artifact'):
        restore(root, data, [], None)
    assert original.read_text() == 'new experiment'


def test_keep_local_records_no_backup_and_preserves_original(repo):
    root, data, destination = repo
    assert keep_local(root, data) == 1
    assert data['files'][-1]['storage'] == 'local'
    assert not data['files'][-1].get('archive_uri')
    assert not destination.exists()
    assert not check(root, data)


def test_retained_file_cannot_be_silently_ignored(repo):
    root, data, destination = repo
    archive(root, data, destination)
    (root / '.gitignore').write_text('*.csv\n')
    assert any('Retained file is ignored' in e for e in check(root, data))


def test_report_line_links_are_not_mistaken_for_missing_files(repo):
    root, data, destination = repo
    archive(root, data, destination)
    (root / 'reports/example/README.md').write_text('[Lines](summary.csv:1)\n')
    assert not check(root, data)


def test_new_ignored_generated_directory_is_detected(repo):
    root, data, destination = repo
    archive(root, data, destination)
    p = root / 'reports/new/generated/output.csv'
    p.parent.mkdir(parents=True)
    p.write_text('unclassified')
    (root / '.gitignore').write_text('reports/*/generated/\n')
    assert any('Unclassified generated-output directory' in e for e in check(root, data))


def test_wiki_sources_resolve_missing_catalogued_outputs_only(repo, monkeypatch):
    from scripts.wiki import lint

    root, data, destination = repo
    archive(root, data, destination)
    source = root / 'reports/example/rows.csv'
    source.unlink()
    wiki = root / 'wiki/example.md'
    wiki.parent.mkdir()
    wiki.write_text('''---
title: Example
type: source
status: reviewed
updated: 2026-09-28
sources:
  - reports/example/rows.csv
tags:
  - test
---
# Example
## Summary
## Key Points
## Evidence
## Links
''')
    monkeypatch.setattr(lint, 'ROOT', root)
    assert not lint.check_page(wiki, {'example'}, '[[example]]', declared_outputs(data))
    assert any('does not exist' in e for e in lint.check_page(wiki, {'example'}, '[[example]]'))
