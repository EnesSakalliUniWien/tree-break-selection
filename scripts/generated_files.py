#!/usr/bin/env python3
"""Check, archive and restore explicitly classified report files without deleting originals."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path, PurePosixPath
from urllib.parse import unquote, urlparse

ROOT = Path(__file__).resolve().parents[1]
CATALOG = 'reports/artifacts.json'
ROLES = {'source', 'retained_result', 'generated_output'}
BEGIN = '# BEGIN catalogued generated outputs'
END = '# END catalogued generated outputs'


def git(root: Path, *args: str) -> bytes:
    return subprocess.check_output(['git', '-C', str(root), *args])


def sha256(path: Path) -> str:
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def safe_path(root: Path, name: str) -> Path:
    parts = PurePosixPath(name)
    if parts.is_absolute() or '..' in parts.parts or '\\' in name or not name:
        raise ValueError(f'Unsafe repository path: {name!r}')
    path = root / name
    if not path.resolve().is_relative_to(root.resolve()):
        raise ValueError(f'Path escapes repository: {name}')
    if path.is_symlink():
        raise ValueError(f'Symlink is not an artifact file: {name}')
    return path


def load_catalog(root: Path) -> dict:
    data = json.loads((root / CATALOG).read_text())
    if data.get('schema_version') != 1:
        raise ValueError('Unsupported artifact catalog schema')
    if not isinstance(data.get('max_staged_bytes'), int) or data['max_staged_bytes'] <= 0:
        raise ValueError('max_staged_bytes must be positive')
    studies = data['studies']
    for study, record in studies.items():
        if len(PurePosixPath(study).parts) != 2 or not study.startswith('reports/'):
            raise ValueError(f'Expected reports/<study>: {study}')
        safe_path(root, study)
        if not record.get('reproduction_source') or not record.get('description'):
            raise ValueError(f'Missing study provenance: {study}')
        safe_path(root, record['reproduction_source'])
    seen = set()
    for item in data['files']:
        name = item['path']
        safe_path(root, name)
        if name in seen or item['study'] not in studies or not name.startswith(item['study'] + '/'):
            raise ValueError(f'Duplicate or unowned artifact: {name}')
        seen.add(name)
        if item['role'] not in ROLES or not item.get('reason'):
            raise ValueError(f'Missing classification: {name}')
        if item['role'] != 'source':
            if not re.fullmatch('[0-9a-f]{64}', item.get('sha256', '')):
                raise ValueError(f'Missing SHA-256: {name}')
            if not isinstance(item.get('bytes'), int) or item['bytes'] < 0:
                raise ValueError(f'Missing byte count: {name}')
        uri = item.get('archive_uri')
        if item['role'] == 'generated_output' and item.get('storage', 'pending') not in {'pending', 'local', 'archive'}:
            raise ValueError(f'Unknown output storage: {name}')
        if uri is not None and (item['role'] != 'generated_output' or not uri.startswith('file:///')):
            raise ValueError(f'Unsupported archive URI: {name}')
        if item.get('storage') == 'archive' and not uri:
            raise ValueError(f'Missing archive URI: {name}')
    return data


def verify_file(path: Path, item: dict) -> None:
    if not path.is_file() or path.stat().st_size != item['bytes'] or sha256(path) != item['sha256']:
        raise ValueError(f'Missing or changed artifact: {path}')


def archive_path(item: dict, archive_root: Path | None = None) -> Path:
    if archive_root is not None:
        return archive_root / 'sha256' / item['sha256'][:2] / item['sha256']
    uri = item.get('archive_uri')
    if not uri:
        raise ValueError(f'Archive destination pending: {item["path"]}')
    parsed = urlparse(uri)
    if parsed.scheme != 'file' or parsed.netloc:
        raise ValueError(f'Unsupported local archive URI: {uri}')
    return Path(unquote(parsed.path))


def ignore_block(data: dict) -> str:
    names = sorted(i['path'].removeprefix('reports/') for i in data['files']
                   if i['role'] == 'generated_output' and i.get('storage') in {'local', 'archive'})
    # Escape Git pattern metacharacters so each rule identifies one literal file.
    escaped = ['/' + re.sub(r'([\\*?\[\]])', r'\\\1', name) for name in names]
    return '\n'.join([BEGIN, *escaped, END]) + '\n'


def atomic_text(path: Path, text: str) -> None:
    with tempfile.NamedTemporaryFile('w', dir=path.parent, delete=False) as stream:
        temporary = Path(stream.name)
        stream.write(text)
    os.replace(temporary, path)


def copy_verified(source: Path, destination: Path, item: dict) -> None:
    verify_file(source, item)
    if destination.exists():
        verify_file(destination, item)
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=destination.parent, delete=False) as stream:
        temporary = Path(stream.name)
    try:
        shutil.copyfile(source, temporary)
        verify_file(temporary, item)
        # An exclusive hard link prevents replacing an existing archive object.
        try:
            os.link(temporary, destination)
        except FileExistsError:
            verify_file(destination, item)
    finally:
        temporary.unlink(missing_ok=True)


def archive(root: Path, data: dict, destination: Path) -> int:
    destination = destination.expanduser().resolve()
    if destination.is_relative_to(root.resolve()):
        raise ValueError('Archive destination must be outside the repository')
    tracked = set(git(root, 'ls-files', '-z').decode().split('\0'))
    items = [i for i in data['files'] if i['role'] == 'generated_output']
    for item in items:
        if item['path'] in tracked:
            raise ValueError(f'Refusing to ignore an already tracked file: {item["path"]}')
        source = safe_path(root, item['path'])
        if not source.exists() and item.get('archive_uri'):
            source = archive_path(item)
        copy_verified(source, archive_path(item, destination), item)
    # Publish archive locations and ignore rules only after every copy verifies.
    for item in items:
        item['archive_uri'] = archive_path(item, destination).as_uri()
        item['storage'] = 'archive'
    save_catalog(root, data)
    return len(items)


def save_catalog(root: Path, data: dict) -> None:
    atomic_text(root / CATALOG, json.dumps(data, indent=2) + '\n')
    ignore = root / 'reports/.gitignore'
    existing = ignore.read_text() if ignore.exists() else ''
    if BEGIN in existing:
        existing = re.sub(re.escape(BEGIN) + r'.*?' + re.escape(END) + r'\n?', '', existing, flags=re.S)
    atomic_text(ignore, existing.rstrip() + '\n\n' + ignore_block(data))


def keep_local(root: Path, data: dict) -> int:
    """Explicitly retain outputs locally without claiming they are backed up."""
    tracked = set(git(root, 'ls-files', '-z').decode().split('\0'))
    items = [i for i in data['files'] if i['role'] == 'generated_output' and i.get('storage', 'pending') == 'pending']
    for item in items:
        if item['path'] in tracked:
            raise ValueError(f'Refusing to ignore an already tracked file: {item["path"]}')
        verify_file(safe_path(root, item['path']), item)
    for item in items:
        item['storage'] = 'local'
    save_catalog(root, data)
    return len(items)


def restore(root: Path, data: dict, names: list[str], archive_root: Path | None) -> int:
    items = {i['path']: i for i in data['files'] if i['role'] == 'generated_output'}
    selected = names or list(items)
    for name in selected:
        if name not in items:
            raise ValueError(f'Not a catalogued archived output: {name}')
        item = items[name]
        copy_verified(archive_path(item, archive_root), safe_path(root, name), item)
    return len(selected)


def declared_outputs(data: dict) -> set[str]:
    return {i['path'] for i in data['files'] if i['role'] == 'generated_output'
            and i.get('storage') in {'local', 'archive'}}


def check(root: Path, data: dict, *, verify_archive: bool = False, archive_root: Path | None = None) -> list[str]:
    errors = []
    items = {i['path']: i for i in data['files']}
    tracked = set(git(root, 'ls-files', '-z').decode().split('\0'))
    for name, item in items.items():
        path = safe_path(root, name)
        try:
            if item['role'] == 'source':
                if not path.is_file():
                    raise ValueError(f'Missing source: {name}')
            elif item['role'] == 'retained_result' or path.exists():
                verify_file(path, item)
            if item['role'] == 'generated_output':
                if name in tracked:
                    raise ValueError(f'Archived output is tracked: {name}')
                if item.get('storage', 'pending') == 'pending':
                    raise ValueError(f'Storage decision pending: {name}')
                if verify_archive and item.get('storage') == 'archive':
                    verify_file(archive_path(item, archive_root), item)
        except ValueError as exc:
            errors.append(str(exc))
    for study, record in data['studies'].items():
        if not safe_path(root, record['reproduction_source']).is_file():
            errors.append(f'Missing reproduction documentation: {study}')
        for path in (root / study).rglob('*'):
            if path.is_file() and '__pycache__' not in path.parts and path.name != '.DS_Store':
                name = path.relative_to(root).as_posix()
                if name not in items:
                    errors.append(f'Unclassified report file: {name}')
    controls = {'reports/README.md', 'reports/.gitignore', CATALOG}
    visible = git(root, 'ls-files', '--others', '--exclude-standard', '-z', 'reports').decode().split('\0')
    for name in visible:
        if name and name not in controls and name not in items:
            errors.append(f'Unclassified new report file: {name}')
    for directory in (root / 'reports').glob('*/generated'):
        study = directory.parent.relative_to(root).as_posix()
        if study in data['studies'] or study == 'reports/audits':
            continue  # Managed studies checked above; old audit output is pre-existing.
        if any(p.is_file() for p in directory.rglob('*')):
            errors.append(f'Unclassified generated-output directory: {study}/generated')
    staged = git(root, 'diff', '--cached', '--name-only', '--diff-filter=ACMR', '-z').decode().split('\0')
    new_staged = set(git(root, 'diff', '--cached', '--name-only', '--diff-filter=A', '-z').decode().split('\0'))
    for name in filter(None, staged):
        size = int(git(root, 'cat-file', '-s', ':' + name))
        if size > data['max_staged_bytes']:
            errors.append(f'Oversized staged file ({size} bytes): {name}')
        if name.startswith('reports/') and name not in controls and name not in items and name in new_staged:
            errors.append(f'Unclassified staged report file: {name}')
    ignore = root / 'reports/.gitignore'
    if any(i.get('storage') in {'local', 'archive'} for i in items.values()):
        if not ignore.exists() or ignore_block(data) not in ignore.read_text():
            errors.append('reports/.gitignore does not match catalogued archived outputs')
    ignored_result = subprocess.run(
        ['git', '-C', str(root), 'check-ignore', '--no-index', '-z', '--stdin'],
        input=('\0'.join(items) + '\0').encode(), capture_output=True, check=False,
    )
    if ignored_result.returncode not in (0, 1):
        raise ValueError('git check-ignore failed: ' + ignored_result.stderr.decode())
    ignored = set(ignored_result.stdout.decode().split('\0'))
    for name, item in items.items():
        if item['role'] != 'generated_output' and name in ignored:
            errors.append(f'Retained file is ignored: {name}')
        if name in declared_outputs(data) and name not in ignored:
            errors.append(f'Generated output is not ignored: {name}')
    # Local Markdown report links may resolve through the archive catalog on a fresh clone.
    archived = declared_outputs(data)
    for name, item in items.items():
        path = root / name
        if item['role'] != 'source' or path.suffix != '.md' or not path.exists():
            continue
        text = re.sub(r'```.*?```', '', path.read_text(), flags=re.S)
        for link in re.findall(r'\]\(([^)]+)\)', text):
            link = unquote(link.split('#')[0].strip('<>'))
            link = re.sub(r':\d+(?:-\d+)?$', '', link)
            if not link or re.match(r'^[a-zA-Z][a-zA-Z0-9+.-]*:', link):
                continue
            target = (path.parent / link).resolve()
            if target.is_relative_to(root.resolve()) and not target.exists():
                rel = target.relative_to(root.resolve()).as_posix()
                if rel not in archived:
                    errors.append(f'Unresolved report reference in {name}: {link}')
    return sorted(set(errors))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    sub = parser.add_subparsers(dest='command', required=True)
    audit = sub.add_parser('check')
    audit.add_argument('--verify-archive', action='store_true')
    audit.add_argument('--archive-root', type=Path)
    save = sub.add_parser('archive')
    save.add_argument('--archive-root', type=Path, required=True)
    sub.add_parser('keep-local')
    recover = sub.add_parser('restore')
    recover.add_argument('paths', nargs='*')
    recover.add_argument('--archive-root', type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    try:
        data = load_catalog(root)
        if args.command == 'check':
            errors = check(root, data, verify_archive=args.verify_archive, archive_root=args.archive_root)
            if errors:
                print('\n'.join(errors))
                return 1
            local = sum(i.get('storage') == 'local' for i in data['files'])
            missing = sum(not (root / n).exists() for n in declared_outputs(data))
            print(f'generated-files check passed: {len(data["files"])} classified files; '
                  f'{local} local-only outputs; {missing} outputs absent locally; '
                  f'archive bytes {"verified" if args.verify_archive else "not checked"}')
        elif args.command == 'archive':
            print(f'Archived and verified {archive(root, data, args.archive_root)} files; originals preserved')
        elif args.command == 'keep-local':
            print(f'Keeping {keep_local(root, data)} outputs locally and ignoring them in Git; no backup created')
        else:
            print(f'Restored or verified {restore(root, data, args.paths, args.archive_root)} files')
    except (ValueError, KeyError, OSError, subprocess.CalledProcessError) as exc:
        print(f'generated-files failed: {exc}')
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
