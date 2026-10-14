"""Run a module while retaining the relocated study's original evidence paths.

The scientific protocol pins absolute evidence names. A symlink preserves file
access but changes Path.resolve(), so this launcher aliases only the authorized
physical output directory back to its original name. The migration manifest is
SHA-256 pinned before any alias is installed; its file hashes describe cutover,
not subsequent mutable worker checkpoints and logs.
"""

import argparse
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import re
import runpy
import sys


SOURCE_ROOT = Path('/home/sji/Documents/Programming/Python/embryogenesis-simulation/outputs')
TARGET_ROOT = Path('/mnt/data/Documents/Python/embryogenesis-simulation/outputs')
_SHA256 = re.compile(r'[0-9a-fA-F]{64}\Z')


def _relative_name(value):
    if not isinstance(value, str):
        raise ValueError('Migration names must be relative strings')
    name = Path(value)
    if name.is_absolute() or '..' in name.parts or not name.parts:
        raise ValueError('Invalid relative migration name: '+value)
    return name


def validate_record(record_path, record_sha256, source_root, target_root):
    """Check the pinned cutover record and the current original-root symlink."""
    if not isinstance(record_sha256, str) or not _SHA256.fullmatch(record_sha256):
        raise ValueError('A 64-character migration record SHA-256 is required')
    data = Path(record_path).read_bytes()
    if hashlib.sha256(data).hexdigest() != record_sha256.lower():
        raise ValueError('Migration record SHA-256 differs')
    record = json.loads(data)
    if record.get('passed') is not True:
        raise ValueError('Migration record has not passed verification')
    if (record.get('source_root') != str(source_root) or
            record.get('target_root') != str(target_root)):
        raise ValueError('Migration record does not name the authorized roots')
    files = record.get('files')
    if not isinstance(files, dict) or not files:
        raise ValueError('Migration record requires its verified file registry')
    for name, entry in files.items():
        _relative_name(name)
        if (not isinstance(entry, dict) or
                not isinstance(entry.get('sha256'), str) or
                not _SHA256.fullmatch(entry['sha256']) or
                not isinstance(entry.get('source_stat'), dict) or
                not isinstance(entry.get('destination_stat'), dict)):
            raise ValueError('Incomplete migration file evidence: '+name)
    symlinks = record.get('symlinks', {})
    if not isinstance(symlinks, dict):
        raise ValueError('Invalid migration symlink registry')
    for name, target in symlinks.items():
        _relative_name(name)
        if isinstance(target, dict):
            if (not isinstance(target.get('source_stat'), dict) or
                    not isinstance(target.get('destination_stat'), dict)):
                raise ValueError('Incomplete migration symlink evidence: '+name)
            target = target.get('target')
        if not isinstance(target, str):
            raise ValueError('Invalid migration symlink target: '+name)
    if not source_root.is_symlink() or not target_root.is_dir():
        raise ValueError('Original outputs must be a symlink to the destination directory')
    physical_target = target_root.resolve(strict=True)
    if source_root.resolve(strict=True) != physical_target:
        raise ValueError('Original outputs symlink does not point to the authorized destination')
    return record, physical_target


@contextmanager
def validated_alias(record_path, record_sha256, *, source_root=SOURCE_ROOT, target_root=TARGET_ROOT):
    """Install the alias temporarily after validating the pinned migration."""
    source_root, target_root = Path(source_root), Path(target_root)
    record, physical_target = validate_record(record_path, record_sha256, source_root, target_root)
    native_resolve = Path.resolve

    def resolve(path, strict=False):
        resolved = native_resolve(path, strict=strict)
        try:
            relative = resolved.relative_to(physical_target)
        except ValueError:
            return resolved
        return source_root / relative

    Path.resolve = resolve
    try:
        yield record
    finally:
        Path.resolve = native_resolve


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record', required=True, type=Path)
    parser.add_argument('--record-sha256', required=True)
    parser.add_argument('-m', dest='module', required=True)
    tokens = sys.argv[1:] if argv is None else list(argv)
    # Everything after the module name belongs to that module, including an
    # initial option or another -m. argparse REMAINDER alone rejects some such
    # initial options, so split at the wrapper's -m before parsing its flags.
    module_index = tokens.index('-m') if '-m' in tokens else len(tokens)
    args = parser.parse_args(tokens[:module_index+2])
    module_args = tokens[module_index+2:]
    old_argv, old_path = sys.argv, sys.path[:]
    with validated_alias(args.record, args.record_sha256,
                         source_root=SOURCE_ROOT, target_root=TARGET_ROOT):
        try:
            sys.argv = [args.module, *module_args]
            # Match python -m: imports start in the current working directory,
            # and runpy supplies the module's __main__, __package__ and argv[0].
            sys.path[0] = os.getcwd()
            runpy.run_module(args.module, run_name='__main__', alter_sys=True)
        finally:
            sys.argv = old_argv
            sys.path[:] = old_path


if __name__ == '__main__':
    main()
