"""CPU-only checks for the external output relocation launcher."""

import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import pytest


SCRIPT = Path(__file__).parent.parent / 'scripts' / 'with_output_path_alias.py'
SPEC = importlib.util.spec_from_file_location('output_path_alias_launcher', SCRIPT)
launcher = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(launcher)


def migration(tmp_path):
    source, target = tmp_path/'original'/'outputs', tmp_path/'destination'/'outputs'
    source.parent.mkdir(); target.mkdir(parents=True)
    source.symlink_to(target, target_is_directory=True)
    (target/'nested').mkdir(); (target/'nested'/'state').write_bytes(b'checkpoint')
    (target/'state-link').symlink_to('nested/state')
    manifest = dict(passed=True, source_root=str(source), target_root=str(target),
        files={'nested/state': dict(sha256=hashlib.sha256(b'checkpoint').hexdigest(),
            source_stat={'size': 10}, destination_stat={'size': 10})},
        symlinks={'state-link': dict(target='nested/state', source_stat={'mode': 41471},
            destination_stat={'mode': 41471})})
    record = tmp_path/'manifest.json'
    record.write_text(json.dumps(manifest))
    sha = hashlib.sha256(record.read_bytes()).hexdigest()
    return source, target, record, sha


def test_alias_preserves_original_paths_and_leaves_unrelated_paths_unchanged(tmp_path):
    source, target, record, sha = migration(tmp_path)
    outside = tmp_path/'outside'; outside.mkdir()
    sibling = target.parent/'outputs-other'; sibling.mkdir()
    native = Path.resolve
    with launcher.validated_alias(record, sha, source_root=source, target_root=target):
        assert source.resolve() == source
        assert target.resolve() == source
        assert (source/'nested'/'state').resolve(strict=True) == source/'nested'/'state'
        assert (target/'nested'/'state').resolve() == source/'nested'/'state'
        assert outside.resolve() == native(outside)
        assert sibling.resolve() == native(sibling)
    assert Path.resolve is native
    assert source.resolve() == target


def test_native_strict_behavior_and_symlink_destinations_are_preserved(tmp_path):
    source, target, record, sha = migration(tmp_path)
    outside = tmp_path/'outside'; outside.mkdir()
    (outside/'normal').write_text('outside')
    (target/'inside-link').symlink_to(target/'nested')
    (target/'escape-link').symlink_to(outside)
    (outside/'normal-link').symlink_to(outside/'normal')
    (outside/'output-link').symlink_to(target/'nested'/'state')
    with launcher.validated_alias(record, sha, source_root=source, target_root=target):
        assert (source/'inside-link'/'state').resolve(strict=True) == source/'nested'/'state'
        assert (source/'escape-link'/'normal').resolve(strict=True) == outside/'normal'
        assert (outside/'normal-link').resolve(strict=True) == outside/'normal'
        assert (outside/'output-link').resolve(strict=True) == source/'nested'/'state'
        assert (source/'missing').resolve(strict=False) == source/'missing'
        with pytest.raises(FileNotFoundError):
            (source/'missing').resolve(strict=True)
        with pytest.raises(FileNotFoundError):
            (source/'escape-link'/'missing').resolve(strict=True)


@pytest.mark.parametrize('change', ['sha', 'passed', 'root', 'symlink', 'registry'])
def test_alias_is_never_installed_with_unverified_migration_evidence(tmp_path, change):
    source, target, record, sha = migration(tmp_path)
    if change == 'sha':
        sha = '0'*64
    elif change == 'symlink':
        source.unlink(); source.symlink_to(tmp_path)
    else:
        manifest = json.loads(record.read_text())
        if change == 'passed':
            manifest['passed'] = False
        elif change == 'root':
            manifest['target_root'] = str(tmp_path/'wrong')
        else:
            manifest['files']['../escape'] = manifest['files'].pop('nested/state')
        record.write_text(json.dumps(manifest)); sha = hashlib.sha256(record.read_bytes()).hexdigest()
    native = Path.resolve
    with pytest.raises(ValueError):
        with launcher.validated_alias(record, sha, source_root=source, target_root=target):
            pytest.fail('Unverified alias installed')
    assert Path.resolve is native


def test_alias_is_restored_after_module_exception(tmp_path):
    source, target, record, sha = migration(tmp_path)
    native = Path.resolve
    with pytest.raises(RuntimeError):
        with launcher.validated_alias(record, sha, source_root=source, target_root=target):
            raise RuntimeError('module failed')
    assert Path.resolve is native


@pytest.mark.parametrize('starts_with_option', [False, True])
def test_launcher_keeps_module_argv_and_package_semantics(tmp_path, monkeypatch, starts_with_option):
    source, target, record, sha = migration(tmp_path)
    package_name = f'relocation_launcher_test_package_{int(starts_with_option)}'
    package = tmp_path/package_name; package.mkdir()
    (package/'__init__.py').write_text('')
    result = tmp_path/'module-result.json'
    (package/'entry.py').write_text(
        'import json, sys\nfrom pathlib import Path\n'
        f'Path({str(result)!r}).write_text(json.dumps(dict(argv=sys.argv, name=__name__, '
        'package=__package__, spec=__spec__.name, output=str(Path(sys.argv[-1]).resolve()))))\n')
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(launcher, 'SOURCE_ROOT', source)
    monkeypatch.setattr(launcher, 'TARGET_ROOT', target)
    old_argv, old_path, native = sys.argv, sys.path[:], Path.resolve
    module_args = (['--output', str(source)] if starts_with_option else
                   [str(result), '--output', str(source)])
    launcher.main(['--record', str(record), '--record-sha256', sha, '-m',
        package_name+'.entry', *module_args])
    report = json.loads(result.read_text())
    assert report['argv'][0] == str(package/'entry.py')
    assert report['argv'][1:] == module_args
    assert report['name'] == '__main__'
    assert report['package'] == package_name
    assert report['spec'] == package_name+'.entry'
    assert report['output'] == str(source)
    assert sys.argv is old_argv and sys.path == old_path and Path.resolve is native
