"""Verify the new storage root uses the existing guarded logical-path alias."""
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

spec=importlib.util.spec_from_file_location('singular_output_alias',
    Path(__file__).parent.parent/'scripts/with_study_output_path_alias.py')
wrapper=importlib.util.module_from_spec(spec);spec.loader.exec_module(wrapper)


def test_singular_output_launcher_preserves_names_and_restores_alias(tmp_path, monkeypatch):
    source=tmp_path/'original/output';target=tmp_path/'disk/output'
    source.parent.mkdir();target.mkdir(parents=True);source.symlink_to(target)
    (target/'state').write_bytes(b'state')
    record=tmp_path/'manifest.json'
    record.write_text(json.dumps(dict(passed=True,source_root=str(source),target_root=str(target),
        files={'state':dict(sha256=hashlib.sha256(b'state').hexdigest(),source_stat={},destination_stat={})})))
    digest=hashlib.sha256(record.read_bytes()).hexdigest()
    monkeypatch.setattr(wrapper,'SOURCE_ROOT',source);monkeypatch.setattr(wrapper,'TARGET_ROOT',target)
    package=tmp_path/'singular_output_fixture';package.mkdir();(package/'__init__.py').write_text('')
    result=tmp_path/'result.json'
    (package/'entry.py').write_text('from pathlib import Path\nimport json,sys\n'
        f'Path({str(result)!r}).write_text(json.dumps(str(Path(sys.argv[1]).resolve())))\n')
    monkeypatch.chdir(tmp_path);native=Path.resolve
    wrapper.main(['--record',str(record),'--record-sha256',digest,'-m',
        'singular_output_fixture.entry',str(source/'state')])
    assert json.loads(result.read_text())==str(source/'state')
    assert Path.resolve is native and source.resolve()==target
    assert wrapper.base.SOURCE_ROOT.name=='outputs'
    with pytest.raises(ValueError,match='SHA-256'):
        wrapper.main(['--record',str(record),'--record-sha256','0'*64,'-m','singular_output_fixture.entry'])
    assert Path.resolve is native and wrapper.base.SOURCE_ROOT.name=='outputs'
