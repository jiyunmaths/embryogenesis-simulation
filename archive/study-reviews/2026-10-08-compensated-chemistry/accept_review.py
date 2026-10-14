"""Freeze independently accepted backend evidence; do not relax any gate."""
import hashlib,json
from pathlib import Path
from datetime import datetime,timezone

def sha(file):return hashlib.sha256(Path(file).read_bytes()).hexdigest()
def read(file):return json.loads(Path(file).read_text())
def main():
    root=Path('output/studies/compensated-chemistry-validation');archive=Path(__file__).parent
    r=read(archive/'independent_review.json');p=read(root/'protocol.json');s=read(root/'summary.json')
    assert r['passed'] and r['inputs_unchanged'] and s['passed'] and read(root/'status.json')['state']=='completed'
    assert r['protocol_sha256']==sha(root/'protocol.json') and all(x['passed'] for x in r['prefix'])
    assert r['diagnosed_control']['passed'] and r['method']==p['solver']['method']==s['method']
    files={**r['inputs_sha256'],**p['source_sha256'],**p['input_sha256']}
    for file,h in files.items():assert sha(file)==h,file
    for folder in (archive,Path('archive/study-launches/2026-10-08-compensated-chemistry')):
        for file in folder.rglob('*'):
            if file.is_file() and '__pycache__' not in file.parts:files[str(file.resolve())]=sha(file)
    for file in (root/'timing.json',root/'status.json',root/'summary.json',root/'prefix/summary.json'):
        files[str(file.resolve())]=sha(file)
    result=dict(passed=True,method=r['method'],accepted_utc=datetime.now(timezone.utc).isoformat(),scope='Isolated frozen-chemistry arithmetic only, exact existing delayed inputs. No blanket validation of moving chemistry or changed biology.',protocol_sha256=r['protocol_sha256'],independent_review=str((archive/'independent_review.json').resolve()),independent_review_sha256=sha(archive/'independent_review.json'),counts=r['counts'],diagnosed_control=r['diagnosed_control'],unit_tests=13,criteria=p['criteria'],files_sha256=files)
    out=Path('docs/compensated_chemistry_validation.json');assert not out.exists();out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('passed','method','counts','diagnosed_control')},indent=2))

if __name__=='__main__':main()
