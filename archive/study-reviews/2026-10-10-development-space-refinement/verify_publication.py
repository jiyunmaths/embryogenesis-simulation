"""Verify unchanged scientific pins and local links after publishing the review."""
from pathlib import Path
import datetime
import hashlib
import json
import re
import time

HERE=Path(__file__).parent
started=time.monotonic()
def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()
rpath=HERE/'independent_review.json'
r=json.loads(rpath.read_text())
assert r['review_completed'] and r['scoped_screen_accepted']
assert not r['qualification_accepted'] and not r['full_development_qualified']
for p,expected in r['consumed_sha256'].items():
    assert sha(p)==expected, f'Changed scientific pin: {p}'
a=json.loads(Path('docs/development_space_refinement_assessment.json').read_text())
assert a['independent_review']['sha256']==sha(rpath)
for key in ('space_rows','time_row','paths','protocol_sha256','scoped_screen_accepted','new_scientific_histories'):
    assert a[key]==r[key], key
fpath=Path('docs/figures/development-space-refinement-summary.json')
f=json.loads(fpath.read_text())
for p,expected in f['inputs_sha256'].items():
    assert sha(p)==expected, p
texts={
    Path('docs/development_space_refinement_assessment.md'):Path('docs/development_space_refinement_assessment.md').read_text(),
    Path('docs/development_space_refinement_launch.md'):Path('docs/development_space_refinement_launch.md').read_text(),
}
for name,prefix in [('README.md','**Feedback-off spatial and finest-grid time screens completed:**'),('docs/plan.md','**Feedback-off spatial/finer-grid time screen completed and independently confirmed:**')]:
    p=Path(name)
    texts[p]=next(line for line in p.read_text().splitlines() if line.startswith(prefix))
checked_links=[]
for p,content in texts.items():
    for href in re.findall(r'\]\(([^)]+)\)',content):
        if href.startswith(('https:','http:','#')):continue
        target=(p.parent/href.split('#')[0]).resolve()
        assert target.exists(), f'Missing link {p}: {href}'
        checked_links.append({'document':str(p),'target':href})
artifacts=[rpath,HERE/'review_saved_arrays.py',HERE/'review.log',HERE/'plot_summary.py',Path('docs/development_space_refinement_assessment.md'),Path('docs/development_space_refinement_assessment.json'),Path('docs/figures/development-space-refinement-summary.png'),Path('docs/figures/development-space-refinement-summary.svg'),fpath,Path('README.md'),Path('docs/plan.md'),Path('docs/development_space_refinement_launch.md')]
result={'verified_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'unchanged_consumed_pins':len(r['consumed_sha256']),'preserved_evidence_pins':r['preserved_count'],'links_checked':checked_links,'artifacts_sha256':{str(p):sha(p) for p in artifacts},'new_jobs_launched':False,'scoped_screen_accepted':True,'full_development_qualified':False,'verification_seconds':time.monotonic()-started}
(HERE/'final-verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['unchanged_consumed_pins','preserved_evidence_pins','scoped_screen_accepted','full_development_qualified','verification_seconds']}))
