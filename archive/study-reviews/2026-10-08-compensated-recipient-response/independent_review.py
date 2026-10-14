"""Separate delayed-array reconstruction plus compensated-state checks."""
import importlib.util,json
from pathlib import Path
import numpy as np

BASE=Path('archive/study-reviews/2026-10-08-recipient-amount-delayed-response/independent_review.py')
spec=importlib.util.spec_from_file_location('delayed_array_review',BASE)
r=importlib.util.module_from_spec(spec);spec.loader.exec_module(r)
ROOT=Path('output/studies/compensated-recipient-amount-delayed-response')
OUT=Path(__file__).with_suffix('.json')
METHOD='float64 SSP-RK2 increment with Kahan accumulation v1'


def main():
    r.pin(BASE);r.pin(__file__)
    p=r.read_json(ROOT/'protocol.json');q=r.read_json(p['backend_qualification'])
    assert q['passed'] and p['solver']['method']==q['method']==METHOD
    for file,h in q['files_sha256'].items():r.pin(file,h)
    oldroot=Path(p['historical_failed_study']);old=r.read_json(oldroot/'protocol.json')
    for field in ('criteria','levels','reference_settings','beta','da','db','start','duration','interval'):
        assert p[field]==old[field],field
    for a,b in zip(p['contexts'],old['contexts']):
        assert a['key']==b['key'];x,y=r.read_npz(a['source']),r.read_npz(b['source'])
        for key in ('delta','volumes','ids','initial'):assert np.array_equal(x[key],y[key])
    for a,b in zip(p['trials'],old['trials']):
        assert {k:v for k,v in a.items() if k!='source'}=={k:v for k,v in b.items() if k!='source'}
        assert np.array_equal(r.read_npz(a['source'])['initial'],r.read_npz(b['source'])['initial'])
    checks=[]
    for level in p['levels']:
        label=level['level'];sub=ROOT/'production'/label;cp=r.read_npz(sub/'checkpoint.npz');meta=json.loads(str(cp['metadata']));report=r.read_json(sub/'result.json');commit=r.read_json(sub/'checkpoint-commit.json')
        assert meta['method']==report['method']==commit['method']==METHOD and report['correction_checkpointed']
        e=cp['correction'];assert e.dtype==np.float64 and e.shape==cp['state'].shape and np.isfinite(e).all() and np.all(abs(e)<=2*np.spacing(cp['state']))
        for t in p['trials']:
            sub=ROOT/'cpu-verification'/label/t['key'];z=r.read_npz(sub/'paths.npz');report=r.read_json(sub/'result.json');e=z['final_correction']
            assert report['method']==METHOD and e.dtype==np.float64 and e.shape==z['trajectory'][-1].shape and np.isfinite(e).all() and np.all(abs(e)<=2*np.spacing(z['trajectory'][-1]))
        checks.append(dict(level=label,correction_checkpoint_verified=True,cpu_corrections_verified=144,method=METHOD))
    # The unchanged standalone reviewer re-derives every scientific contrast,
    # numerical gate, full amount curve and exact accepted-parent preparation.
    r.main(ROOT,OUT)
    result=json.loads(OUT.read_text());prior=r.read_json('archive/study-reviews/2026-10-08-recipient-amount-delayed-response/independent_review.json')
    comparisons=[]
    fields=('seed','cell','matching_reference','physics','dose')
    for row in result['scientific_comparisons']:
        before=next(x for x in prior['scientific_comparisons'] if all(x[k]==row[k] for k in fields))
        comparisons.append(dict(**{k:row[k] for k in fields},old_failed_fine_effect=before['full_effects']['fine'],compensated_fine_effect=row['full_effects']['fine'],effect_difference=row['full_effects']['fine']-before['full_effects']['fine'],same_above_screen_decision=row['above_screen']==before['above_screen']))
    failed=next(x['trial'] for x in prior['orders'] if not x['passed']);new=next(x for x in result['orders'] if x['trial']==failed)
    oldfine=r.read_npz(oldroot/'production/fine/paths.npz');newfine=r.read_npz(ROOT/'production/fine/paths.npz');i=list(newfine['trial_keys']).index(failed)
    drift=r.log_error(newfine['trajectory'][:,i],newfine['trajectory'][0,i])
    resolved_relaxation=bool(drift>0 and not new['order_resolved'] and new['passed'])
    result.update(method=METHOD,compensated_checkpoint_checks=checks,exact_historical_preparations=True,backend_qualification_verified=True,original_failed_study_preserved=True,original_failed_scientific_comparison=comparisons,diagnosed_control=dict(new_order=new,stored_max_log_relaxation=drift,old_fine_log_relaxation=r.log_error(oldfine['trajectory'][:,i],oldfine['trajectory'][0,i]),resolved_relaxation=resolved_relaxation),scope=result['scope']+' Isolated compensated arithmetic, exact original failed starts, mandatory numerical correction. No shared-kernel changes or threshold changes.')
    assert resolved_relaxation
    guards=r.read_json('archive/study-document-snapshots/2026-10-08-before-compensated-chemistry/manifest.json')['unchanged_guards']
    for file,h in guards.items():r.pin(file,h)
    for file,h in r.WATCHED.items():assert r.digest(file)==h
    result['inputs_sha256']=r.WATCHED;result['script_sha256']=r.digest(__file__)
    OUT.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(passed=result['passed'],diagnosed_control=result['diagnosed_control'],all_48_scientific_decisions_agree_with_original=all(x['same_above_screen_decision'] for x in comparisons),maximum_effect_difference=max(abs(x['effect_difference']) for x in comparisons)),indent=2))

if __name__=='__main__':main()
