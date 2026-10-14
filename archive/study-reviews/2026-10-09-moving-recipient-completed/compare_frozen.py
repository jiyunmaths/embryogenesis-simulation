"""Descriptive comparison with the exact matching frozen subset, after review."""
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).parent
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    review_file = HERE/'independent_review.json'
    review = json.loads(review_file.read_text())
    assert review['passed'] and review['accepted_at_tested_scope']
    p = json.loads(Path('output/studies/moving-recipient-amount-response/protocol.json').read_text())
    frozen_file = Path(p['accepted_delayed_parent'])/'summary.json'
    assert sha(frozen_file) == p['input_sha256'][str(frozen_file.resolve())]
    frozen = json.loads(frozen_file.read_text())
    assert frozen['passed']
    effects = frozen['response_summaries']['fine']['effects']
    rows = []
    for moving in review['response_comparisons']:
        corresponding = [r for r in effects if
            (r['seed'],r['cell'],r['dose'],r['matching_reference'],r['kind'],r['condition']) ==
            (moving['seed'],moving['cell'],moving['dose'],'S','chemical','sham_physics')]
        assert len(corresponding) == 1
        old = corresponding[0]
        a,b = old['normalized_response_shift_rms'],moving['effects']['fine']
        assert (a>.01) == old['response_shift_detected']
        rows.append(dict(history=moving['seed'],cell=moving['cell'],dose=moving['dose'],
            frozen_E=a,moving_E=b,fractional_E_change=b/a-1,
            same_screen_decision=bool((a>.01)==(b>.01))))
    assert len(rows)==12
    result = dict(passed=True,matching_reference='S',initial_physics='S',rows=rows,
        decisions_agree=sum(r['same_screen_decision'] for r in rows),comparisons=12,
        independent_histories=3,new_histories=0,
        review_sha256=sha(review_file),frozen_summary_sha256=sha(frozen_file),script_sha256=sha(__file__),
        interpretation='Descriptive matched-subset comparison, not a new moving-versus-frozen equivalence gate or isolation of geometry mediation. Both begin from the accepted t510 matched preparations; frozen physics remains fixed while the new study co-evolves mechanics, polarity, volumes, dilution and transport. Target chemistry is not rematched at t570.')
    (HERE/'frozen-comparison.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
