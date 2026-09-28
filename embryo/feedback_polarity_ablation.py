"""Complete the matched component study with tension + adhesion, no polar tension."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
from pathlib import Path
import shutil
from . import feedback_long as original
from .resolution import write_json

ARM = 'tension_adhesion'
COEFFICIENTS = (.25, .35, 0.)


def initialize(checkpoint, signals):
    sim = original.initialize(checkpoint, signals, 'full')
    sim.config.polarity_tension = 0.
    sim.config.validate()
    return sim


def prepare(source, output):
    source, output = Path(source).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError('Choose a fresh directory')
    p = json.loads((source/'protocol.json').read_text())
    status = json.loads((source/'status.json').read_text())
    if status['state'] != 'completed' or status['completed'] != 10:
        raise ValueError('Requires completed ten-arm reference study')
    if original.digest(source/'protocol.json') != status['protocol_sha256']:
        raise ValueError('Reference protocol changed')
    inputs = {**p['source_sha256'], p['checkpoint']:p['checkpoint_sha256'],
              str(source/'initial_states.npz'):p['initial_states_sha256']}
    for path, expected in inputs.items():
        if original.digest(path) != expected:
            raise ValueError('Reference input changed: '+path)
    for state in ['formation', 'persistence']:
        for arm in original.ARMS:
            result = json.loads((source/f'{state}_{arm}'/'result.json').read_text())
            if not result['completed'] or not result['quality_pass']:
                raise ValueError('Reference failed quality checks')
    if tuple(p['arms']['full']) != (.25, .35, .35):
        raise ValueError('Unexpected full-feedback coefficients')
    output.mkdir(parents=True)
    shutil.copyfile(source/'initial_states.npz', output/'initial_states.npz')
    p['arms'] = {ARM: list(COEFFICIENTS)}
    p['jobs'] = [dict(state=s, arm=ARM) for s in ['formation','persistence']]
    p['reference'] = str(source)
    reference_files = [source/'protocol.json', source/'comparison.json']
    reference_files += [source/f'{s}_{a}'/name for s in ['formation','persistence']
                        for a in original.ARMS for name in ['history.json','result.json']]
    p['reference_sha256'] = {str(f):original.digest(f) for f in reference_files}
    p['source_sha256'][str(Path(__file__).resolve())] = original.digest(Path(__file__))
    p['scope'] = ('Exact starting geometry and chemical arrays from completed component study. '
                  'Full coupling with only polarity_tension set to zero. Polarity dynamics and '
                  'chemical modulation of polarity remain active. Tests necessity of polar mechanical '
                  'action for sustained-contrast suppression in this one matched setting, not general necessity.')
    p['assessment'] = dict(primary='formation late minimum log-activator SD > 0.1 over time 63–78',
                           secondary='persistence using the identical frozen-prepared reference state',
                           no_retuning=True)
    write_json(output/'protocol.json', p)
    write_json(output/'status.json', dict(state='prepared', protocol_sha256=original.digest(output/'protocol.json')))
    return p


def worker(args):
    # Extend only this worker process's arm registry. Production source files and
    # the running developmental studies remain unchanged. The mature-geometry
    # initializer, every-step quality audit, and resume mechanism are reused.
    if ARM in original.ARMS:
        raise ValueError('Unexpected arm registration')
    original.ARMS[ARM] = COEFFICIENTS
    try:
        return original.worker(args)
    finally:
        del original.ARMS[ARM]


def assess(formation, reference_full, reference_baseline):
    if not all(r['quality_pass'] for r in [formation,reference_full,reference_baseline]):
        return 'inconclusive_numerical_quality'
    if reference_full['persistent'] or not reference_baseline['persistent']:
        return 'inconclusive_reference_outcomes'
    return ('removing_polar_mechanics_restores_formation' if formation['persistent'] else
            'suppression_persists_without_polar_mechanics')


def compare(output, results):
    output = Path(output); p = json.loads((output/'protocol.json').read_text())
    reference = Path(p['reference'])
    new = {r['job']['state']:r for r in results}
    refs = {s:{a:json.loads((reference/f'{s}_{a}'/'result.json').read_text())
               for a in original.ARMS} for s in ['formation','persistence']}
    verdict = assess(new['formation'],refs['formation']['full'],refs['formation']['baseline'])
    interpretation = {
        'removing_polar_mechanics_restores_formation':
            'Removing polarity mechanics rescues the declared formation outcome with tension and adhesion retained. Supports its necessity for suppression by the tested full combination on this initial geometry and horizon.',
        'suppression_persists_without_polar_mechanics':
            'Tension and adhesion together still suppress the declared formation outcome without polarity mechanics. Polarity mechanics is not necessary for suppression in this tested setting.',
    }.get(verdict, 'No causal verdict: reference outcomes or numerical quality do not support the comparison.')
    report = dict(verdict=verdict,interpretation=interpretation,new_arms=new,reference=refs,
                  scope=p['scope'])
    write_json(output/'comparison.json',report)
    lines=['# Tension + adhesion without polarity mechanics','',interpretation,'',
           '| Arm | Formation: late minimum log SD | Persistence: late minimum log SD |',
           '|---|---:|---:|']
    for arm in ['baseline','tension','adhesion',ARM,'polarity','full']:
        f = new['formation'] if arm == ARM else refs['formation'][arm]
        q = new['persistence'] if arm == ARM else refs['persistence'][arm]
        lines.append(f"| {arm} | {f['late_min_log_activator_sd']:.6g} | {q['late_min_log_activator_sd']:.6g} |")
    lines += ['', 'The threshold is 0.1 throughout time 63–78. This is one mature geometry and one chemical perturbation seed; it does not establish general necessity, inherited identity, or continuum consistency.']
    (output/'RESULTS.md').write_text('\n'.join(lines)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2,2,figsize=(11,7),layout='constrained')
    for col,state in enumerate(['formation','persistence']):
        for arm in ['baseline','full','polarity',ARM]:
            folder = output if arm == ARM else reference
            h = json.loads((folder/f'{state}_{arm}'/'history.json').read_text())
            t = [r['metrics']['time'] for r in h]
            axes[0,col].plot(t,[r['attribute_std'][0] for r in h],label=arm)
            axes[1,col].plot(t,[r['spectrum']['maximum_spatial_growth'] for r in h],label=arm)
        axes[0,col].set(title=state,ylabel='SD of log activator');axes[0,col].axhline(.1,color='gray',ls=':')
        axes[1,col].set(xlabel='Model time',ylabel='Frozen homogeneous growth');axes[1,col].axhline(0,color='gray',ls=':')
        axes[0,col].legend(fontsize=8)
    fig.savefig(output/'comparison.png',dpi=160);plt.close(fig)
    return report


def run(output):
    output = Path(output);p = json.loads((output/'protocol.json').read_text())
    ph = original.digest(output/'protocol.json')
    if ph != json.loads((output/'status.json').read_text())['protocol_sha256']:
        raise ValueError('Protocol changed')
    expected = {**p['source_sha256'], **p['reference_sha256'],
                p['checkpoint']:p['checkpoint_sha256'],
                str(output/'initial_states.npz'):p['initial_states_sha256']}
    for path,h in expected.items():
        if original.digest(path) != h:
            raise ValueError('Frozen input changed: '+path)
    rows,failures = [],[]
    write_json(output/'status.json',dict(state='running',protocol_sha256=ph,completed=0,total=2))
    with ProcessPoolExecutor(max_workers=2) as pool:
        futures = {pool.submit(worker,(str(output),job)):job for job in p['jobs']}
        for f in as_completed(futures):
            try: rows.append(f.result())
            except Exception as error: failures.append(dict(job=futures[f],error=str(error)))
            write_json(output/'status.json',dict(state='running',protocol_sha256=ph,completed=len(rows),total=2,failures=failures))
    if not failures: compare(output,rows)
    write_json(output/'status.json',dict(state='failed' if failures else 'completed',protocol_sha256=ph,completed=len(rows),total=2,failures=failures))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['prepare','run'])
    parser.add_argument('--source',type=Path,default=Path('outputs/feedback-long'))
    parser.add_argument('--output',type=Path,default=Path('outputs/feedback-polarity-ablation'))
    args=parser.parse_args()
    if args.command == 'prepare': prepare(args.source,args.output)
    else: run(args.output)
