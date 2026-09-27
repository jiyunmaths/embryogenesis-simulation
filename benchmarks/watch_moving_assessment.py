"""Write the assessment automatically once the concurrent pilot finishes."""
import argparse
import json
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent))
from assess_moving_causal import assess
from embryo.resolution import write_json


def watch(output, controller_pid):
    output = Path(output)
    target = output / 'assessment-status.json'
    write_json(target, {'state': 'waiting', 'controller_pid': controller_pid})
    try:
        while True:
            status = json.loads((output / 'status.json').read_text())
            if status['state'] == 'completed':
                break
            if status['state'] == 'failed' or not Path(f'/proc/{controller_pid}').exists():
                raise RuntimeError('Simulation controller stopped before completion; inspect finish.log and worker logs')
            time.sleep(10)
        write_json(target, {'state': 'assessing'})
        result = assess(output)
        write_json(target, {'state': 'completed', 'report': str(output / 'ASSESSMENT.md'),
                   'all_quality_pass': result['original_comparison']['all_quality_pass']})
        print((output / 'ASSESSMENT.md').read_text(), flush=True)
    except Exception as error:
        write_json(target, {'state': 'failed', 'error': str(error)})
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--controller-pid', type=int, required=True)
    args = parser.parse_args()
    watch(args.output, args.controller_pid)
