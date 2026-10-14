"""Run a recorded module using the relocated study tree's original path names.

Uses the existing tested, SHA-pinned migration alias without changing any
scientific source. The older plural outputs launcher remains unchanged.
"""
import importlib.util
from pathlib import Path

SOURCE_ROOT = Path('/home/sji/Documents/Programming/Python/embryogenesis-simulation/output/studies')
TARGET_ROOT = Path('/mnt/data/Documents/Python/embryogenesis-simulation/output/studies')
spec = importlib.util.spec_from_file_location('study_output_alias_base',
    Path(__file__).with_name('with_output_path_alias.py'))
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)


def main(argv=None):
    old_source, old_target = base.SOURCE_ROOT, base.TARGET_ROOT
    base.SOURCE_ROOT, base.TARGET_ROOT = SOURCE_ROOT, TARGET_ROOT
    try:
        return base.main(argv)
    finally:
        base.SOURCE_ROOT, base.TARGET_ROOT = old_source, old_target


if __name__ == '__main__':
    main()
