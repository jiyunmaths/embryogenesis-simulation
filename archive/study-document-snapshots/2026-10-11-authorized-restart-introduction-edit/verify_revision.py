from pathlib import Path
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import re
import sys
import numpy as np

ROOT = Path('/home/sji/Documents/Programming/Python/embryogenesis-simulation')
MAN = ROOT / 'manuscript'
ARCH = MAN / 'archive/before-restart-introduction-2026-10-11'

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def intro(text):
    return text.split('# Introduction\n', 1)[1].split('# Methods\n', 1)[0]

before = (ARCH / 'sections.md').read_text()
after = (MAN / 'sections.md').read_text()
assert before.split('# Introduction\n')[0] == after.split('# Introduction\n')[0]
assert before.split('# Methods\n', 1)[1] == after.split('# Methods\n', 1)[1]
for item in json.loads((ARCH / 'manifest.json').read_text())['files']:
    assert sha(ROOT / item['snapshot']) == item['sha256']
assert sha(MAN / 'references.bib') == sha(ARCH / 'references.bib')
assert sha(MAN / 'main.pdf') == sha(ARCH / 'main.pdf')
assert set(re.findall(r'@([A-Za-z0-9_:-]+)', intro(before))) == set(re.findall(r'@([A-Za-z0-9_:-]+)', intro(after)))
assert 'resulting frozen network' in intro(after)
assert 'may induce a relaxation transient' in intro(after)
assert 'cannot distinguish baseline mechanical relaxation' in intro(after)
for name in ['sections.md', 'README.md']:
    assert all(line == line.rstrip() for line in (MAN / name).read_text().splitlines())
html = (MAN / 'introduction_methods.html').read_text()
md = (MAN / 'introduction_methods.md').read_text()
tex = (MAN / 'sections.tex').read_text()
assert 'citeproc-not-found' not in html
assert 'may induce a relaxation transient' in md and 'may induce a relaxation transient' in tex

sources = ['embryo/phase_carry_polarity.py', 'embryo/polarity_robustness.py',
           'embryo/parameter_robustness_moving.py', 'embryo/parameter_robustness.py',
           'embryo/feedback_survival.py', 'embryo/attribute_development.py',
           'embryo/polarity_conductance_controls.py', 'embryo/moving_initiation_controls.py']
code = (ROOT / 'embryo/parameter_robustness_moving.py').read_text()
assert "g['branch'] != 'switch_on'" in code
assert 'host.activator, host.inhibitor = state' in code
assert "not host.config.feedback" in code
survival = (ROOT / 'embryo/feedback_survival.py').read_text()
assert "if arm=='switch_on':" in survival
assert "sim.attribute_mode='direct';sim.config.feedback=True" in survival

carry_path = ROOT / 'outputs/phase-carry-polarity/protocol.json'
carry = json.loads(carry_path.read_text())
assert {(j['family'], j['point']['chi']) for j in carry['jobs']} == {('uniform', 0.), ('uniform', .35)}
assert carry['accepted_configs']['fine']['feedback']
assert carry['accepted_configs']['fine']['fate_tension'] == .25
assert carry['accepted_configs']['fine']['fate_adhesion'] == .35
starts = []
for seed in (7, 8, 9):
    job = next(j for j in carry['jobs'] if j['seed'] == seed and j['level'] == 'fine' and j['point']['chi'] == 0.)
    with np.load(job['source'], allow_pickle=False) as z:
        metadata = json.loads(str(z['metadata'].item()))
        old_a = z['activator'].copy()
        ids = z['ids'].copy()
        mode = str(z['attribute_mode'].item())
    with np.load(job['chemical_file'], allow_pickle=False) as z:
        reset_a = z['uniform'][0].copy()
        assert np.array_equal(ids, z['ids'])
    c = metadata['config']
    assert metadata['time'] == 150. and c['feedback'] and mode == 'direct'
    assert c['polarity_tension'] == .35
    assert np.std(np.log(old_a)) > .1 and np.std(np.log(reset_a)) < .01
    old_r, reset_r = np.tanh(old_a - 1.), np.tanh(reset_a - 1.)
    tension_change = c['fate_tension'] * (reset_r - old_r)
    adhesion_change = c['fate_adhesion'] * (np.outer(reset_r, reset_r) - np.outer(old_r, old_r))
    np.fill_diagonal(adhesion_change, 0.)
    assert abs(tension_change).max() > 0 and abs(adhesion_change).max() > 0
    starts.append(dict(history=seed, source=job['source'], source_sha256=sha(job['source']),
                      chemical_file=job['chemical_file'], chemical_file_sha256=sha(job['chemical_file']),
                      time=metadata['time'], feedback=c['feedback'], mode=mode,
                      inherited_chi=c['polarity_tension'],
                      original_log_activator_sd=float(np.std(np.log(old_a))),
                      reset_log_activator_sd=float(np.std(np.log(reset_a))),
                      tension_and_adhesion_changed_by_reset=True,
                      other_preparation_changes='The assay also changes D_b from the inherited .4 to .55; chi=0 removes the inherited .35 directional action. These are declared matched interventions, not uninterrupted development.'))

evidence_paths = ['docs/phase_carry_polarity_assessment.md',
                  'docs/polarity_conductance_controls_assessment.md',
                  'docs/moving_initiation_controls_assessment.md',
                  'docs/phase_carry_maintenance_assessment.md',
                  'docs/moving_geometry_initiation_controls.md']
evidence = {p: sha(ROOT / p) for p in evidence_paths}

sys.path.insert(0, str(ROOT))
spec = importlib.util.spec_from_file_location('study_alias', ROOT / 'scripts/with_study_output_path_alias.py')
alias = importlib.util.module_from_spec(spec)
spec.loader.exec_module(alias)
with alias.base.validated_alias(
    ROOT / 'archive/output-relocation/2026-10-10-output/studies-manifest.json',
    '8fb2545d8bae35a1a673ec3a43d3cc191148f03cac58892f29bd368182a8aff6',
    source_root=alias.SOURCE_ROOT, target_root=alias.TARGET_ROOT,
):
    from embryo import qualify_repeated_development as q
    active_path = ROOT / 'output/studies/fresh-zygote-pilot/protocol.json'
    active = q.read(active_path)
    q.verify(active)

utc = datetime.now(timezone.utc).isoformat()
manifest = json.loads((MAN / 'evidence_manifest.json').read_text())
original_manifest = json.loads((ARCH / 'evidence_manifest.json').read_text())
for key in ['reports', 'table_values_verified', 'history_registry', 'claim_precision_audit']:
    assert manifest[key] == original_manifest[key]
changes = {}
for name in ['sections.md', 'README.md']:
    key = 'manuscript/' + name
    changes[key] = dict(before=sha(ARCH / name), after=sha(MAN / name))
    manifest['files_sha256'][key] = sha(MAN / name)
manifest['editorial_revisions'].append(dict(
    utc=utc,
    scope='User-authorized introduction-only mechanical/restart qualification and formation-maintenance connection; no new numerical results.',
    archive=str(ARCH.relative_to(ROOT)), manuscript_hash_updates=changes,
    verification='manuscript/restart_introduction_verification.json',
    original_numerical_evidence_retained=True, manuscript_wide_evidence_migration_pending=True,
))
(MAN / 'evidence_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')

# Preserve the earlier non-input manuscript pins explicitly; leave the scientific manifest intact.
parent = Path(active['preserved_manifest'])
pins = json.loads(parent.read_text())['preserved_sha256']
relocations = []
for file, expected in pins.items():
    path = Path(file)
    if path.parent == MAN and path.exists() and sha(path) != expected:
        saved = MAN / 'archive/before-carry-introduction-2026-10-11' / path.name
        assert sha(saved) == expected, path
        assert file not in active['source_sha256'] and file not in active['input_sha256']
        relocations.append(dict(original_path=file, expected_sha256=expected,
                                preserved_bytes=str(saved), edited_sha256=sha(path)))
record_root = ROOT / 'archive/study-document-snapshots/2026-10-11-authorized-restart-introduction-edit'
record_root.mkdir(exist_ok=False)
(record_root / 'record.json').write_text(json.dumps(dict(
    utc=utc,
    authorization='User explicitly requested introduction refinements to mechanical attribution, mature reset confounding, the missing control and the coexistence explanation.',
    parent_preservation_manifest=str(parent), parent_preservation_manifest_sha256=sha(parent),
    non_input_document_relocations=relocations, source_and_input_pins_unchanged=True,
    original_numerical_evidence_unchanged=True,
), indent=2) + '\n')

record = dict(
    checked_utc=utc, checks_passed=True,
    scope='Introduction-only refinement plus guide/evidence notes and regenerated Markdown/HTML/LaTeX. No new simulation or PDF build.',
    archive=str(ARCH.relative_to(ROOT)),
    abstract_and_methods_onward_unchanged=True, bibliography_unchanged=True, original_pdf_unchanged=True,
    citation_keys_unchanged=True, build_passed=True,
    preparation_chain='No-feedback patterned development to t=90; switch_on feedback/patterned continuation to t=150; frozen endpoint study selects switch_on graph; moving robustness copies that checkpoint; polarity robustness copies/retimes it; carry initiation copies/retimes it and loads near-uniform chemistry.',
    source_sha256={p: sha(ROOT / p) for p in sources},
    preparation_checks=starts,
    qualified_assessment_sha256=evidence,
    missing_control_scope='No matched all-mechanical-feedback-off branch in this carry-qualified mature near-uniform initiation series. Historical no-carry keep_off developed-state continuations and fresh developmental feedback-off controls are different preparations, not this control.',
    interpretation='Reset-induced material changes are demonstrated from saved starts and the constitutive law; an actual causal contribution of a restart transient, and the decomposition into passive versus activity-dependent remodeling, remain unmeasured. Total conductance decline is not a universal initiation threshold.',
    coexistence_scope='Linear uniform stability concerns a frozen operator; a finite-amplitude patterned attractor may coexist. No full moving-system stability theorem is asserted.',
    active_study=dict(protocol_sha256=sha(active_path), source_hashes_checked=len(active['source_sha256']), input_hashes_checked=len(active['input_sha256']), unchanged=True),
    historical_preservation_record=str((record_root / 'record.json').relative_to(ROOT)),
    manuscript_sha256={str((MAN / name).relative_to(ROOT)): sha(MAN / name) for name in ['sections.md', 'sections.tex', 'README.md', 'introduction_methods.md', 'introduction_methods.html', 'evidence_manifest.json']},
)
(MAN / 'restart_introduction_verification.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps(dict(checks_passed=True, source_histories_checked=len(starts),
                     reset_changes_tension_and_adhesion=True, abstract_and_methods_onward_unchanged=True,
                     bibliography_and_pdf_unchanged=True,
                     active_source_and_input_pins=len(active['source_sha256']) + len(active['input_sha256'])), indent=2))
