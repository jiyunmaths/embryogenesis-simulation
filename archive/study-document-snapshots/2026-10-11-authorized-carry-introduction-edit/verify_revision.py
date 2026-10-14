from pathlib import Path
from datetime import datetime, timezone
from collections import Counter
import hashlib
import importlib.util
import json
import re
import sys

ROOT = Path('/home/sji/Documents/Programming/Python/embryogenesis-simulation')
MAN = ROOT / 'manuscript'
ARCHIVE = MAN / 'archive/before-carry-introduction-2026-10-11'

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def introduction(text):
    return text.split('# Introduction\n', 1)[1].split('# Methods\n', 1)[0]

before = (ARCHIVE / 'sections.md').read_text()
after = (MAN / 'sections.md').read_text()
assert before.split('# Introduction\n')[0] == after.split('# Introduction\n')[0]
assert before.split('# Methods\n', 1)[1] == after.split('# Methods\n', 1)[1]
old_intro, new_intro = introduction(before), introduction(after)
old_citations = set(re.findall(r'@([A-Za-z0-9_:-]+)', old_intro))
new_citations = set(re.findall(r'@([A-Za-z0-9_:-]+)', new_intro))
bib = (MAN / 'references.bib').read_text()
old_bib = (ARCHIVE / 'references.bib').read_text()
entries = re.findall(r'^@[A-Za-z]+\{([^,]+),', bib, re.M)
assert len(entries) == len(set(entries))
assert set(re.findall(r'@([A-Za-z0-9_:-]+)', after)) <= set(entries)
assert old_citations <= new_citations
for key in ['othmer1971', 'petit2017']:
    bib = re.sub(r'^@article\{' + key + r',\n.*?^\}\n', '', bib, flags=re.M | re.S)
assert bib == old_bib, 'Previously registered bibliography entries changed'

for item in json.loads((ARCHIVE / 'manifest.json').read_text())['files']:
    assert digest(ROOT / item['snapshot']) == item['sha256'], item['snapshot']
for path in ['main.pdf', 'main.tex', 'build.py']:
    assert digest(MAN / path) == digest(ARCHIVE / path)

html = (MAN / 'introduction_methods.html').read_text()
md = (MAN / 'introduction_methods.md').read_text()
tex = (MAN / 'sections.tex').read_text()
assert 'citeproc-not-found' not in html
assert not re.search(r'\[@[A-Za-z]', md)
for key in new_citations:
    assert 'id="ref-' + key + '"' in html, key
assert 'Histories 7 and 8 fail to initiate at either setting.' in md
assert 'Histories 7 and 8 fail to initiate at either setting.' in tex
assert 'instability and dynamic pattern in cellular networks' in md.lower()
assert 'theory of turing patterns on time varying networks' in md.lower()
assert 'One matched mature starting geometry' not in new_intro

support = [
    ('docs/phase_carry_polarity_assessment.md',
     'Only history 9 has matched polarity-specific initiation suppression; histories 7 and 8 fail with either directional-tension setting.'),
    ('docs/polarity_conductance_controls_assessment.md',
     'Preserving initial symmetric G while retaining actual M(t), mechanics and dilution restores formation in all three histories at both settings; spectra support a contribution, not complete mediation or full moving stability.'),
    ('docs/phase_carry_maintenance_assessment.md',
     'Developed chemistry persists for 240 elapsed units at both settings in all three histories; all developed-start frozen endpoints support local chemical coexistence.'),
    ('docs/phase_carry_exchange_response_assessment.md',
     'Carry-qualified mature exchange responses are donor-nearer in all three histories, without exact concentration retention or donor equivalence.'),
    ('docs/moving_pulse_time_amount_response_assessment.md',
     'Immediate moving context-response differences at matched initial recipient chemistry, species amounts and physical fields.'),
    ('docs/moving_recipient_amount_response_assessment.md',
     'Selective delayed moving response differences after an initially amount-matched preparation; amounts and local chemistry evolve again before the challenge.'),
]
scientific_evidence = [{'path': path, 'sha256': digest(ROOT / path), 'claim_scope': scope}
                       for path, scope in support]

citations = [
    dict(key='othmer1971', doi='10.1016/0022-5193(71)90154-8',
         sources=['https://www.sciencedirect.com/science/article/pii/0022519371901548', 'https://pubmed.ncbi.nlm.nih.gov/5571122/'],
         checked='Publisher abstract and bibliographic record: Othmer and Scriven, J. Theor. Biol. 32(3), 507–537 (1971); instability in cellular compartment networks. Full text was not accessible.'),
    dict(key='petit2017', doi='10.1103/PhysRevLett.119.148301',
         sources=['https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.119.148301', 'https://arxiv.org/abs/1705.08025'],
         checked='APS metadata and abstract plus author preprint abstract: Petit, Lauwens, Fanelli and Carletti; PRL 119(14), 148301 (2017); imposed time-varying networks. No new full moving-system stability claim follows.'),
    dict(key='tseng2024', doi='10.1371/journal.pbio.3002636',
         sources=['https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3002636'],
         checked='Primary article: connexin/gap-junction interventions alter feather patterning; the identity of the transported inhibitory factor is unresolved.'),
    dict(key='liu2022', doi='10.1038/s41467-022-28149-3',
         sources=['https://www.nature.com/articles/s41467-022-28149-3'],
         checked='Primary article: Nodal activity spreads by production induced in neighboring cells; relay production is not symmetric conservative exchange.'),
    dict(key='autorino2026', doi='10.1038/s41556-026-01954-4',
         sources=['https://www.nature.com/articles/s41556-026-01954-4'],
         checked='Primary article: adhesion/rigidity interventions alter extracellular porosity, Nodal dispersal and Lefty timing; mechanism differs from cell-contact conductances.'),
    dict(key='yu2024', doi='10.1126/sciadv.adn0172',
         sources=['https://pmc.ncbi.nlm.nih.gov/articles/PMC11062584/'],
         checked='Primary full text: 3D vertex model with reactions/diffusion, mechanics and reversible network reconnection; also allows cell division. This contradicts a blanket fixed-network/continuum-only novelty claim.'),
    dict(key='mercker2013', doi='10.1371/journal.pone.0082617',
         sources=['https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0082617'],
         checked='Primary full text: morphogen expression coupled to curvature on a continuous deforming cell sheet.'),
    dict(key='nissen2018', doi='10.7554/eLife.38407',
         sources=['https://pmc.ncbi.nlm.nih.gov/articles/PMC6286147/'],
         checked='Primary full text: interacting cell particles and polarity-dependent forces, with changing neighborhood assignments; not a fixed mechanical neighborhood claim.'),
    dict(key='kuang2023', doi='10.1038/s41540-023-00265-w',
         sources=['https://www.nature.com/articles/s41540-023-00265-w'],
         checked='Primary article: phase-field multicellular morphology and cleavage; phase-field representation is established methodology.'),
    dict(key='vangorder2021', doi='10.1007/s00285-021-01552-y',
         sources=['https://arxiv.org/abs/1904.09683'],
         checked='Author preprint abstract: evolving manifolds and non-autonomous linearized dynamics; evolving domain specified rather than this coupled cell mechanics.'),
    dict(key='nakao2010', doi='10.1038/nphys1651',
         sources=['https://arxiv.org/abs/1005.1986'],
         checked='Author preprint: discrete network patterning, multistability and hysteresis.'),
    dict(key='wolfrum2012', doi='10.1016/j.physd.2012.05.002',
         sources=['https://oa.tib.eu/renate/handle/123456789/4020'],
         checked='Author preprint bibliographic record and abstract: differentiated-node branches before uniform instability, multistability and hysteresis.'),
    dict(key='hajali2025', doi='10.1371/journal.pcsy.0000044',
         sources=['https://journals.plos.org/complexsystems/article?id=10.1371/journal.pcsy.0000044'],
         checked='Primary article: spectral gaps and graph interventions affect pattern diversity; unweighted graph setting differs from changing conductances and cell volumes.'),
]

# Check pinned scientific inputs without stepping, launching, or stopping a simulation.
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
    protocol_path = ROOT / 'output/studies/fresh-zygote-pilot/protocol.json'
    assert digest(protocol_path) == '728f584676ea4eeeee530656015aec61931e9296ff095eff8b7c94f9832a13ee'
    protocol = q.read(protocol_path)
    q.verify(protocol)

original_evidence = json.loads((ARCHIVE / 'evidence_manifest.json').read_text())
current_evidence = json.loads((MAN / 'evidence_manifest.json').read_text())
for key in ['reports', 'table_values_verified', 'history_registry', 'claim_precision_audit']:
    assert current_evidence[key] == original_evidence[key]
update_hashes = {}
for name in ['README.md', 'sections.md', 'references.bib']:
    key = 'manuscript/' + name
    update_hashes[key] = {'before': digest(ARCHIVE / name), 'after': digest(MAN / name)}
    current_evidence['files_sha256'][key] = digest(MAN / name)
utc = datetime.now(timezone.utc).isoformat()
current_evidence.setdefault('editorial_revisions', []).append(dict(
    utc=utc,
    scope='Introduction-only switch to carry-qualified mature conductance, maintenance, exchange and response evidence; verified literature additions. Abstract and Methods onward retain original-method analyses.',
    archive=str(ARCHIVE.relative_to(ROOT)),
    manuscript_hash_updates=update_hashes,
    introduction_evidence=scientific_evidence,
    original_numerical_evidence_retained=True,
    manuscript_wide_evidence_migration_pending=True,
))
(MAN / 'evidence_manifest.json').write_text(json.dumps(current_evidence, indent=2) + '\n')

record = dict(
    checked_utc=utc,
    scope='User-authorized introduction and bibliography revision, guide qualification and regenerated Markdown/HTML/LaTeX previews. No PDF rebuild or quantitative evidence replacement.',
    archive=str(ARCHIVE.relative_to(ROOT)),
    abstract_unchanged=True,
    methods_results_discussion_appendix_unchanged=True,
    original_pdf_unchanged=True,
    original_bibliography_entries_unchanged=True,
    original_numerical_evidence_manifest_fields_unchanged=True,
    scientific_unit='Three reused histories 7, 8, 9; nested branches, graphs, cells and perturbations are not additional histories.',
    formation_scope='Mature 16-cell continuations with carry; inherited original developmental geometry; not fresh carry-qualified development.',
    introduction_evidence=scientific_evidence,
    citation_checks=citations,
    introduction_citation_keys=sorted(new_citations),
    added_introduction_citation_keys=sorted(new_citations - old_citations),
    added_bibliography_keys=['othmer1971', 'petit2017'],
    introduction_words_before=len(old_intro.split()),
    introduction_words_after=len(new_intro.split()),
    build='python manuscript/build.py: passed, no unresolved citation warnings; all introduction references present in generated HTML bibliography.',
    active_study=dict(protocol_sha256=digest(protocol_path), source_hashes_checked=len(protocol['source_sha256']), input_hashes_checked=len(protocol['input_sha256']), unchanged=True),
    manuscript_wide_evidence_migration_pending=['Abstract', 'Methods: numerical method and assay preparations/windows', 'Results tables and figures', 'Discussion alignment', 'New submission PDF'],
    manuscript_sha256={str((MAN / name).relative_to(ROOT)): digest(MAN / name) for name in ['sections.md', 'sections.tex', 'references.bib', 'introduction_methods.md', 'introduction_methods.html', 'README.md', 'evidence_manifest.json']},
    checks_passed=True,
)
(MAN / 'carry_introduction_verification.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps(dict(checked_utc=utc, checks_passed=True, abstract_and_methods_onward_unchanged=True,
                     original_pdf_unchanged=True, introduction_citations=len(new_citations),
                     added_citations=sorted(new_citations - old_citations),
                     evidence_assessments=len(scientific_evidence),
                     active_source_and_input_pins=len(protocol['source_sha256']) + len(protocol['input_sha256'])), indent=2))
