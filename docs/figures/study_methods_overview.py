"""Draw the study's method schematic without reading trajectories or using GPU.

Run with the project Python environment. The diagram is a schematic, not data.
"""
from pathlib import Path
import hashlib
import json
import os

os.environ.setdefault('MPLCONFIGDIR', '/tmp/embryogenesis-methods-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, PathPatch
from matplotlib.path import Path as MplPath
from matplotlib.colors import to_rgb
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
IMAGES = ROOT/'docs/images'
PDFS = ROOT/'output/pdf'
IMAGES.mkdir(parents=True, exist_ok=True)
PDFS.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'mathtext.fontset': 'dejavusans',
                     'svg.fonttype': 'none', 'pdf.fonttype': 42,
                     'savefig.facecolor': 'white'})
W, H = 1700, 1840
fig = plt.figure(figsize=(17, 18.4), facecolor='white')
ax = fig.add_axes([0, 0, 1, 1])
ax.set(xlim=(0, W), ylim=(H, 0), aspect='equal')
ax.axis('off')
INK, MUTED, BORDER = '#162d3c', '#49606d', '#cbd6dc'
BLUE, TEAL, PURPLE, ORANGE = '#2069a1', '#137c76', '#78519e', '#b76528'
TEXTS = []

def txt(x, y, s, size=13.2, color=INK, weight='normal', ha='left', va='top', parent=None):
    artist = ax.text(x, y, s, fontsize=size, color=color, weight=weight,
                     ha=ha, va=va, zorder=6)
    TEXTS.append((artist, parent))
    return artist

def tint(color, strength=.07):
    return tuple((1-strength)+strength*np.asarray(to_rgb(color)))

def box(x, y, w, h, title, color=BLUE):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0,rounding_size=9',
                              facecolor=tint(color), edgecolor=BORDER, linewidth=.8, zorder=2))
    ax.plot([x+14, x+14], [y+16, y+h-16], color=color, lw=2.5, zorder=3)
    bounds = (x+24, y+8, x+w-13, y+h-8)
    txt(x+29, y+15, title, size=15.1, weight='bold', parent=bounds)
    return bounds

def lines(x, y, strings, size=13.2, gap=25, parent=None, color=INK):
    for i, string in enumerate(strings):
        txt(x, y+gap*i, string, size=size, parent=parent, color=color)

def section(letter, y, title, subtitle=None):
    txt(45, y, letter, 19, weight='bold', color=BLUE)
    txt(83, y, title, 18.2, weight='bold')
    if subtitle:
        txt(1655, y+3, subtitle, 12, color=MUTED, ha='right')

def arrow(start, end, color=MUTED, lw=1.5, rad=0, style='-|>', scale=15):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle=style, mutation_scale=scale,
                                connectionstyle=f'arc3,rad={rad}', color=color,
                                linewidth=lw, shrinkA=0, shrinkB=0, zorder=4))

def elbow(points, color=MUTED, lw=1.4):
    vertices = np.asarray(points)
    ax.add_patch(PathPatch(MplPath(vertices, [MplPath.MOVETO]+[MplPath.LINETO]*(len(points)-1)),
                           fill=False, edgecolor=color, linewidth=lw, zorder=3))
    arrow(points[-2], points[-1], color=color, lw=lw)

def sphere(x, y, radius, color='#b2d1df'):
    ax.add_patch(Circle((x, y), radius, fc=color, ec='#65889a', lw=.65, zorder=3))
    for k in range(12):
        r = radius*(.77-.035*k)
        c = np.asarray(to_rgb(color))*(1-.025*k)+np.ones(3)*.025*k
        ax.add_patch(Circle((x-radius*.17, y-radius*.18), r, fc=c, ec='none', zorder=3))
    ax.add_patch(Circle((x-radius*.27, y-radius*.3), radius*.14, fc='white', ec='none', alpha=.55, zorder=3))

txt(45, 30, 'Studying chemical organization in a changing 3D embryo', 27.5, weight='bold')
txt(45, 83, 'Methods overview | Activator-inhibitor feedback, conservative communication and polarity-dependent mechanics', 14, color=MUTED)
ax.plot([45, 1655], [119, 119], color=BORDER, lw=1)

# A. Biological representation and supplied developmental rules.
section('A', 136, 'Founder, cleavage and continuous cell attributes')
b = box(45, 171, 1610, 170, 'One zygote develops into individually deformable cells', TEAL)
centers = [109, 283, 451, 620, 789]
for n, cx in zip([1, 2, 4, 8, 16], centers):
    radius = 29/n**(1/3)
    if n == 1:
        points = [(cx, 251)]
    elif n == 2:
        points = [(cx-22, 251), (cx+22, 251)]
    elif n == 4:
        points = [(cx+dx, 251+dy) for dx in (-18, 18) for dy in (-18, 18)]
    elif n == 8:
        points = [(cx+dx+layer*10, 251+dy-layer*10) for layer in (-.5, .5) for dx in (-18, 18) for dy in (-18, 18)]
    else:
        points = [(cx+dx, 239+dy) for dx in (-33, -11, 11, 33) for dy in (-33, -11, 11, 33)]
    for x, y in sorted(points, key=lambda p:p[1]):
        sphere(x, y, radius)
    txt(cx, 290, f'{n} cell'+('s' if n>1 else ''), 12.8, ha='center')
for x0, x1 in zip(centers[:-1], centers[1:]):
    arrow((x0+58, 251), (x1-58, 251), color=TEAL)
lines(930, 215, [r'Cell state: $\phi_i(\mathbf{x},t),\ a_i,\ b_i,\ \mathbf{p}_i$',
                r'Founder chemistry: $a=b=1$; polarity: $\mathbf{p}=0$',
                'Shape attributes: axis ratio and asphericity',
                'Context: position, lineage, volumes and contacts'], gap=24, parent=b)
txt(74, 316, 'Cleavage: longest-axis spindle; gradual ring constriction; neck/overlap-gated abscission; amount-balanced inheritance.', 12.7, parent=b)

# B. Actual dependency graph: primary loop and the two polarity inputs.
section('B', 365, 'Coupled model: geometry changes communication, and activity changes mechanics')
X = [45, 595, 1145]
top, bottom, cw, ch = 401, 655, 510, 185
g = box(X[0], top, cw, ch, 'Cell geometry and contact extraction', TEAL)
lines(74, 448, [r'$h(\phi)=\phi^2(3-2\phi),\quad q(\phi)=\phi^2(1-\phi)^2$',
                r'$V_i=\int h(\phi_i)\,dV$',
                'Diffuse shell overlaps estimate contact areas.',
                '3D voxels resolve shape; each cell is well mixed.'], size=13.3, gap=28, parent=g)
t = box(X[1], top, cw, ch, 'Conservative chemical transport', TEAL)
lines(624, 449, [r'$g_{ij}=\widehat{\mathcal{A}}_{ij}/\ell_{ij},\quad M=\mathrm{diag}(V_i)$',
                r'$K=\mathrm{diag}(G\mathbf{1})-G,\quad \Delta_V=-M^{-1}K$',
                r'$\sum_i V_i(\Delta_V c)_i=0$',
                'Area/distance closure needs geometric validation.'], size=13.3, gap=28, parent=t)
c = box(X[2], top, cw, ch, 'Activator-inhibitor reactions', BLUE)
lines(1174, 447, [r'$\dot a_i=a_i^2/b_i-a_i+D_a(\Delta_Va)_i-a_i\dot V_i/V_i$',
                 r'$\dot b_i=\beta(a_i^2-b_i)+D_b(\Delta_Vb)_i-b_i\dot V_i/V_i$'], size=13.1, gap=33, parent=c)
for x, label, color in [(1232, 'a', BLUE), (1570, 'b', PURPLE)]:
    ax.add_patch(Circle((x, 548), 18, fc=tint(color,.14), ec=color, lw=1, zorder=4))
    txt(x, 548, label, 15, ha='center', va='center')
arrow((1253, 538), (1547, 538), color=BLUE)
txt(1400, 517, 'induces inhibitor', 11.5, ha='center')
arrow((1548, 563), (1254, 563), color=PURPLE, style='-[,widthB=0.45,lengthB=0')
txt(1400, 567, 'inhibits activator', 11.5, ha='center')
txt(1180, 511, 'self-activates', 11.2)
arrow((1214, 541), (1220, 565), color=BLUE, rad=.95, scale=11)

pol = box(X[0], bottom, cw, ch, 'Apical-basal polarity', PURPLE)
lines(74, 701, [r'$\dot{\mathbf{p}}_i=\alpha\frac{2a_i}{1+a_i}\mathbf{q}_i+\eta(\Delta_{\rm rw}\mathbf{p})_i$',
                r'$\qquad-(\mu+\|\mathbf{p}_i\|^2)\mathbf{p}_i$',
                'Exposure cue points toward free cortex.',
                r'$\Delta_{\rm rw}$ averages orientation only.'], size=13.1, gap=28, parent=pol)
mech = box(X[1], bottom, cw, ch, 'Shape mechanics and volume dilution', ORANGE)
lines(624, 701, [r'$E=E_{\rm interface}+E_{\rm volume}+E_{\rm overlap}-E_{\rm adhesion}$',
                 r'$\partial_t\phi_i=-\delta E/\delta\phi_i$',
                 'Overdamped force; coefficients frozen per step.',
                 r'After motion: $c_i\leftarrow c_i V_i^{\rm old}/V_i^{\rm new}$'], size=13.1, gap=28, parent=mech)
mat = box(X[2], bottom, cw, ch, 'Continuous material response', ORANGE)
lines(1174, 701, [r'$r_i=\tanh(a_i-1),\quad \gamma_i^0=\gamma_0(1+c_\gamma r_i)$',
                 r'$A_{ij}=A_0(1+c_A r_i r_j)$',
                 r'$\gamma_i(\mathbf{x})=\gamma_i^0[1-\chi\mathbf{p}_i\!\cdot\!\widehat{\mathbf{r}}_i(\mathbf{x})]$',
                 'Activity and polarity alter tension / adhesion.'], size=13.1, gap=28, parent=mat)
arrow((555, 491), (595, 491), TEAL)
arrow((1105, 491), (1145, 491), TEAL)
arrow((1460, 586), (1460, 655), BLUE)
txt(1478, 610, 'activity', 12, color=MUTED)
arrow((153, 586), (153, 655), PURPLE)
txt(175, 629, 'exposure cue', 11.8, color=MUTED)
elbow([(1300, 586), (1300, 616), (445, 616), (445, 655)], BLUE)
txt(895, 591, 'activity modulates polarity', 11.7, ha='center', color=MUTED)
arrow((555, 747), (595, 747), PURPLE)
arrow((1145, 747), (1105, 747), ORANGE)
elbow([(840, 840), (840, 864), (29, 864), (29, 491), (45, 491)], ORANGE)
txt(433, 870, 'Updated shape and volume feed the next step', 12.2, ha='center', color=MUTED)
txt(955, 870, r'Uniform chemistry still permits an exposure-driven polarity cue.', 12.2, color=MUTED)
txt(45, 901, r'Frozen-mode screen: $S=M^{-1/2}KM^{-1/2}\ \longrightarrow\ \lambda_k$; growth from $\max\,\mathrm{Re}\,\mathrm{eig}[J-\lambda_k\,\mathrm{diag}(D_a,D_b)]$.', 13.1)

# C. Six families of comparisons. No observed outcome is drawn.
section('C', 942, 'Experiments distinguish formation, persistence and context-sensitive behavior')
items = [
    ('Formation versus maintenance', BLUE, [
        'Near-uniform chemistry versus developed patterns',
        'Feedback / tension / adhesion / polarity controls',
        'Frozen spectra and endpoint coexistence assays']),
    ('Exchange and network allocation', TEAL, [
        'Conservatively exchange chemistry; retain geometry',
        'Vary surrounding allocation and initial dosage',
        'Match recipient chemistry + initial species amounts']),
    ('Perturbation-response behavior', BLUE, [
        'Pulse one cell; subtract its own continuing control',
        'Signed responses, response area and recovery',
        'Fractional and common removed-amount challenges']),
    ('Environmental dependence', PURPLE, [
        'Isolate cells or clamp a shared chemical reservoir',
        'Then allow a finite shared reservoir to evolve',
        'Separate frozen assays; no live extracellular field']),
    ('Robustness and mechanism controls', ORANGE, [
        'Diffusivity ratios and polarity-tension contrasts',
        'Frozen screens, then selected moving continuations',
        'Preserve conductances or remove dilution']),
    ('Three separate measured outcomes', TEAL, [
        'Cell attributes: differences, persistence and recovery',
        'Signaling: modes and chemical response profiles',
        'Geometry: axis ratio, asphericity and axis persistence']),
]
for index, (title, color, body) in enumerate(items):
    x, y = X[index%3], 978+(index//3)*154
    bounds = box(x, y, cw, 137, title, color)
    lines(x+29, y+51, body, size=12.45, gap=25, parent=bounds)

# D. Current prospective protocol, explicitly pending.
section('D', 1294, 'Current moving confirmation: match once, let the tissue evolve, then challenge', 'ONGOING')
timeline = [
    (45, 415, 't = 510 | Initial matching', TEAL, [
        'Identical recipient activator / inhibitor',
        'Equal initial surrounding species amounts',
        'Shared initial sham geometry and polarity',
        '6 recipients x 2 chemical allocations']),
    (505, 500, 't = 510-570 | Free reorganization', ORANGE, [
        'Chemistry, transport, polarity and mechanics evolve',
        'Volumes and amounts may diverge after matching',
        '12 backgrounds x two timesteps = 24 paths',
        'Use each actual fine endpoint for the next stage']),
    (1050, 605, 't = 570-630 | Delayed response', BLUE, [
        'Own unpulsed control, -10% pulse, common-amount pulse',
        'Same actual fine start at both timesteps; no rematching',
        '12 backgrounds x three arms x two timesteps = 72 paths',
        'Compare signed two-species response waveforms']),
]
for x, w, title, color, body in timeline:
    bounds = box(x, 1331, w, 157, title, color)
    lines(x+29, 1380, body, size=12.25, gap=25, parent=bounds)
arrow((465, 1410), (500, 1410))
arrow((1010, 1410), (1045, 1410))
txt(45, 1498, r'96 new continuations; $\Delta t=0.00375,\ 0.001875$; three existing histories (7, 8, 9). Recipients, doses and timesteps are nested tests.', 12.45)

# E. Computational and evidential gates.
section('E', 1540, 'Numerical validation and evidence handling')
checks = [
    ('Equations and conservation', TEAL, [
        'Amount-balanced exchange / dilution',
        'Positive chemical integration',
        'Frozen graph stability and recovery']),
    ('Precision and integration', BLUE, [
        'SSP-RK2 chemical substeps',
        'Phase carry + chemical compensation',
        'First-order full-system splitting']),
    ('Independent comparisons', PURPLE, [
        'Actual-start native CPU / GPU checks',
        'DOP853 and Radau frozen references',
        'Timestep, response and restart gates']),
    ('Reproducibility and quality', ORANGE, [
        'Pinned protocols, sources and inputs',
        'Full-state restart with corrections / RNG',
        'Volume, radius, clipping, boundary gates']),
]
for index, (title, color, body) in enumerate(checks):
    x = 45+410*index
    bounds = box(x, 1578, 380, 132, title, color)
    lines(x+29, 1628, body, size=11.7, gap=24, parent=bounds)
txt(45, 1723, 'CPU: development / cleavage and native references. GPU: resident float64 PyTorch chemistry + custom CUDA mechanics / geometry / polarity.', 12.2)
ax.plot([45, 1655], [1756, 1756], color=BORDER, lw=.9)
lines(45, 1768, [
    'Scope: dimensionless, well-mixed regulators; permissive nutrients / energy; no biomass growth between divisions; no cavity mechanics.',
    'Geometry closure, full developmental convergence, cell autonomy and inheritance remain separate requirements.',
    'Schematic cells and arrows, not simulation output. Two chemical species do not define two cell identities.'
], size=12.1, gap=21, color=MUTED)

# Detect text outside its assigned panel or the canvas before publishing.
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
outside = []
for artist, parent in TEXTS:
    bounds = artist.get_window_extent(renderer).transformed(ax.transData.inverted())
    x0, x1 = sorted((bounds.x0, bounds.x1))
    y0, y1 = sorted((bounds.y0, bounds.y1))
    limit = parent or (15, 12, W-15, H-12)
    if x0 < limit[0]-1 or y0 < limit[1]-1 or x1 > limit[2]+1 or y1 > limit[3]+1:
        outside.append(dict(text=artist.get_text(), actual=[x0,y0,x1,y1], allowed=limit))
if outside:
    raise RuntimeError(json.dumps(outside, indent=2))

stem = 'study-methods-overview'
fig.savefig(IMAGES/(stem+'.svg'))
fig.savefig(IMAGES/(stem+'.png'), dpi=240)
fig.savefig(PDFS/(stem+'.pdf'), metadata={'Title': 'Methods for chemical organization in a changing 3D embryo',
                                        'Subject': 'Schematic of the model, assays, numerical validation and interpretation limits'})
preview = ROOT/'tmp/pdfs/study-methods-overview-preview.png'
preview.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(preview, dpi=100)
sources = ['README.md','docs/model.md','docs/attribute_development.md','docs/live_transport.md',
           'docs/graph_signaling.md','docs/cytokinesis.md','docs/compensated_chemistry.md',
           'docs/moving_recipient_amount_response.md','docs/phase_carry_polarity_assessment.md',
           'docs/attribute_common_environment.md','docs/attribute_finite_reservoir.md']
sha = lambda file: hashlib.sha256(Path(file).read_bytes()).hexdigest()
record = dict(schematic=True, observed_results_inferred=False, gpu_used=False,
              source_sha256={name:sha(ROOT/name) for name in sources},
              generator_sha256=sha(__file__), text_inside_panels=True,
              outputs={str(file.relative_to(ROOT)):sha(file) for file in
                       [IMAGES/(stem+'.svg'), IMAGES/(stem+'.png'), PDFS/(stem+'.pdf')]})
(IMAGES/(stem+'.json')).write_text(json.dumps(record, indent=2)+'\n')
print(json.dumps({'outputs': list(record['outputs']), 'text_inside_panels':True,'preview':str(preview)}, indent=2))
