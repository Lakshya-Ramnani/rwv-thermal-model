"""
Cross-section schematic through one radial channel of the RWV insert.

Pure illustration, no physics. Produces Figure 1 of the paper.

Author: Lakshya Ramnani
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

# Palette: fill / edge pairs, then accents
DRUG = ('#D4E4F2', '#3D6E94')
TI = ('#C8C5BC', '#5B5750')
F127 = ('#C5BDE0', '#544981')
MEDIUM = ('#FAE3B8', '#9C7220')
TISSUE = ('#C7E2D4', '#3F7355')
FRIT, TEXT, LABEL, HEAT, COLD, DIM = '#3A3A38', '#1F1F1D', '#3D3A36', '#A8362B', '#2E5E8C', '#6B6760'
SERIF = dict(family='serif')

# Geometry, figure units
Y_BOT, Y_TOP = 10, 60
CH_BOT, CH_TOP = 28, 42
CH_MID = (CH_BOT + CH_TOP) / 2
X_DRUG_L, X_TI_L, X_TI_R, X_TIS_R = 5, 28, 80, 102
X_F127_R = X_TI_L + 20
X_FRIT_L, X_FRIT_R = X_F127_R, X_F127_R + 1.5

fig, ax = plt.subplots(figsize=(11, 7), dpi=300)
ax.set_xlim(0, 116)
ax.set_ylim(-12, 78)
ax.set_aspect('equal')
ax.axis('off')


def box(x0, x1, y0, y1, colors, lw=1.0, z=1, **kw):
    fill, edge = colors
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor=fill,
                           edgecolor=edge, linewidth=lw, zorder=z, **kw))


def label(x, y, text, size=11, color=TEXT, weight='normal', style='normal'):
    ax.text(x, y, text, ha='center', va='center', fontsize=size, color=color,
            fontweight=weight, style=style, **SERIF)


def dim(x1, x2, y, text):
    for x in (x1, x2):
        ax.plot([x, x], [y - 0.7, y + 0.7], color=DIM, linewidth=0.7)
    ax.annotate('', xy=(x2, y), xytext=(x1, y),
                arrowprops=dict(arrowstyle='<|-|>', color=DIM, linewidth=0.7, mutation_scale=8))
    ax.text((x1 + x2) / 2, y - 1.6, text, ha='center', va='top',
            fontsize=8.5, color=DIM, **SERIF)


def arrow(x1, y1, x2, y2, color):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle='-|>', mutation_scale=15,
                                 color=color, linewidth=1.8, shrinkA=0, shrinkB=0, zorder=8))


def badge(x, y, n, color, radius=1.8, size=9):
    ax.add_patch(plt.Circle((x, y), radius, facecolor=color, edgecolor='white',
                            linewidth=1.3, zorder=9))
    ax.text(x, y, str(n), ha='center', va='center', fontsize=size,
            color='white', fontweight='bold', zorder=10)


# --- Regions ---------------------------------------------------------------
box(X_DRUG_L, X_TI_L, Y_BOT, Y_TOP, DRUG)
box(X_TI_L, X_TI_R, Y_BOT, Y_TOP, TI, lw=1.2)
box(X_TI_L, X_TI_R, CH_BOT, CH_TOP, ('white', 'none'), z=2)          # clear the bore
box(X_TI_L, X_F127_R, CH_BOT, CH_TOP, F127, z=3)
box(X_FRIT_R, X_TI_R, CH_BOT, CH_TOP, MEDIUM, lw=0.8, z=3, linestyle=(0, (4, 2)))
box(X_TI_R, X_TIS_R, Y_BOT, Y_TOP, TISSUE)

# Sintered frit hatching
for i in range(5):
    y = CH_BOT + (CH_TOP - CH_BOT) * (i + 0.5) / 5
    ax.plot([X_FRIT_L, X_FRIT_R], [y, y], color=FRIT, linewidth=0.7, zorder=4)
for x in (X_FRIT_L, X_FRIT_R):
    ax.plot([x, x], [CH_BOT, CH_TOP], color=FRIT, linewidth=1.2, zorder=4)

# Channel walls
for y in (CH_BOT, CH_TOP):
    ax.plot([X_TI_L, X_TI_R], [y, y], color=TI[1], linewidth=1.4, zorder=5)

# Spheroids, kept clear of the tissue-chamber label
rng = np.random.default_rng(7)
placed = 0
for _ in range(100):                                     # bounded rejection sampling
    if placed == 7:
        break
    cx = rng.uniform(X_TI_R + 3.5, X_TIS_R - 3.5)
    cy = rng.uniform(Y_BOT + 3, Y_TOP - 3)
    if (cx - (X_TI_R + 12)) ** 2 + (cy - 36) ** 2 < 64:
        continue
    ax.add_patch(plt.Circle((cx, cy), rng.uniform(1.3, 2.1), facecolor=TISSUE[1],
                            alpha=0.55, linewidth=0, zorder=2))
    placed += 1

# --- Region labels ---------------------------------------------------------
x_drug_c = (X_DRUG_L + X_TI_L) / 2
label(x_drug_c, 38, 'Drug\nchamber', weight='bold')
label(x_drug_c, 31, '4 $\\degree$C', 9, DRUG[1], style='italic')

x_ti_c = (X_TI_L + X_TI_R) / 2
label(x_ti_c, 53, 'Ti-6Al-4V wall', weight='bold')
label(x_ti_c, 49, '5 mm thick', 8.5, LABEL, style='italic')

label(X_TI_R + 12, 36, 'Tissue\nchamber', weight='bold')
label(X_TI_R + 12, 30, '37 $\\degree$C medium\n+ spheroids', 8.5, TISSUE[1], style='italic')

# --- Leader labels ---------------------------------------------------------
for text, x, y_anchor, y_text, color, size, style in [
    ('F-127 plug', (X_TI_L + X_F127_R) / 2, CH_BOT, 19, F127[1], 10, 'normal'),
    ('sintered Ti frit', (X_FRIT_L + X_FRIT_R) / 2, CH_TOP, 67, FRIT, 9, 'normal'),
    ('warm medium\n(fills channel)', (X_FRIT_R + X_TI_R) / 2, CH_BOT, 19, MEDIUM[1], 9, 'italic'),
]:
    ax.annotate(text, xy=(x, y_anchor), xytext=(x, y_text), ha='center', fontsize=size,
                color=color, style=style, **SERIF,
                arrowprops=dict(arrowstyle='-', color=color, linewidth=0.7))

# --- Dimensions ------------------------------------------------------------
dim(X_TI_L, X_F127_R, 8, '2 mm')
dim(X_FRIT_R, X_TI_R, 8, '3 mm')

x_dim = X_TIS_R + 3
for y in (CH_BOT, CH_TOP):
    ax.plot([X_TIS_R + 0.5, x_dim + 0.5], [y, y], color=DIM, linewidth=0.5, alpha=0.6)
ax.annotate('', xy=(x_dim, CH_TOP), xytext=(x_dim, CH_BOT),
            arrowprops=dict(arrowstyle='<|-|>', color=DIM, linewidth=0.7, mutation_scale=8))
ax.text(x_dim + 1, CH_MID, '2 mm\nchannel\ndiameter', fontsize=7.5, color=DIM,
        va='center', ha='left', **SERIF)

# --- Heat paths ------------------------------------------------------------
arrow(33, 47, 33, CH_TOP + 0.5, HEAT)                    # 1, radial in from above
arrow(43, 23, 43, CH_BOT - 0.5, HEAT)                    # 1, radial in from below
badge(33, 47.5, 1, HEAT)
arrow(60, CH_MID, X_FRIT_R + 0.3, CH_MID, HEAT)          # 2, axial through the frit
badge(64, CH_MID, 2, HEAT)
arrow(X_TI_L + 5, CH_MID, X_TI_L + 0.3, CH_MID, COLD)    # 3, axial loss to the drug
badge(X_TI_L + 6.5, CH_MID, 3, COLD)

# --- Legend ----------------------------------------------------------------
ax.text(5, -2.5, 'Heat transfer paths:', fontsize=10, color=TEXT, fontweight='bold', **SERIF)
for i, (n, color, text) in enumerate([
    (1, HEAT, 'Cylindrical Ti walls heat F-127 radially  (dominant path)'),
    (2, HEAT, 'Warm medium heats F-127 axially through the frit'),
    (3, COLD, 'F-127 loses heat axially to cold drug'),
]):
    y = -5.2 - 3 * i
    badge(7, y, n, color, radius=1.4, size=8)
    ax.text(9.6, y, text, ha='left', va='center', fontsize=9, color=TEXT, **SERIF)

# --- Orientation and title -------------------------------------------------
ax.annotate('', xy=(115, 67), xytext=(108, 67),
            arrowprops=dict(arrowstyle='->', color=LABEL, linewidth=1.0))
ax.text(111.5, 64.5, 'z (axial)', fontsize=8.5, color=LABEL, va='center',
        ha='center', style='italic', **SERIF)
ax.text(58, 73, 'Cross-section through one radial channel', ha='center', va='bottom',
        fontsize=12.5, color=TEXT, fontweight='bold', **SERIF)

fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(f'figure_cross_section.{ext}', dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print("Saved figure_cross_section.png and .pdf")
