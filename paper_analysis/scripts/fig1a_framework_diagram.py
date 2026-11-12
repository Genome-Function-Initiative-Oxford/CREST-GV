"""Fig. 1a -- CREST-GV enrichment framework (schematic).

Draws the conceptual pipeline described in the Methods and in the Fig. 1
legend: trait-associated variants (N) and background variants (G) are
intersected with the open chromatin regions of each cell type (reference
atlas or a user-provided collection); the observed vs. background overlap
proportion gives a per-cell-type enrichment Z-score, which is standardised
within its data collection and converted to a right-tailed P-value.

This is a schematic, not a data plot -- it carries no numbers of its own.
Sized to the portrait aspect ratio compose_figures.py allocates to panel a
(height equal to panel b's full height, width ~0.35x panel b's width).

Run from the repo root:
    python scripts/fig1a_framework_diagram.py
"""

import os
import sys

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

sys.path.insert(0, os.path.dirname(__file__))
from palette import CATEGORICAL, INK_PRIMARY, INK_SECONDARY, SURFACE, apply_base_style

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG_DIR = os.path.join(ROOT, "figures")

apply_base_style()

FIG_W, FIG_H = 3.6, 5.85  # portrait, matches the box compose_figures.py gives panel a
fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")
fig.patch.set_facecolor(SURFACE)
ax.set_facecolor(SURFACE)


def box(xc, yc, w, h, text, sub=None, facecolor="#ffffff", edgecolor=INK_SECONDARY, fontsize=8.3):
    patch = FancyBboxPatch(
        (xc - w / 2, yc - h / 2), w, h,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        linewidth=1.1, edgecolor=edgecolor, facecolor=facecolor, zorder=2,
    )
    ax.add_patch(patch)
    if sub:
        ax.text(xc, yc + h * 0.14, text, ha="center", va="center",
                 fontsize=fontsize, fontweight="bold", color=INK_PRIMARY, zorder=3)
        ax.text(xc, yc - h * 0.26, sub, ha="center", va="center",
                 fontsize=fontsize - 1.3, color=INK_SECONDARY, zorder=3, wrap=True)
    else:
        ax.text(xc, yc, text, ha="center", va="center",
                 fontsize=fontsize, fontweight="bold", color=INK_PRIMARY, zorder=3)
    return (xc, yc - h / 2), (xc, yc + h / 2)  # (bottom, top) anchor points


def arrow(p0, p1, color=INK_SECONDARY):
    ax.add_patch(FancyArrowPatch(
        p0, p1, arrowstyle="-|>", mutation_scale=9,
        linewidth=1.2, color=color, zorder=1, shrinkA=1, shrinkB=1,
    ))


# panel letter
ax.text(0.03, 0.985, "a", transform=ax.transAxes, fontsize=13, fontweight="bold", va="top", ha="left")

# --- stage 1: two inputs, side by side ---
y1 = 0.90
(n_bottom, n_top) = box(0.265, y1, 0.44, 0.115,
                         r"Trait variants $\mathcal{N}$",
                         "GWAS lead + proxies,\nor user-provided set",
                         facecolor="#eaf1fc", edgecolor=CATEGORICAL["blue"], fontsize=7.8)
(g_bottom, g_top) = box(0.735, y1, 0.44, 0.115,
                         r"Background $\mathcal{G}$",
                         "1000 Genomes\nProject",
                         facecolor="#f3f2ef", edgecolor=INK_SECONDARY, fontsize=7.8)

# --- stage 2: intersect with open chromatin ---
y2 = 0.685
(chrom_bottom, chrom_top) = box(0.5, y2, 0.82, 0.125,
                                 r"Open chromatin $C_{ij}$",
                                 "cell type $i$, collection $j$\n(reference atlas or upload)",
                                 facecolor="#e6f6ef", edgecolor=CATEGORICAL["aqua"])
arrow(n_bottom, (0.5 - 0.12, chrom_top[1]), color=CATEGORICAL["blue"])
arrow(g_bottom, (0.5 + 0.12, chrom_top[1]), color=INK_SECONDARY)

# --- stage 3: enrichment Z-score ---
y3 = 0.47
(z_bottom, z_top) = box(0.5, y3, 0.82, 0.125,
                         r"Enrichment score $Z_{ij}$",
                         "observed vs. background\noverlap proportion",
                         facecolor="#fdf0e8", edgecolor=CATEGORICAL["orange"])
arrow(chrom_bottom, z_top, color=CATEGORICAL["aqua"])

# --- stage 4: standardisation ---
y4 = 0.255
(std_bottom, std_top) = box(0.5, y4, 0.82, 0.125,
                             r"Standardised $\tilde{Z}_{ij}$",
                             "re-scaled within\ndata collection $j$",
                             facecolor="#eeecf8", edgecolor=CATEGORICAL["violet"])
arrow(z_bottom, std_top, color=CATEGORICAL["orange"])

# --- stage 5: P-value ---
y5 = 0.065
(p_bottom, p_top) = box(0.5, y5, 0.6, 0.095,
                         r"$P$-value (right-tailed)",
                         facecolor="#fbeaea", edgecolor=CATEGORICAL["red"], fontsize=9)
arrow(std_bottom, p_top, color=CATEGORICAL["violet"])

out_base = os.path.join(FIG_DIR, "Fig1a_framework_diagram")
fig.savefig(out_base + ".png", dpi=300, bbox_inches="tight", facecolor=SURFACE)
fig.savefig(out_base + ".pdf", bbox_inches="tight", facecolor=SURFACE)
print("wrote", out_base + ".png/.pdf")
