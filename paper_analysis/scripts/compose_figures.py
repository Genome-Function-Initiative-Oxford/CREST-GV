"""Assemble per-panel images into the two composite main-text figures.

Per Simone: a figure's lettered panels should sit together on one composited
page/image, not just be labelled as parts of the same figure across separate
pages. This script builds those composites from the already-generated panel
images -- it doesn't recompute anything, just arranges existing PNGs.

Fig. 1a-d -> Fig1_abcd_composite.{png,pdf}:
    a, Fig1a_framework_diagram.png (schematic built by
       fig1a_framework_diagram.py; panel letter baked in).
    b, WEB_crestgv_overview.png (raw web-app screenshot montage, placeholder
       -- labelled 'b' here since the source image has no panel label of its
       own).
    c-d, Fig2_standardisation.png (already a single image with both panels
         labelled and stacked by fig2_standardisation.py).

Fig. 2 (the cross-trait heatmap, Fig3_cross_trait_heatmap.{png,pdf}) is its
own standalone figure, built by fig3_cross_trait_heatmap.py -- not composited
here, per Simone: its tall, many-row layout doesn't sit well alongside other
panels, and it's its own figure number now rather than a panel of Fig. 1.

Fig. 3a-f -> Fig3_abcdef_composite.{png,pdf}: stacks Fig4_per_trait_rankings.png
(a-d, already labelled) above Fig5_upload_example.png (e-f, already labelled).

Run from the repo root:
    python scripts/compose_figures.py
"""

import os
import sys

import matplotlib.image as mpimg
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG_DIR = os.path.join(ROOT, "figures")


def stack_rows(images, out_base, width_in=10):
    """Stack full-width image rows vertically, preserving each row's own aspect ratio."""
    widths = [im.shape[1] for im in images]
    heights = [im.shape[0] for im in images]
    W = max(widths)
    heights_scaled = [h * (W / w) for h, w in zip(heights, widths)]

    fig = plt.figure(figsize=(width_in, width_in * sum(heights_scaled) / W))
    gs = fig.add_gridspec(len(images), 1, height_ratios=heights_scaled, hspace=0.04)
    for i, im in enumerate(images):
        ax = fig.add_subplot(gs[i])
        ax.imshow(im)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis("off")

    fig.savefig(out_base + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(out_base + ".pdf", bbox_inches="tight")
    print("wrote", out_base + ".png/.pdf")


def main():
    # --- Fig. 1a-d ---
    a_img = mpimg.imread(os.path.join(FIG_DIR, "Fig1a_framework_diagram.png"))
    b_img = mpimg.imread(os.path.join(FIG_DIR, "WEB_crestgv_overview.png"))
    cd_img = mpimg.imread(os.path.join(FIG_DIR, "Fig2_standardisation.png"))

    b_h, b_w = b_img.shape[0], b_img.shape[1]
    a_h, a_w = a_img.shape[0], a_img.shape[1]
    # scale panel a to panel b's height, then express its width as a fraction of b_w
    a_w_frac = (a_w * (b_h / a_h)) / b_w
    row1_w = b_w * (1 + a_w_frac)
    row1_h = b_h
    row2_w, row2_h = cd_img.shape[1], cd_img.shape[0]

    W = max(row1_w, row2_w)
    row1_h_scaled = row1_h * (W / row1_w)
    row2_h_scaled = row2_h * (W / row2_w)

    fig = plt.figure(figsize=(10, 10 * (row1_h_scaled + row2_h_scaled) / W))
    gs = fig.add_gridspec(2, 1, height_ratios=[row1_h_scaled, row2_h_scaled], hspace=0.04)

    gs_row1 = gs[0].subgridspec(1, 2, width_ratios=[a_w_frac, 1], wspace=0.02)
    ax_a = fig.add_subplot(gs_row1[0, 0])
    ax_b = fig.add_subplot(gs_row1[0, 1])
    ax_cd = fig.add_subplot(gs[1])

    ax_a.imshow(a_img, aspect="auto")
    ax_a.set_xticks([])
    ax_a.set_yticks([])
    ax_a.axis("off")

    ax_b.imshow(b_img)
    ax_b.text(0.02, 0.97, "b", transform=ax_b.transAxes, fontsize=13, fontweight="bold",
              va="top", ha="left", color="white",
              bbox=dict(boxstyle="square,pad=0.15", facecolor="black", alpha=0.55, linewidth=0))
    ax_b.set_xticks([])
    ax_b.set_yticks([])
    for spine in ax_b.spines.values():
        spine.set_visible(False)

    ax_cd.imshow(cd_img)
    ax_cd.set_xticks([])
    ax_cd.set_yticks([])
    ax_cd.axis("off")

    out_base = os.path.join(FIG_DIR, "Fig1_abcd_composite")
    fig.savefig(out_base + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(out_base + ".pdf", bbox_inches="tight")
    print("wrote", out_base + ".png/.pdf")

    # --- Fig. 2a-f ---
    stack_rows(
        [
            mpimg.imread(os.path.join(FIG_DIR, "Fig4_per_trait_rankings.png")),
            mpimg.imread(os.path.join(FIG_DIR, "Fig5_upload_example.png")),
        ],
        os.path.join(FIG_DIR, "Fig3_abcdef_composite"),
    )


if __name__ == "__main__":
    main()
