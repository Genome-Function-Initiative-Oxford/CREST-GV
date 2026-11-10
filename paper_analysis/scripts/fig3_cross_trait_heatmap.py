"""Fig. 3 -- CREST-GV enrichment across four GWAS traits.

Standardised scores (Z-tilde) for the cell types significant (Benjamini-
Hochberg FDR q < 0.10, within each trait, across the 390 kept cell types) in
at least one of: heart rate, red blood cell traits, type 1 diabetes, type 2
diabetes. Rows grouped by the trait in which each cell type is significant.

Run from the repo root:
    python scripts/fig3_cross_trait_heatmap.py
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm

sys.path.insert(0, os.path.dirname(__file__))
from palette import apply_base_style, DIVERGING_NEG, DIVERGING_POS, DIVERGING_MID, INK_MUTED
from stats_utils import bh_fdr_mask

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCORES_FILE = os.path.join(ROOT, "data", "scores_390kept.tsv")
PVALS_FILE = os.path.join(ROOT, "data", "pvals_390kept.tsv")

QTHR = 0.10
TRAITS = ["heart_rate", "red_blood_cell_traits", "type_1_diabetes_mellitus", "type_2_diabetes_mellitus"]
TRAIT_LABELS = ["Heart rate", "Red blood cell traits", "Type 1 diabetes", "Type 2 diabetes"]
TRAIT_COLS = [f"{t}_score" for t in TRAITS]


def main():
    apply_base_style()

    scores = pd.read_csv(SCORES_FILE, sep="\t", index_col=0)[TRAIT_COLS]
    pvals = pd.read_csv(PVALS_FILE, sep="\t", index_col=0)[TRAIT_COLS]
    scores.columns = TRAIT_LABELS
    pvals.columns = TRAIT_LABELS

    is_sig = pd.DataFrame({t: bh_fdr_mask(pvals[t], QTHR) for t in TRAIT_LABELS})
    sig_in_any = is_sig.any(axis=1)
    cell_types = scores.index[sig_in_any]
    print(f"{len(cell_types)} cell types significant (BH-FDR q<{QTHR}) in >=1 of the 4 traits")

    sub = scores.loc[cell_types]
    sub_sig = is_sig.loc[cell_types]

    # Group rows by the trait each cell type is "most" significant for: its
    # primary trait is the one with the lowest P among traits it's significant in.
    primary_trait = sub_sig.apply(
        lambda row: pvals.loc[row.name, TRAIT_LABELS][row[row].index].idxmin() if row.any() else None,
        axis=1,
    )
    order = (
        pd.DataFrame({"primary": primary_trait, "score": sub.max(axis=1)})
        .assign(trait_rank=lambda d: d["primary"].map({t: i for i, t in enumerate(TRAIT_LABELS)}))
        .sort_values(["trait_rank", "score"], ascending=[True, False])
        .index
    )
    sub = sub.loc[order]

    fig_h = max(6, 0.22 * len(sub))
    fig, ax = plt.subplots(figsize=(7, fig_h))

    vmax = np.nanmax(np.abs(sub.values))
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    cmap = plt.cm.colors.LinearSegmentedColormap.from_list(
        "diverging", [DIVERGING_NEG, DIVERGING_MID, DIVERGING_POS]
    )

    im = ax.imshow(sub.values, aspect="auto", cmap=cmap, norm=norm)
    ax.set_xticks(range(len(TRAIT_LABELS)))
    ax.set_xticklabels(TRAIT_LABELS, rotation=30, ha="right")
    ax.set_yticks(range(len(sub)))
    ax.set_yticklabels(sub.index, fontsize=6)
    ax.grid(False)

    cbar = fig.colorbar(im, ax=ax, shrink=0.5, pad=0.02)
    cbar.set_label(r"Standardised score ($\tilde{Z}$)", color=INK_MUTED)

    ax.set_title(f"Cell types significant at BH-FDR q < {QTHR} in >=1 trait", loc="left", fontsize=11)
    fig.tight_layout()

    out_base = os.path.join(ROOT, "figures", "Fig3_cross_trait_heatmap")
    fig.savefig(out_base + ".png")
    fig.savefig(out_base + ".pdf")
    print("wrote", out_base + ".png/.pdf")

    sub.to_csv(os.path.join(ROOT, "tables", "Fig3_heatmap_values.csv"))


if __name__ == "__main__":
    main()
