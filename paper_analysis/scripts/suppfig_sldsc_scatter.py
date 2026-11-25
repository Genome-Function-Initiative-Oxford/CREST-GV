"""Supplementary Fig. -- CREST-GV vs S-LDSC concordance on the final atlas.

(a) Spearman rho (S-LDSC --h2-cts coefficient vs standardised CREST-GV
    Z-score), one point per trait, restricted to the 390 QC-kept profiles.
(b) Scatter of S-LDSC coefficient vs standardised CREST-GV Z-score, pooling
    all 4 traits x 390 cell types, coloured by trait.

Run from the repo root (after paper_analysis/sldsc/scripts/compare_with_crestgv.py):
    python scripts/suppfig_sldsc_scatter.py
"""

import os
import sys

import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(__file__))
from palette import apply_base_style, CATEGORICAL, GRIDLINE, INK_MUTED

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H2_CTS_DIR = os.path.join(ROOT, "sldsc", "h2_cts_paper_atlas")

TRAIT_COLORS = {
    "heart_rate": CATEGORICAL["blue"],
    "red_blood_cell_traits": CATEGORICAL["red"],
    "type_1_diabetes_mellitus": CATEGORICAL["aqua"],
    "type_2_diabetes_mellitus": CATEGORICAL["orange"],
}
TRAIT_LABELS = {
    "heart_rate": "Heart rate",
    "red_blood_cell_traits": "Red blood cell traits",
    "type_1_diabetes_mellitus": "Type 1 diabetes",
    "type_2_diabetes_mellitus": "Type 2 diabetes",
}


def main():
    apply_base_style()

    pooled = []
    rhos = {}
    for trait in TRAIT_COLORS:
        df = pd.read_csv(os.path.join(H2_CTS_DIR, f"{trait}.vs_crestgv.tsv"), sep="\t")
        df["trait"] = trait
        pooled.append(df)
        rho, p = spearmanr(df["Coefficient"], df["CREST_GV_norm"])
        rhos[trait] = (rho, p)
    pooled = pd.concat(pooled, ignore_index=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.5))

    # Panel a: per-trait rho
    trait_order = sorted(rhos, key=lambda t: rhos[t][0], reverse=True)
    for i, trait in enumerate(trait_order):
        rho, p = rhos[trait]
        ax1.scatter(i, rho, color=TRAIT_COLORS[trait], s=110, edgecolor="white", linewidth=0.8, zorder=3)
        ax1.text(i, rho + 0.015, f"$P$={p:.1e}", ha="center", va="bottom", fontsize=7, color=INK_MUTED)
    ax1.axhline(0, color=GRIDLINE, linewidth=1)
    ax1.set_xticks(range(len(trait_order)))
    ax1.set_xticklabels([TRAIT_LABELS[t] for t in trait_order], rotation=20, ha="right")
    ax1.set_ylabel("Spearman rho (S-LDSC coefficient vs CREST-GV standardised Z)")
    ax1.set_ylim(0, max(r for r, _ in rhos.values()) * 1.35)
    ax1.set_title("a   Concordance by trait (390 QC-kept cell types)", loc="left", fontsize=11)

    # Panel b: pooled scatter
    for trait, color in TRAIT_COLORS.items():
        sub = pooled[pooled["trait"] == trait]
        ax2.scatter(
            sub["Coefficient"], sub["CREST_GV_norm"],
            color=color, s=16, alpha=0.5, edgecolor="none", label=TRAIT_LABELS[trait],
        )
    ax2.axhline(0, color=GRIDLINE, linewidth=1)
    ax2.axvline(0, color=GRIDLINE, linewidth=1)
    ax2.set_xlabel("S-LDSC per-cell-type coefficient")
    ax2.set_ylabel("Standardised CREST-GV Z-score")
    ax2.set_title("b   Pooled across all 4 traits x 390 cell types", loc="left", fontsize=11)
    ax2.legend(frameon=False, fontsize=8)

    fig.tight_layout()
    out_base = os.path.join(ROOT, "figures", "SuppFig_SLDSC_scatter")
    fig.savefig(out_base + ".png")
    fig.savefig(out_base + ".pdf")
    print("wrote", out_base + ".png/.pdf")


if __name__ == "__main__":
    main()
