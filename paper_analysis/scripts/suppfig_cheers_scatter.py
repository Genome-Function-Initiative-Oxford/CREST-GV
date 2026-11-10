"""Supplementary Fig. -- CREST-GV vs CHEERS concordance on the final atlas.

(a) Per-collection x per-trait Spearman rho (CHEERS p-value vs CREST-GV's
    standardised p-value), one point per of the 44 collection x trait pairs,
    grouped by collection.
(b) Scatter of -log10(CHEERS p-value) vs -log10(CREST-GV standardised
    p-value) pooling all pairs, coloured by trait.

Run from the repo root (after cheers_benchmark_compare.py):
    python scripts/suppfig_cheers_scatter.py
"""

import glob
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from palette import apply_base_style, CATEGORICAL, GRIDLINE, INK_MUTED

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BENCHMARK_DIR = os.path.join(ROOT, "cheers_benchmark")

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

    summary = pd.read_csv(os.path.join(BENCHMARK_DIR, "summary.tsv"), sep="\t")
    collection_order = (
        summary.groupby("collection")["spearman_rho_vs_norm"].median().sort_values(ascending=False).index.tolist()
    )

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Panel a: per-collection rho, points coloured by trait
    for i, coll in enumerate(collection_order):
        sub = summary[summary["collection"] == coll]
        for _, row in sub.iterrows():
            ax1.scatter(i, row["spearman_rho_vs_norm"], color=TRAIT_COLORS[row["trait"]], s=70,
                        edgecolor="white", linewidth=0.6, zorder=3)
        med = sub["spearman_rho_vs_norm"].median()
        ax1.plot([i - 0.3, i + 0.3], [med, med], color=INK_MUTED, linewidth=2, zorder=2)

    ax1.axhline(0, color=GRIDLINE, linewidth=1)
    ax1.set_xticks(range(len(collection_order)))
    ax1.set_xticklabels(collection_order, rotation=40, ha="right")
    ax1.set_ylabel("Spearman rho (CHEERS vs CREST-GV standardised p-value)")
    ax1.set_title("a   Concordance by collection (horizontal bar = median)", loc="left", fontsize=11)

    handles = [plt.Line2D([0], [0], marker="o", linestyle="", color=c, markersize=8, label=TRAIT_LABELS[t])
               for t, c in TRAIT_COLORS.items()]
    ax1.legend(handles=handles, loc="lower left", frameon=False, fontsize=8)

    # Panel b: pooled -log10(p) scatter across all 44 pairs
    all_pairs = []
    for f in sorted(glob.glob(os.path.join(BENCHMARK_DIR, "per_pair", "*.tsv"))):
        trait = "_".join(os.path.basename(f).replace(".tsv", "").split("__")[1:])
        df = pd.read_csv(f, sep="\t")
        df["trait"] = trait
        all_pairs.append(df)
    pooled = pd.concat(all_pairs, ignore_index=True)
    pooled = pooled[(pooled["CHEERS_pval"] > 0) & (pooled["CREST-GV_pval_norm"] > 0)]

    for trait, color in TRAIT_COLORS.items():
        sub = pooled[pooled["trait"] == trait]
        ax2.scatter(
            -np.log10(sub["CHEERS_pval"]), -np.log10(sub["CREST-GV_pval_norm"]),
            color=color, s=18, alpha=0.5, edgecolor="none", label=TRAIT_LABELS[trait],
        )
    ax2.set_xlabel("-log10(CHEERS P)")
    ax2.set_ylabel("-log10(CREST-GV standardised P)")
    ax2.set_title("b   Pooled across all 44 collection x trait pairs", loc="left", fontsize=11)
    ax2.legend(frameon=False, fontsize=8)

    fig.tight_layout()
    out_base = os.path.join(ROOT, "figures", "SuppFig_CHEERS_scatter")
    fig.savefig(out_base + ".png")
    fig.savefig(out_base + ".pdf")
    print("wrote", out_base + ".png/.pdf")


if __name__ == "__main__":
    main()
