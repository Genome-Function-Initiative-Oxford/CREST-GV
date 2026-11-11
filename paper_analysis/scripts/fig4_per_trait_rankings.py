"""Fig. 3a-d -- Cell type rankings for the four GWAS traits.

Top and bottom 10 standardised enrichment scores (Z-tilde) across all
collections for (a) heart rate, (b) red blood cell traits, (c) type 1
diabetes, (d) type 2 diabetes. Bars coloured by sign; asterisk marks
Benjamini-Hochberg FDR q < 0.10 (computed over all 390 kept cell types per
trait, not just the top/bottom 20 shown).

Run from the repo root:
    python scripts/fig4_per_trait_rankings.py
"""

import os
import sys

import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from palette import apply_base_style, DIVERGING_POS, DIVERGING_NEG, GRIDLINE
from stats_utils import bh_fdr_mask

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCORES_FILE = os.path.join(ROOT, "data", "scores_390kept.tsv")
PVALS_FILE = os.path.join(ROOT, "data", "pvals_390kept.tsv")

QTHR = 0.10
TRAITS = [
    ("heart_rate", "a   Heart rate"),
    ("red_blood_cell_traits", "b   Red blood cell traits"),
    ("type_1_diabetes_mellitus", "c   Type 1 diabetes"),
    ("type_2_diabetes_mellitus", "d   Type 2 diabetes"),
]


def main():
    apply_base_style()

    scores = pd.read_csv(SCORES_FILE, sep="\t", index_col=0)
    pvals = pd.read_csv(PVALS_FILE, sep="\t", index_col=0)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    all_tables = {}

    for ax, (trait, label) in zip(axes.flat, TRAITS):
        col = f"{trait}_score"
        s = scores[col].sort_values(ascending=False)
        top = s.head(10)
        bottom = s.tail(10).sort_values(ascending=False)
        combined = pd.concat([top, bottom])
        sig_full = bh_fdr_mask(pvals[col], QTHR)  # BH-FDR over all 390 cell types, not just the displayed 20
        sig = sig_full.loc[combined.index]

        colors = [DIVERGING_POS if v >= 0 else DIVERGING_NEG for v in combined.values]
        y_pos = range(len(combined))[::-1]
        ax.barh(y_pos, combined.values, color=colors, height=0.65)
        ax.axvline(0, color=GRIDLINE, linewidth=1)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(combined.index, fontsize=7)
        ax.set_title(label, loc="left", fontsize=11)
        ax.set_xlabel(r"$\tilde{Z}$")

        for y, v, is_sig in zip(y_pos, combined.values, sig.values):
            if is_sig:
                offset = 0.05 * (abs(combined.values).max())
                x = v + offset if v >= 0 else v - offset
                ax.text(x, y, "*", ha="left" if v >= 0 else "right", va="center", fontsize=12)

        all_tables[trait] = combined.to_frame(col).assign(significant=sig.values)

    fig.suptitle(f"Top/bottom 10 standardised scores per trait (* = BH-FDR q < {QTHR})", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    out_base = os.path.join(ROOT, "figures", "Fig4_per_trait_rankings")
    fig.savefig(out_base + ".png")
    fig.savefig(out_base + ".pdf")
    print("wrote", out_base + ".png/.pdf")

    for trait, table in all_tables.items():
        table.to_csv(os.path.join(ROOT, "tables", f"Fig4_{trait}_top_bottom10.csv"))


if __name__ == "__main__":
    main()
