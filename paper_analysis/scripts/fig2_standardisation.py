"""Fig. 1c-d -- Cross-collection standardisation of enrichment scores.

(c) Raw Z_ij per collection, for the cell types retained after QC.
(d) Standardised Z-tilde_ij (same cell types), after within-collection
    standardisation.

Trait shown: red_blood_cell_traits (one of the four traits highlighted in the
Results; chosen here for consistency with Fig. 1e and Fig. 2, which use the
same four traits).

Run from the repo root:
    python scripts/fig2_standardisation.py
"""

import os
import sys

import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from palette import apply_base_style, SERIES_1, SERIES_2, INK_MUTED, GRIDLINE

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TRAITS_DIR = "/ceph/project/Wellcome_Discovery/sriva/Git/CREST-GV_v2_backup/ALL_GWAS_CATALOG_TRAITS"
QC_FILE = os.path.join(ROOT, "data", "CREST-GV_qc_summary_all.csv")
SCORES_FILE = os.path.join(ROOT, "data", "scores_390kept.tsv")

TRAIT = "red_blood_cell_traits"
TRAIT_COL = f"{TRAIT}_score"

# Collection display order (matches Table 1)
COLLECTIONS = [
    "cad", "catlas_fetal", "catlas_adult", "mpal", "super_pbmc", "calderon",
    "ludwig2019", "erythoid_d7_d10_d13_d17", "h1_hescs", "immune_cell", "greenleaf_brain",
]


def main():
    apply_base_style()

    qc = pd.read_csv(QC_FILE, sep="\t")
    kept = qc[qc["qc_label"].isin(["Good", "Decent"])]
    kept_by_collection = kept.groupby("collection")["cell_type"].apply(set).to_dict()

    scores = pd.read_csv(SCORES_FILE, sep="\t", index_col=0)

    raw_by_collection, norm_by_collection = [], []
    for coll in COLLECTIONS:
        stats = pd.read_csv(
            os.path.join(TRAITS_DIR, TRAIT, f"{coll}_statistics_CREST-GV.csv"),
            sep="\t", index_col=0,
        )
        kept_cts = kept_by_collection[coll]
        present = [ct for ct in stats.index if ct in kept_cts]
        raw_by_collection.append(stats.loc[present, "CREST-GV"].values)
        norm_keys = [f"{coll}_{ct}" for ct in present]
        norm_by_collection.append(scores.loc[norm_keys, TRAIT_COL].values)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    positions = range(1, len(COLLECTIONS) + 1)

    for ax, data, color, title in [
        (ax1, raw_by_collection, SERIES_1, "c   Raw enrichment score ($Z_{ij}$)"),
        (ax2, norm_by_collection, SERIES_2, r"d   Standardised score ($\tilde{Z}_{ij}$)"),
    ]:
        bp = ax.boxplot(
            data, positions=positions, widths=0.6, patch_artist=True, showfliers=True,
            medianprops=dict(color="#0b0b0b", linewidth=1.5),
            flierprops=dict(marker="o", markersize=3, markerfacecolor=color,
                             markeredgecolor="none", alpha=0.5),
            boxprops=dict(facecolor=color, edgecolor=color, alpha=0.35, linewidth=1.2),
            whiskerprops=dict(color=color, linewidth=1.2),
            capprops=dict(color=color, linewidth=1.2),
        )
        ax.axhline(0, color=GRIDLINE, linewidth=1, zorder=0)
        ax.set_title(title, loc="left", fontsize=11)
        ax.set_ylabel("Enrichment score")

    ax2.set_xticks(positions)
    ax2.set_xticklabels(COLLECTIONS, rotation=40, ha="right")
    fig.suptitle(f"Trait shown: {TRAIT.replace('_', ' ')}", fontsize=10, color=INK_MUTED, y=0.995)
    fig.tight_layout()

    out_base = os.path.join(ROOT, "figures", "Fig2_standardisation")
    fig.savefig(out_base + ".png")
    fig.savefig(out_base + ".pdf")
    print("wrote", out_base + ".png/.pdf")

    # Also dump the underlying numbers for the figure legend / caption text.
    summary = []
    for coll, raw, norm in zip(COLLECTIONS, raw_by_collection, norm_by_collection):
        summary.append({
            "collection": coll, "n_cell_types": len(raw),
            "raw_min": raw.min(), "raw_max": raw.max(),
            "norm_min": norm.min(), "norm_max": norm.max(),
        })
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(os.path.join(ROOT, "tables", "Fig2_ranges_by_collection.csv"), index=False)
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
