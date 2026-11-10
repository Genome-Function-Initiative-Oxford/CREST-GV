"""Supplementary Fig. -- QC score distributions.

Histogram of the chromatin-bigwig-qc prediction score across all 511 collected
tracks, with the Good/Decent/Bad class thresholds (Q1 = 4754.02, median =
6471.38 of the model's fitting-data distribution) marked and the three zones
shaded with the reserved status palette.

Run from the repo root:
    python scripts/suppfig_qc_distributions.py
"""

import os
import sys

import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from palette import apply_base_style, QC_LABEL_COLOR, INK_MUTED, SERIES_1

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
QC_FILE = os.path.join(ROOT, "data", "CREST-GV_qc_summary_all.csv")

Q1 = 4754.02
MEDIAN = 6471.38


def main():
    apply_base_style()

    qc = pd.read_csv(QC_FILE, sep="\t")
    qc = qc[qc["prediction"] != "NA"].copy()
    qc["prediction"] = qc["prediction"].astype(float)

    fig, ax = plt.subplots(figsize=(8, 5))

    xmin, xmax = qc["prediction"].min() - 200, qc["prediction"].max() + 200
    ax.axvspan(xmin, Q1, color=QC_LABEL_COLOR["Bad"], alpha=0.08, zorder=0)
    ax.axvspan(Q1, MEDIAN, color=QC_LABEL_COLOR["Decent"], alpha=0.10, zorder=0)
    ax.axvspan(MEDIAN, xmax, color=QC_LABEL_COLOR["Good"], alpha=0.08, zorder=0)

    ax.hist(qc["prediction"], bins=40, color=SERIES_1, edgecolor="white", linewidth=0.5)

    ax.axvline(Q1, color=QC_LABEL_COLOR["Bad"], linewidth=1.5, linestyle="--")
    ax.axvline(MEDIAN, color=QC_LABEL_COLOR["Good"], linewidth=1.5, linestyle="--")

    ymax = ax.get_ylim()[1]
    ax.text(xmin + (Q1 - xmin) / 2, ymax * 0.95, "Bad", ha="center", color=QC_LABEL_COLOR["Bad"], fontsize=10, weight="bold")
    ax.text(Q1 + (MEDIAN - Q1) / 2, ymax * 0.95, "Decent", ha="center", color="#8a6a00", fontsize=10, weight="bold")
    ax.text(MEDIAN + (xmax - MEDIAN) / 2, ymax * 0.95, "Good", ha="center", color=QC_LABEL_COLOR["Good"], fontsize=10, weight="bold")

    ax.set_xlabel("Predicted constant promoter + CTCF site recovery score")
    ax.set_ylabel("Number of tracks")
    ax.set_title(
        f"QC score distribution, n={len(qc)} tracks "
        f"(Q1={Q1:.0f}, median={MEDIAN:.0f})",
        loc="left", fontsize=11,
    )
    fig.tight_layout()

    out_base = os.path.join(ROOT, "figures", "SuppFig_QC_distributions")
    fig.savefig(out_base + ".png")
    fig.savefig(out_base + ".pdf")
    print("wrote", out_base + ".png/.pdf")

    counts = qc["qc_label"].value_counts()
    print(counts)
    counts.to_csv(os.path.join(ROOT, "tables", "SuppFig_QC_label_counts.csv"))


if __name__ == "__main__":
    main()
