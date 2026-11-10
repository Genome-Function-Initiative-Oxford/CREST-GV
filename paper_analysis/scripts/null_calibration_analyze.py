"""Analyze null calibration: is CREST-GV's P-value well-calibrated?

For each of the 4 highlighted traits, 100 random genome-wide null variant
sets (matched in size to the real trait) were run through CREST-GV against
all 11 atlas collections (null_calibration_run.py). Under a well-calibrated
test, P-values should be uniformly distributed on [0,1], so the empirical
false-positive rate at a nominal threshold alpha should equal alpha, and a
QQ-plot of observed vs expected -log10(P) should follow the y=x line.

Calibrates the STANDARDISED P-value (the one the manuscript's P<\\pthr
significance criterion actually uses), not the raw per-collection P-value
null_calibration_run.py recorded. Re-derives it offline from the already-
collected raw Z: for each (trait, replicate, collection), restandardises
using only the QC-kept cell types in that collection (mean/std computed over
just those), exactly matching how the real atlas was re-standardised -- no
need to rerun CREST-GV itself.

Restricted throughout to the 390 QC-kept cell types (the published atlas).

Run from the repo root (after null_calibration_run.py):
    python scripts/null_calibration_analyze.py
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
from palette import apply_base_style, CATEGORICAL, GRIDLINE, INK_MUTED, SERIES_1

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NULL_DIR = os.path.join(ROOT, "null_calibration")
QC_FILE = os.path.join(ROOT, "data", "CREST-GV_qc_summary_all.csv")

TRAITS = ["heart_rate", "red_blood_cell_traits", "type_1_diabetes_mellitus", "type_2_diabetes_mellitus"]
TRAIT_LABELS = {
    "heart_rate": "Heart rate", "red_blood_cell_traits": "Red blood cell traits",
    "type_1_diabetes_mellitus": "Type 1 diabetes", "type_2_diabetes_mellitus": "Type 2 diabetes",
}
TRAIT_COLORS = {
    "heart_rate": CATEGORICAL["blue"], "red_blood_cell_traits": CATEGORICAL["red"],
    "type_1_diabetes_mellitus": CATEGORICAL["aqua"], "type_2_diabetes_mellitus": CATEGORICAL["orange"],
}
ALPHAS = [0.05, 0.01]


def main():
    apply_base_style()

    qc = pd.read_csv(QC_FILE, sep="\t")
    kept = qc[qc["qc_label"].isin(["Good", "Decent"])]
    kept_keys = set(kept["collection"] + "_" + kept["cell_type"])

    all_data = {}
    fpr_rows = []
    per_collection_rows = []
    at_least_one_rows = []

    for trait in TRAITS:
        df = pd.read_csv(os.path.join(NULL_DIR, f"{trait}_null_results.tsv"), sep="\t")
        df["key"] = df["collection"] + "_" + df["cellType"]
        df = df[df["key"].isin(kept_keys)].copy()

        # Re-standardise within (replicate, collection), using only the kept
        # cell types present -- same recipe as the real 390-profile atlas.
        grp = df.groupby(["replicate", "collection"])["Z"]
        mean_z = grp.transform("mean")
        std_z = grp.transform(lambda s: s.std(ddof=0))  # match np.std default (ddof=0), as in crestgv core
        z_tilde = (df["Z"] - mean_z) / std_z
        df["pval"] = stats.norm.sf(z_tilde)
        df["pval"] = df["pval"].clip(lower=1e-300)  # guard log10(0)
        all_data[trait] = df

        n_tests = len(df)
        for alpha in ALPHAS:
            fpr = (df["pval"] < alpha).mean()
            fpr_rows.append({"trait": trait, "alpha": alpha, "empirical_fpr": fpr, "n_tests": n_tests})

        for coll, sub in df.groupby("collection"):
            for alpha in ALPHAS:
                fpr = (sub["pval"] < alpha).mean()
                per_collection_rows.append({
                    "trait": trait, "collection": coll, "alpha": alpha,
                    "empirical_fpr": fpr, "n_tests": len(sub),
                })

        for alpha in ALPHAS:
            per_rep_any_hit = df.groupby("replicate")["pval"].apply(lambda p: (p < alpha).any())
            at_least_one_rows.append({
                "trait": trait, "alpha": alpha,
                "p_at_least_one_false_hit": per_rep_any_hit.mean(),
                "n_cell_types_tested": df.groupby("replicate").size().iloc[0],
            })

    fpr_df = pd.DataFrame(fpr_rows)
    per_collection_df = pd.DataFrame(per_collection_rows)
    at_least_one_df = pd.DataFrame(at_least_one_rows)

    os.makedirs(os.path.join(ROOT, "tables"), exist_ok=True)
    fpr_df.to_csv(os.path.join(ROOT, "tables", "NullCalibration_fpr_by_trait.csv"), index=False)
    per_collection_df.to_csv(os.path.join(ROOT, "tables", "NullCalibration_fpr_by_collection.csv"), index=False)
    at_least_one_df.to_csv(os.path.join(ROOT, "tables", "NullCalibration_atleastone_hit.csv"), index=False)

    print("=== Empirical false-positive rate vs nominal alpha ===")
    print(fpr_df.pivot(index="trait", columns="alpha", values="empirical_fpr").to_string())
    print("\n=== P(>=1 false hit across the 390-profile atlas in one null replicate) ===")
    print(at_least_one_df.pivot(index="trait", columns="alpha", values="p_at_least_one_false_hit").to_string())

    # --- Figure: (a) QQ-plot pooling all traits, (b) empirical FPR vs nominal alpha per trait ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

    for trait, color in TRAIT_COLORS.items():
        pvals = np.sort(all_data[trait]["pval"].values)
        n = len(pvals)
        expected = -np.log10(np.arange(1, n + 1) / (n + 1))
        observed = -np.log10(pvals)
        # subsample for plotting speed/clarity (51,100 points per trait)
        idx = np.linspace(0, n - 1, 3000).astype(int)
        ax1.scatter(expected[idx], observed[idx], s=4, alpha=0.4, color=color, label=TRAIT_LABELS[trait])

    lims = [0, max(ax1.get_xlim()[1], ax1.get_ylim()[1])]
    ax1.plot(lims, lims, color=GRIDLINE, linewidth=1.5, linestyle="--", zorder=0)
    ax1.set_xlabel("Expected -log10(P) under a uniform null")
    ax1.set_ylabel("Observed -log10(P)")
    ax1.set_title("a   QQ-plot: null-replicate P-values vs uniform expectation", loc="left", fontsize=10)
    ax1.legend(frameon=False, fontsize=8, loc="upper left")

    width = 0.35
    x = np.arange(len(TRAITS))
    for i, alpha in enumerate(ALPHAS):
        vals = [fpr_df[(fpr_df.trait == t) & (fpr_df.alpha == alpha)]["empirical_fpr"].iloc[0] for t in TRAITS]
        ax2.bar(x + (i - 0.5) * width, vals, width=width, label=f"nominal P<{alpha}",
                color=SERIES_1 if i == 0 else CATEGORICAL["orange"])
        ax2.axhline(alpha, color=(SERIES_1 if i == 0 else CATEGORICAL["orange"]), linewidth=1.2, linestyle=":")

    ax2.set_xticks(x)
    ax2.set_xticklabels([TRAIT_LABELS[t] for t in TRAITS], rotation=20, ha="right")
    ax2.set_ylabel("Empirical false-positive rate")
    ax2.set_title("b   Empirical FPR vs nominal alpha (dotted = expected)", loc="left", fontsize=10)
    ax2.legend(frameon=False, fontsize=8)

    fig.tight_layout()
    out_base = os.path.join(ROOT, "figures", "SuppFig_null_calibration")
    fig.savefig(out_base + ".png")
    fig.savefig(out_base + ".pdf")
    print(f"\nwrote {out_base}.png/.pdf")


if __name__ == "__main__":
    main()
