"""Recompute the 390-QC-kept-profile standardised score/p-value matrices, fixing a
standard-deviation convention bug found while re-verifying provenance.

Background: the existing `scores.xlsx`/`pvals.xlsx` (390 kept profiles x 1,125
traits, used throughout paper_analysis) were built by re-standardising the raw,
unstandardised CREST-GV Z-score within each collection using only its QC-kept
(Good/Decent) cell types -- correct in principle, but the standard deviation used
was the *population* std (ddof=0), while the package's own internal
standardisation (`crestgv.core`'s `__add_statistics`, which produced the original
`CREST-GV_norm` column against the full 511-profile set) uses pandas' default
*sample* std (ddof=1). Confirmed by reproducing a specific cell value exactly:
catlas_adult_A_Cardiomyocyte/heart_rate matches the saved 3.489502989041492 only
with ddof=0; ddof=1 gives 3.450513 instead. The gap scales with 1/sqrt(1-1/n), so
it is small for large collections (<2% for n>=19) but not negligible for the two
smallest (h1_hescs n=4: ~15.5%; ludwig2019 n=8: ~6.9%).

This script produces THREE outputs for comparison, all restricted to the 390 kept
profiles:
  1. scores_390kept_ddof1.xlsx / pvals_390kept_ddof1.xlsx
     The corrected re-standardisation: raw CREST-GV Z, restricted to kept cell
     types per collection, standardised with ddof=1 (matching the package's own
     convention). This is the candidate replacement for the current scores.xlsx.
  2. scores_390kept_from_norm.xlsx / pvals_390kept_from_norm.xlsx
     The original CREST-GV_norm/CREST-GV_pval_norm columns (standardised against
     the full, unfiltered 511-profile set), simply filtered down to the 390 kept
     rows -- no re-standardisation at all. For the retained rows this must be
     numerically identical to scores_unfiltered_511.xlsx (sanity check), and
     comparing it against output (1) isolates the effect of restricting the
     standardisation *basis* to kept-only, separately from the ddof fix.

Reads the raw per-trait, per-collection statistics CSVs directly (12,375 reads:
11 collections x 1,125 traits) rather than reusing the already-aggregated
unfiltered-511 file for output (1), since that file only stores the already-
standardised CREST-GV_norm, not the raw Z needed to re-derive a from-scratch
ddof=1 standardisation restricted to kept cell types.

Run from the repo root (paper_analysis/):
    python scripts/recompute_390kept_ddof1.py
"""

import os

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TRAITS_DIR = "/ceph/project/Wellcome_Discovery/sriva/Git/CREST-GV_v2_backup/ALL_GWAS_CATALOG_TRAITS"
QC_FILE = os.path.join(ROOT, "data", "CREST-GV_qc_summary_all.csv")
OUT_DIR = os.path.join(ROOT, "data")

COLLECTIONS = [
    "cad", "calderon", "catlas_fetal", "catlas_adult",
    "erythoid_d7_d10_d13_d17", "greenleaf_brain", "h1_hescs",
    "immune_cell", "ludwig2019", "mpal", "super_pbmc",
]


def main():
    qc = pd.read_csv(QC_FILE, sep="\t")
    kept = qc[qc["qc_label"].isin(["Good", "Decent"])]
    kept_by_collection = kept.groupby("collection")["cell_type"].apply(list).to_dict()

    old_unfiltered = pd.read_excel(
        os.path.join(TRAITS_DIR, "scores_unfiltered_511.xlsx"), index_col=0
    )
    traits = [c[: -len("_score")] for c in old_unfiltered.columns]
    print(f"{len(traits)} traits to process across {len(COLLECTIONS)} collections "
          f"({len(traits) * len(COLLECTIONS)} file reads)")

    all_kept_index = [f"{coll}_{ct}" for coll in COLLECTIONS for ct in kept_by_collection[coll]]

    ddof1_scores = pd.DataFrame(index=all_kept_index, columns=traits, dtype=float)
    ddof1_pvals = pd.DataFrame(index=all_kept_index, columns=traits, dtype=float)

    for ti, trait in enumerate(traits):
        if ti % 100 == 0:
            print(f"  trait {ti}/{len(traits)}: {trait}")
        for coll in COLLECTIONS:
            path = os.path.join(TRAITS_DIR, trait, f"{coll}_statistics_CREST-GV.csv")
            df = pd.read_csv(path, sep="\t", index_col=0)
            kept_cts = kept_by_collection[coll]
            sub = df.loc[df.index.isin(kept_cts), "CREST-GV"]
            mean, std = sub.mean(), sub.std(ddof=1)
            z = (sub - mean) / std
            p = stats.norm.sf(z)
            idx = [f"{coll}_{ct}" for ct in sub.index]
            ddof1_scores.loc[idx, trait] = z.values
            ddof1_pvals.loc[idx, trait] = p

    ddof1_scores.columns = [f"{t}_score" for t in traits]
    ddof1_pvals.columns = [f"{t}_score" for t in traits]

    ddof1_scores.to_excel(os.path.join(OUT_DIR, "scores_390kept_ddof1.xlsx"))
    ddof1_pvals.to_excel(os.path.join(OUT_DIR, "pvals_390kept_ddof1.xlsx"))
    ddof1_scores.to_csv(os.path.join(OUT_DIR, "scores_390kept_ddof1.tsv"), sep="\t")
    ddof1_pvals.to_csv(os.path.join(OUT_DIR, "pvals_390kept_ddof1.tsv"), sep="\t")
    print("wrote scores_390kept_ddof1.{xlsx,tsv} and pvals_390kept_ddof1.{xlsx,tsv}")

    # Output 2: original CREST-GV_norm, just filtered to kept rows, no recompute.
    old_pvals_unfiltered = pd.read_excel(
        os.path.join(TRAITS_DIR, "pvals_unfiltered_511.xlsx"), index_col=0
    )
    from_norm_scores = old_unfiltered.loc[old_unfiltered.index.isin(all_kept_index)]
    from_norm_pvals = old_pvals_unfiltered.loc[old_pvals_unfiltered.index.isin(all_kept_index)]
    from_norm_scores.to_excel(os.path.join(OUT_DIR, "scores_390kept_from_norm.xlsx"))
    from_norm_pvals.to_excel(os.path.join(OUT_DIR, "pvals_390kept_from_norm.xlsx"))
    print("wrote scores_390kept_from_norm.xlsx and pvals_390kept_from_norm.xlsx "
          f"({len(from_norm_scores)} rows, filtered from the unfiltered-511 file)")


if __name__ == "__main__":
    main()
