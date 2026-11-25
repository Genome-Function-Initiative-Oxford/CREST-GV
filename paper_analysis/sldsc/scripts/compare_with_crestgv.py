#!/usr/bin/env python
"""Compare S-LDSC --h2-cts cell-type-specific coefficients (restricted to the
paper's 11-collection, 390-profile QC-kept atlas) against CREST-GV's own
standardised per-cell-type enrichment scores, for the same 4 traits.

Joins paper_analysis/sldsc/h2_cts_paper_atlas/<trait>.cell_type_results.txt
against local_data/crest_output/<trait>/<collection>_statistics_CREST-GV.csv
(CREST-GV_norm column), matching by collection+cell type name after
translating between the two naming conventions, and reports Spearman rank
correlation per trait.
"""
import argparse
import os

import pandas as pd
from scipy.stats import spearmanr

CREST_GV_ROOT_DEFAULT = "../../local_data/crest_output"
H2_CTS_DIR_DEFAULT = "h2_cts_paper_atlas"

# h2_cts collection name -> crest_output collection file prefix
COLLECTION_TO_CRESTOUT = {
    "CAD": "cad",
    "calderon": "calderon",
    "catlas_adult": "catlas_adult",
    "catlas_fetal": "catlas_fetal",
    "Days7_10_13_17": "erythoid_d7_d10_d13_d17",
    "greenleaf_brain": "greenleaf_brain",
    "H1_hESCs": "h1_hescs",
    "immune_cell": "immune_cell",
    "ludwig2019": "ludwig2019",
    "MPAL": "mpal",
    "super_PBMC": "super_pbmc",
}

TRAITS = ["heart_rate", "type_1_diabetes_mellitus", "type_2_diabetes_mellitus", "red_blood_cell_traits"]


def load_crestgv_scores(trait, crestgv_root):
    rows = {}
    trait_dir = os.path.join(crestgv_root, trait)
    for h2cts_coll, crestout_coll in COLLECTION_TO_CRESTOUT.items():
        path = os.path.join(trait_dir, f"{crestout_coll}_statistics_CREST-GV.csv")
        df = pd.read_csv(path, sep="\t")
        for _, row in df.iterrows():
            name = f"{h2cts_coll}__{row['cellType']}"
            rows[name] = row["CREST-GV_norm"]
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--crestgv-root", default=CREST_GV_ROOT_DEFAULT)
    parser.add_argument("--h2-cts-dir", default=H2_CTS_DIR_DEFAULT)
    args = parser.parse_args()

    print(f"{'trait':<26}{'n':>5}{'rho':>10}{'p':>14}")
    for trait in TRAITS:
        sldsc = pd.read_csv(os.path.join(args.h2_cts_dir, f"{trait}.cell_type_results.txt"), sep="\t")
        crestgv = load_crestgv_scores(trait, args.crestgv_root)
        sldsc["CREST_GV_norm"] = sldsc["Name"].map(crestgv)
        missing = sldsc["CREST_GV_norm"].isna().sum()
        if missing:
            print(f"WARNING: {trait}: {missing} cell types missing a CREST-GV match")
        merged = sldsc.dropna(subset=["CREST_GV_norm"])
        rho, p = spearmanr(merged["Coefficient"], merged["CREST_GV_norm"])
        print(f"{trait:<26}{len(merged):>5}{rho:>10.3f}{p:>14.2e}")
        out_path = os.path.join(args.h2_cts_dir, f"{trait}.vs_crestgv.tsv")
        merged.sort_values("Coefficient_P_value").to_csv(out_path, sep="\t", index=False)


if __name__ == "__main__":
    main()
