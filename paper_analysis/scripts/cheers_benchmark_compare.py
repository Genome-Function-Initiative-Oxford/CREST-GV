"""Re-run the CREST-GV vs CHEERS benchmark comparison on the final atlas.

This does NOT re-run CHEERS itself (consensus peaks / bigwig quantification /
CHEERS_normalize.py / CHEERS_computeEnrichment.py): those steps are unaffected
by the catlas_fetal QC-pipeline pairing bug (they read the correct, complete
`peaks/` folder) and already cover all 11 atlas collections — verified by
diffing a sample pair's raw CREST-GV stats byte-for-byte against the canonical
ALL_GWAS_CATALOG_TRAITS copy (identical). Re-executing the bedtools /
bigWigAverageOverBed steps for ~500 tracks would be pure recompute for no
different answer.

What *does* change relative to the original comparison
(`CHEERS/CREST-GV_vs_CHEERS/scripts/compare_with_crestgv.py`):
  1. `pancreatic_pbmc` dropped (11 collections, not 12 -- it isn't part of the
     published atlas).
  2. Restricted to QC-kept (Good/Decent) cell types only, per collection --
     the original script compared against every cell type including Bad ones.
  3. CREST-GV's standardised score/p-value now comes from the final
     390-profile, QC-filtered re-standardisation
     (paper_analysis/data/{scores,pvals}_390kept.tsv), not the old
     full-collection (including Bad tracks) standardisation the CHEERS
     comparison previously used.
  4. Overall / per-trait / per-collection medians reported explicitly (the
     todo item asked for these).
  5. Significance-overlap counts now use Benjamini-Hochberg FDR q<0.10,
     pooled across all 390 kept cell types per trait -- the same family and
     threshold the manuscript's own significance criterion uses (not the
     raw, uncorrected P<0.05 the original comparison used). Applied to BOTH
     methods: CHEERS's own enrichment P-value has no correction of its own
     either (checked CHEERS_computeEnrichment.py directly -- no FDR/
     Bonferroni step anywhere in it), so comparing uncorrected CHEERS against
     corrected CREST-GV would not be a fair "do both methods agree on the
     few formally significant hits" comparison. Spearman rho itself needs no
     rerun: BH-FDR q is a monotonic transform of the raw P-value, so rank
     correlation is identical whichever scale is used.

Run from the repo root:
    python scripts/cheers_benchmark_compare.py
"""

import os
import sys

import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stats_utils import bh_fdr_mask

FDR_Q = 0.10

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CHEERS_PIPELINE = "/ceph/project/Wellcome_Discovery/sriva/Git/CHEERS/CREST-GV_vs_CHEERS"
TRAITS_DIR = "/ceph/project/Wellcome_Discovery/sriva/Git/CREST-GV_v2_backup/ALL_GWAS_CATALOG_TRAITS"
QC_FILE = os.path.join(ROOT, "data", "CREST-GV_qc_summary_all.csv")
SCORES_FILE = os.path.join(ROOT, "data", "scores_390kept.tsv")
PVALS_FILE = os.path.join(ROOT, "data", "pvals_390kept.tsv")
OUT_DIR = os.path.join(ROOT, "cheers_benchmark")
PER_PAIR_DIR = os.path.join(OUT_DIR, "per_pair")

# collection_key -> on-disk collection directory name (pancreatic_pbmc excluded)
COLLECTIONS = {
    "cad": "CAD",
    "calderon": "calderon",
    "catlas_fetal": "catlas_fetal",
    "catlas_adult": "catlas_adult",
    "erythoid_d7_d10_d13_d17": "Days7_10_13_17",
    "h1_hescs": "H1_hESCs",
    "immune_cell": "immune_cell",
    "ludwig2019": "ludwig2019",
    "mpal": "MPAL",
    "super_pbmc": "super_PBMC",
    "greenleaf_brain": "greenleaf_brain",
}

TRAITS = ["heart_rate", "red_blood_cell_traits", "type_1_diabetes_mellitus", "type_2_diabetes_mellitus"]


def load_cheers_pvalues(collection_dir, trait_name):
    path = os.path.join(
        CHEERS_PIPELINE, "results", "enrichment", collection_dir, trait_name,
        f"{trait_name}_disease_enrichment_pValues.txt",
    )
    df = pd.read_csv(path, sep="\t", header=None, names=["cellType", "CHEERS_pval"])
    df["cellType"] = df["cellType"].str.replace("_ReadsInPeaks", "", regex=False)
    return df


def load_crestgv_raw(collection_key, trait_name):
    path = os.path.join(TRAITS_DIR, trait_name, f"{collection_key}_statistics_CREST-GV.csv")
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={df.columns[0]: "cellType"})
    return df[["cellType", "CREST-GV", "CREST-GV_pval"]]


def compare_pair(collection_key, trait_name, kept_cell_types, norm_scores, norm_pvals):
    collection_dir = COLLECTIONS[collection_key]
    col = f"{trait_name}_score"

    cheers_df = load_cheers_pvalues(collection_dir, trait_name)
    crestgv_raw = load_crestgv_raw(collection_key, trait_name)

    norm_df = pd.DataFrame({
        "key": [f"{collection_key}_{ct}" for ct in kept_cell_types],
        "cellType": list(kept_cell_types),
    })
    norm_df["CREST-GV_norm"] = norm_df["key"].map(norm_scores[col])
    norm_df["CREST-GV_pval_norm"] = norm_df["key"].map(norm_pvals[col])
    norm_df = norm_df.drop(columns="key")

    merged = cheers_df.merge(crestgv_raw, on="cellType", how="inner").merge(norm_df, on="cellType", how="inner")
    # restrict to QC-kept cell types only
    merged = merged[merged["cellType"].isin(kept_cell_types)]

    merged.insert(0, "collection", collection_key)
    merged.insert(1, "trait", trait_name)

    summary = {"collection": collection_key, "trait": trait_name, "n_cell_types": len(merged)}
    if len(merged) >= 3:
        rho_raw, _ = spearmanr(merged["CHEERS_pval"], merged["CREST-GV_pval"])
        rho_norm, _ = spearmanr(merged["CHEERS_pval"], merged["CREST-GV_pval_norm"])
        summary["spearman_rho_vs_raw"] = rho_raw
        summary["spearman_rho_vs_norm"] = rho_norm
    else:
        summary["spearman_rho_vs_raw"] = float("nan")
        summary["spearman_rho_vs_norm"] = float("nan")

    return merged, summary


def main():
    os.makedirs(PER_PAIR_DIR, exist_ok=True)

    qc = pd.read_csv(QC_FILE, sep="\t")
    kept = qc[qc["qc_label"].isin(["Good", "Decent"])]
    kept_by_collection = kept.groupby("collection")["cell_type"].apply(set).to_dict()

    norm_scores = pd.read_csv(SCORES_FILE, sep="\t", index_col=0)
    norm_pvals = pd.read_csv(PVALS_FILE, sep="\t", index_col=0)

    summaries = []
    merged_by_trait = {t: [] for t in TRAITS}
    for collection_key in COLLECTIONS:
        for trait_name in TRAITS:
            kept_cts = kept_by_collection[collection_key]
            try:
                merged, summary = compare_pair(collection_key, trait_name, kept_cts, norm_scores, norm_pvals)
            except FileNotFoundError as e:
                print(f"[{collection_key}/{trait_name}] SKIPPED: {e}", file=sys.stderr)
                continue

            summaries.append(summary)
            merged_by_trait[trait_name].append(merged)
            print(
                f"[{collection_key}/{trait_name}] n={summary['n_cell_types']} "
                f"rho(raw)={summary['spearman_rho_vs_raw']:.3f} rho(norm)={summary['spearman_rho_vs_norm']:.3f}"
            )

    # --- BH-FDR significance, pooled across all 11 collections per trait (390 kept
    # cell types) -- the same family/threshold the manuscript's own significance
    # criterion uses, applied to BOTH methods (CHEERS has no correction of its own).
    fdr_overlap_rows = []
    for trait_name in TRAITS:
        pooled = pd.concat(merged_by_trait[trait_name], ignore_index=True)
        pooled["key"] = pooled["collection"] + "_" + pooled["cellType"]
        pooled = pooled.set_index("key")

        cheers_sig = bh_fdr_mask(pooled["CHEERS_pval"], FDR_Q)
        crestgv_sig = bh_fdr_mask(pooled["CREST-GV_pval_norm"], FDR_Q)
        pooled["CHEERS_fdr_significant"] = cheers_sig
        pooled["CREST-GV_fdr_significant"] = crestgv_sig

        fdr_overlap_rows.append({
            "trait": trait_name,
            "n_cell_types": len(pooled),
            "n_significant_CHEERS": int(cheers_sig.sum()),
            "n_significant_CREST-GV": int(crestgv_sig.sum()),
            "n_significant_both": int((cheers_sig & crestgv_sig).sum()),
            "n_significant_either": int((cheers_sig | crestgv_sig).sum()),
        })

        # write per-pair files (now including the FDR flag columns) and
        # per-collection FDR counts back into the main summary
        for collection_key, grp in pooled.reset_index(drop=True).groupby("collection"):
            grp.drop(columns=["collection", "trait"]).sort_values("CHEERS_pval").to_csv(
                os.path.join(PER_PAIR_DIR, f"{collection_key}__{trait_name}.tsv"), sep="\t", index=False
            )
            for s in summaries:
                if s["collection"] == collection_key and s["trait"] == trait_name:
                    s["n_significant_CHEERS_fdr"] = int(grp["CHEERS_fdr_significant"].sum())
                    s["n_significant_CREST-GV_fdr"] = int(grp["CREST-GV_fdr_significant"].sum())
                    s["n_significant_both_fdr"] = int(
                        (grp["CHEERS_fdr_significant"] & grp["CREST-GV_fdr_significant"]).sum()
                    )

    summary_df = pd.DataFrame(summaries)
    summary_df.to_csv(os.path.join(OUT_DIR, "summary.tsv"), sep="\t", index=False)

    fdr_overlap_df = pd.DataFrame(fdr_overlap_rows)
    fdr_overlap_df.to_csv(os.path.join(OUT_DIR, "fdr_significance_overlap.tsv"), sep="\t", index=False)

    overall_median = summary_df["spearman_rho_vs_norm"].median()
    per_trait = summary_df.groupby("trait")["spearman_rho_vs_norm"].median().sort_values(ascending=False)
    per_collection = summary_df.groupby("collection")["spearman_rho_vs_norm"].median().sort_values(ascending=False)

    with open(os.path.join(OUT_DIR, "medians.tsv"), "w") as f:
        f.write(f"overall_median_spearman_rho\t{overall_median:.4f}\n\n")
        f.write("per_trait_median_spearman_rho\n")
        per_trait.to_csv(f, sep="\t", header=["median_spearman_rho"])
        f.write("\nper_collection_median_spearman_rho\n")
        per_collection.to_csv(f, sep="\t", header=["median_spearman_rho"])

    print(f"\n=== {len(summary_df)} collection x trait pairs ===")
    print(f"Overall median Spearman rho (vs standardised CREST-GV p-value): {overall_median:.3f}")
    print("\nPer-trait median:")
    print(per_trait.to_string())
    print("\nPer-collection median:")
    print(per_collection.to_string())
    print(f"\n=== FDR (q<{FDR_Q}) significance overlap, pooled across all 11 collections per trait ===")
    print(fdr_overlap_df.to_string(index=False))
    print(f"\nWrote: {OUT_DIR}/summary.tsv, {OUT_DIR}/medians.tsv, {OUT_DIR}/fdr_significance_overlap.tsv, {PER_PAIR_DIR}/*.tsv")


if __name__ == "__main__":
    main()
