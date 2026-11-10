"""Regenerate benchmarks/cheers_comparison/docs/scatter_data.json and
cheers_vs_crestgv_scatter.html in the CREST-GV GitHub repo from the current
(ddof=1-corrected, 11-collection, pancreatic_pbmc-excluded) CHEERS benchmark data.

Computes real CHEERS Z-scores from the raw observedMeanRank + log files (mean/sd
of the null rank distribution), rather than relying on CHEERS_pval, which
underflows to exactly 0.0 for a large fraction of cell types and is therefore
useless for picking out a true top hit by magnitude.

Run from paper_analysis/:
    python scripts/regenerate_cheers_scatter_docs.py
"""

import glob
import json
import os

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ENRICH_DIR = "/ceph/project/Wellcome_Discovery/sriva/Git/CHEERS/CREST-GV_vs_CHEERS/results/enrichment"
DOCS_DIR = "/project/Wellcome_Discovery/sriva/Git/CREST-GV/benchmarks/cheers_comparison/docs"

COLLECTIONS = {
    "cad": "CAD", "calderon": "calderon", "catlas_fetal": "catlas_fetal",
    "catlas_adult": "catlas_adult", "erythoid_d7_d10_d13_d17": "Days7_10_13_17",
    "h1_hescs": "H1_hESCs", "immune_cell": "immune_cell", "ludwig2019": "ludwig2019",
    "mpal": "MPAL", "super_pbmc": "super_PBMC", "greenleaf_brain": "greenleaf_brain",
}
TRAITS = ["heart_rate", "red_blood_cell_traits", "type_1_diabetes_mellitus", "type_2_diabetes_mellitus"]


def cheers_z_for(coll_dir, trait, cell_types):
    base = os.path.join(ENRICH_DIR, coll_dir, trait)
    params = {}
    for line in open(os.path.join(base, f"{trait}.log")):
        k, v = line.strip().split("\t")
        params[k] = float(v)
    mm, sd = params["Distribution mean"], params["Distribution sd"]
    rank = pd.read_csv(
        os.path.join(base, f"{trait}_disease_enrichment_observedMeanRank.txt"),
        sep="\t", header=None, names=["cellType", "obs"],
    )
    rank["cellType"] = rank["cellType"].str.replace("_ReadsInPeaks", "", regex=False)
    rank["Z"] = (rank["obs"] - mm) / sd
    return rank.set_index("cellType")["Z"].reindex(cell_types)


def main():
    data = {}
    for coll, label in COLLECTIONS.items():
        traits_out = {}
        for trait in TRAITS:
            path = os.path.join(ROOT, "cheers_benchmark", "per_pair", f"{coll}__{trait}.tsv")
            df = pd.read_csv(path, sep="\t")
            cheers_z = cheers_z_for(label, trait, df["cellType"])
            points = []
            for _, row in df.iterrows():
                points.append({
                    "cellType": row["cellType"],
                    "cheersZ": round(float(cheers_z.loc[row["cellType"]]), 4),
                    "crestZ": round(float(row["CREST-GV_norm"]), 4),
                    "cheersP": float(row["CHEERS_pval"]),
                    "crestPNorm": round(float(row["CREST-GV_pval_norm"]), 6),
                })
            from scipy.stats import spearmanr
            rho, _ = spearmanr(df["CHEERS_pval"], df["CREST-GV_pval_norm"])
            traits_out[trait] = {"points": points, "rho": round(float(rho), 3), "n": len(df)}
        data[coll] = {"label": label, "traits": traits_out}

    with open(os.path.join(DOCS_DIR, "scatter_data.json"), "w") as f:
        json.dump(data, f)
    print("wrote", os.path.join(DOCS_DIR, "scatter_data.json"))

    html_path = os.path.join(DOCS_DIR, "cheers_vs_crestgv_scatter.html")
    with open(html_path) as f:
        lines = f.readlines()
    new_line_idx = next(i for i, l in enumerate(lines) if l.startswith("const DATA = "))
    lines[new_line_idx] = f"const DATA = {json.dumps(data)};\n"
    with open(html_path, "w") as f:
        f.writelines(lines)
    print("patched", html_path)


if __name__ == "__main__":
    main()
