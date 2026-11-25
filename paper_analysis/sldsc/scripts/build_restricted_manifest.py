#!/usr/bin/env python
"""Restrict the full 729-cell-type --ref-ld-chr-cts manifest (cts.ldcts) to
the paper's actual atlas: 11 collections, 390 QC-kept profiles (manuscript
\\ncoll{11}/\\nkept{390}), rather than all 15 collections currently on disk.

Source of truth for which profiles are QC-kept:
paper_analysis/data/CREST-GV_qc_summary_all.csv (tab-separated; columns
collection, cell_type, prediction, qc_label). Kept = qc_label != "Bad".
That file uses crest_output-style lowercase collection names, which are
translated to this pipeline's directory-name convention (e.g. cad -> CAD,
h1_hescs -> H1_hESCs) via COLLECTION_TO_CRESTOUT (shared with
compare_with_crestgv.py; kept in sync manually -- small, stable mapping).

Collections not present in the QC file at all (FireCrest, Granulopoiesis,
pancreatic_pbmc, Turner_byDisease) are not part of the paper's atlas and
are dropped entirely.
"""
import argparse
import csv
from collections import Counter

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


def load_kept_celltypes(qc_csv_path):
    """Return {crestout_collection_name: {kept cell_type names}}."""
    kept = {}
    with open(qc_csv_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["qc_label"] != "Bad":
                kept.setdefault(row["collection"], set()).add(row["cell_type"])
    return kept


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--qc-csv", default="../data/CREST-GV_qc_summary_all.csv")
    parser.add_argument("--in-ldcts", default="combined_annot/cts.ldcts")
    parser.add_argument("--out-ldcts", default="combined_annot/cts_paper_atlas.ldcts")
    args = parser.parse_args()

    kept_by_crestout_collection = load_kept_celltypes(args.qc_csv)

    kept_lines = []
    n_dropped = 0
    with open(args.in_ldcts) as f:
        for line in f:
            name = line.split("\t", 1)[0]
            collection, celltype = name.split("__", 1)
            crestout_collection = COLLECTION_TO_CRESTOUT.get(collection)
            if crestout_collection is None:
                n_dropped += 1  # not part of the paper's 11-collection atlas
                continue
            if celltype not in kept_by_crestout_collection.get(crestout_collection, set()):
                n_dropped += 1  # QC-excluded profile
                continue
            kept_lines.append(line)

    print(f"kept: {len(kept_lines)}  dropped: {n_dropped}  total: {len(kept_lines) + n_dropped}")
    with open(args.out_ldcts, "w") as f:
        f.writelines(kept_lines)

    counts = Counter(line.split("\t", 1)[0].split("__", 1)[0] for line in kept_lines)
    for collection, n in sorted(counts.items()):
        print(f"{collection}: {n}")


if __name__ == "__main__":
    main()
