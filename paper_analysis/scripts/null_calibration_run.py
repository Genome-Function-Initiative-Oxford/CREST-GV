"""Run CREST-GV on every null variant replicate, across all 11 atlas collections.

Uses the bug-fixed CREST-GV_newCode package (division-by-zero guards matter
here -- a random null set is exactly the kind of input that could plausibly
hit a degenerate p_succes or zero-variance collection). For each (trait,
replicate, collection) triple, records raw Z and P for every cell type.
Results are appended per-replicate to one TSV per trait (resumable: reps
already present in the output file are skipped on restart).

Run from the repo root:
    python scripts/null_calibration_run.py
"""

import glob
import os
import sys
import time

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NEWCODE_SRC = os.path.join(os.path.dirname(ROOT), "CREST-GV_newCode", "src")
sys.path.insert(0, NEWCODE_SRC)
from crestgv import crestgv  # noqa: E402

NULL_DIR = os.path.join(ROOT, "data", "null_calibration")
OUT_DIR = os.path.join(ROOT, "null_calibration")
os.makedirs(OUT_DIR, exist_ok=True)

TRAITS = ["heart_rate", "red_blood_cell_traits", "type_1_diabetes_mellitus", "type_2_diabetes_mellitus"]
COLLECTIONS = [
    "cad", "catlas_fetal", "catlas_adult", "mpal", "super_pbmc", "calderon",
    "ludwig2019", "erythoid_d7_d10_d13_d17", "h1_hescs", "immune_cell", "greenleaf_brain",
]


def completed_reps(out_path):
    if not os.path.exists(out_path):
        return set()
    return set(pd.read_csv(out_path, sep="\t", usecols=["replicate"])["replicate"].unique())


def main():
    for trait in TRAITS:
        out_path = os.path.join(OUT_DIR, f"{trait}_null_results.tsv")
        done = completed_reps(out_path)

        null_files = sorted(glob.glob(os.path.join(NULL_DIR, trait, "null_rep*.bed")))
        print(f"[{trait}] {len(null_files)} replicates, {len(done)} already done", flush=True)

        t_trait0 = time.time()
        for null_path in null_files:
            rep = int(os.path.basename(null_path).replace("null_rep", "").replace(".bed", ""))
            if rep in done:
                continue

            t0 = time.time()
            rows = []
            for coll in COLLECTIONS:
                cgv = crestgv(genetic=null_path, output="/tmp/null_calibration_scratch", collection_name=coll)
                df = cgv.calculate_enrichment_score()
                for ct, row in df.iterrows():
                    rows.append({
                        "replicate": rep, "collection": coll, "cellType": ct,
                        "Z": row["CREST-GV"], "pval": row["CREST-GV_pval"],
                    })

            rep_df = pd.DataFrame(rows)
            write_header = not os.path.exists(out_path)
            rep_df.to_csv(out_path, sep="\t", mode="a", header=write_header, index=False)

            elapsed = time.time() - t0
            print(f"[{trait}] rep {rep} done in {elapsed:.1f}s ({len(rows)} rows)", flush=True)

        print(f"[{trait}] all replicates done, total {time.time()-t_trait0:.1f}s", flush=True)

    print("DONE.", flush=True)


if __name__ == "__main__":
    main()
