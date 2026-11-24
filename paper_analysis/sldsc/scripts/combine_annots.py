#!/usr/bin/env python
"""Combine all 729 per-cell-type single-column .annot.gz files (one chromosome at a
time) into one wide thin-annot file with 729 columns.

Computing LD scores (ldsc.py --l2) one annotation at a time would repeat the
expensive genotype window-scan 729 times per chromosome. ldsc computes LD scores
for every column in an annot file in a single pass, so combining columns first
and calling --l2 once per chromosome does the same work ~729x fewer times.
"""
import argparse
import gzip
import os

import pandas as pd

ANNOT_ROOT_DEFAULT = "annot"


def discover_celltypes(annot_root):
    triples = []
    for collection in sorted(os.listdir(annot_root)):
        coll_dir = os.path.join(annot_root, collection)
        if not os.path.isdir(coll_dir):
            continue
        celltypes = sorted(
            {f.split(".")[0] for f in os.listdir(coll_dir) if f.endswith(".annot.gz")}
        )
        for celltype in celltypes:
            triples.append((collection, celltype))
    return triples


def combine_chrom(chrom, triples, annot_root, out_path):
    cols = {}
    for collection, celltype in triples:
        path = os.path.join(annot_root, collection, f"{celltype}.{chrom}.annot.gz")
        s = pd.read_csv(path, compression="gzip")["ANNOT"]
        cols[f"{collection}__{celltype}"] = s.values
    df = pd.DataFrame(cols)
    with gzip.open(out_path, "wt") as f:
        df.to_csv(f, sep="\t", index=False)
    return chrom, df.shape


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--annot-root", default=ANNOT_ROOT_DEFAULT)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--chroms", nargs="+", required=True)
    args = parser.parse_args()

    triples = discover_celltypes(args.annot_root)
    print(f"Combining {len(triples)} cell-type annotations", flush=True)

    os.makedirs(args.out_dir, exist_ok=True)
    for chrom in args.chroms:
        out_path = os.path.join(args.out_dir, f"combined.{chrom}.annot.gz")
        if os.path.exists(out_path):
            print(f"chr{chrom}: already exists, skipping", flush=True)
            continue
        chrom, shape = combine_chrom(chrom, triples, args.annot_root, out_path)
        print(f"chr{chrom}: wrote {shape}", flush=True)


if __name__ == "__main__":
    main()
