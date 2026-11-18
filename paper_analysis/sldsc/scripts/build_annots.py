#!/usr/bin/env python
"""Build per-cell-type S-LDSC .annot.gz files from CREST-GV's own LanceOtron peak calls.

One annotation per cell type (not per collection), one file per chromosome, in the
"thin annot" format ldsc expects (single ANNOT column, 0/1, in .bim SNP order) --
same convention as make_annot.py, but batched for speed:
  - each peak bed is sorted+merged once (not once per chromosome)
  - each chromosome's .bim is loaded once and reused across all ~730 cell types
Naive make_annot.py-per-celltype-per-chromosome would be 729 * 22 ~= 16,000 subprocess
calls (~3s each, mostly import/parse overhead) ~= 13h serial; this batches it to
minutes per chromosome, run across chromosomes in parallel.
"""
import argparse
import gzip
import multiprocessing as mp
import os
import sys

import pandas as pd
from pybedtools import BedTool

COLLECTION_ROOT = "/project/Wellcome_Discovery/datashare/CREST-GV_collection"
PLINK_DIR = "/project/Wellcome_Discovery/adelaney/Projects/02_sLDSC/00_ref_files/hg38/1000G_EUR_plink"
PLINK_PREFIX = "1000G.EUR.hg38"
CHROMS = [str(c) for c in range(1, 23)]


def discover_celltypes():
    """Return list of (collection, celltype, bed_path) for every *_L-tron.bed."""
    triples = []
    for collection in sorted(os.listdir(COLLECTION_ROOT)):
        peaks_dir = os.path.join(COLLECTION_ROOT, collection, "peaks")
        if not os.path.isdir(peaks_dir):
            continue
        for fname in sorted(os.listdir(peaks_dir)):
            if fname.endswith("_L-tron.bed"):
                celltype = fname[: -len("_L-tron.bed")]
                triples.append((collection, celltype, os.path.join(peaks_dir, fname)))
    return triples


def premerge_peaks(triples, cache_dir):
    """Sort+merge each peak bed once; cache as plain bed files. Returns dict
    (collection, celltype) -> cached merged bed path."""
    os.makedirs(cache_dir, exist_ok=True)
    cached = {}
    for collection, celltype, bed_path in triples:
        out_path = os.path.join(cache_dir, f"{collection}__{celltype}.merged.bed")
        if not os.path.exists(out_path):
            merged = BedTool(bed_path).sort().merge()
            merged.saveas(out_path)
        cached[(collection, celltype)] = out_path
    return cached


def build_chrom(chrom, triples, cached_beds, out_root):
    bim_path = os.path.join(PLINK_DIR, f"{PLINK_PREFIX}.{chrom}.bim")
    df_bim = pd.read_csv(
        bim_path, sep=r"\s+", usecols=[0, 1, 2, 3], names=["CHR", "SNP", "CM", "BP"]
    )
    iter_bim = [["chr" + str(c), bp - 1, bp] for c, bp in zip(df_bim["CHR"], df_bim["BP"])]
    bimbed = BedTool(iter_bim)

    n_done = 0
    for collection, celltype, _ in triples:
        out_dir = os.path.join(out_root, collection)
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, f"{celltype}.{chrom}.annot.gz")
        if os.path.exists(out_path):
            n_done += 1
            continue
        merged_bed_path = cached_beds[(collection, celltype)]
        annotbed = bimbed.intersect(BedTool(merged_bed_path))
        hit_bp = {x.start + 1 for x in annotbed}
        annot_col = [1 if bp in hit_bp else 0 for bp in df_bim["BP"]]
        with gzip.open(out_path, "wt") as f:
            f.write("ANNOT\n")
            f.write("\n".join(str(v) for v in annot_col))
            f.write("\n")
        n_done += 1
    return chrom, n_done


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-root", required=True)
    parser.add_argument("--cache-dir", required=True)
    parser.add_argument("--nproc", type=int, default=6)
    parser.add_argument("--chroms", nargs="*", default=CHROMS)
    args = parser.parse_args()

    triples = discover_celltypes()
    print(f"Found {len(triples)} cell types across "
          f"{len({c for c, _, _ in triples})} collections", file=sys.stderr)

    print("Pre-merging peak bed files...", file=sys.stderr)
    cached_beds = premerge_peaks(triples, args.cache_dir)

    os.makedirs(args.out_root, exist_ok=True)
    with mp.Pool(args.nproc) as pool:
        results = pool.starmap(
            build_chrom,
            [(c, triples, cached_beds, args.out_root) for c in args.chroms],
        )
    for chrom, n_done in sorted(results, key=lambda x: int(x[0])):
        print(f"chr{chrom}: {n_done}/{len(triples)} annot files present", file=sys.stderr)


if __name__ == "__main__":
    main()
